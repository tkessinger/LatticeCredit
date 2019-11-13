#!/usr/bin/env julia

## LatticeCredit.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating central bank lending and repayment.

module LatticeLending

    using Random

    export LatticeGame, LatticePopulation
    export evolve_sync!, evolve_async!, randomize_lattice!
    export decompose_strategy, strategy_string
    export get_payoffs, get_payoffs_fast, get_group_payoffs
    export pair_PGG_payoffs
	export update_reputations!

    struct LatticeGame
        # new type for "game" parameters
        # this should minimize the amount of passing around that needs to be done

        r::Float64 # synergy parameter
        z::Float64 # interest on loans (from german "zins", meaning "interest")
        v::Float64 # extra value available to lend
        d::Float64 # reduction amount
        κ::Float64 # temperature
        k::Int64 # degree of network

        # constructor
        function LatticeGame(r::Float64, z::Float64)
            return new(r, z, 0.0, 0.1, 0.0, 4)
        end
        function LatticeGame(r::Float64, z::Float64, v::Float64)
            return new(r, z, v, 0.1, 0.0, 4)
        end
        function LatticeGame(r::Float64, z::Float64, v::Float64, d::Float64)
            return new(r, z, v, d, 0.0, 4)
        end
        function LatticeGame(r::Float64, z::Float64, v::Float64, d::Float64, κ::Float64, k::Int64)
            return new(r, z, v, d, κ, k)
        end

    end

    mutable struct LatticePopulation
        # new type for lattices
        # this should make it easier to pass around neighbors, etc.

        game::LatticeGame # game parameters
        N::Int64 # size of lattice
        lattice::Array{Int64, 2} # array of individuals' strategy strings
        fitnesses::Array{Float64, 2} # array of individuals' fitnesses
        reputations::Array{Int64, 2} # array of reputations
        # the above three are the only mutable attributes
        neighbors::Dict{CartesianIndex,Array{CartesianIndex,1}} # list of neighbors for each individual
        neighbor_pairs::Array{Array{CartesianIndex{2},1},1} # list of all pairs of neighbors
        generation::Int64
        reputation_payback::Bool
        reputation_coop::Bool
        verbose::Bool # turn this on for error tracking

        # constructor if game is already specified
        function LatticePopulation(
            game::LatticeGame,
            N::Int64,
            reputation_payback::Bool=false,
            reputation_coop::Bool=false,
            verbose::Bool=false
            )
            #lattice = zeros(Int64, N,N)
            # see above for explanation of next line
            #rand_init = CartesianIndices(lattice)[randperm(N^2)[1:floor(Int64,N^2*init_freq)]]
            #[lattice[x] = true for x in rand_init]
            lattice = zeros(Int64, N, N)
            reputations = ones(Int64, N, N)
            neighbors = get_neighbors(lattice)
            neighbor_pairs = get_neighbor_pairs(lattice, neighbors)
            #fitnesses = zeros(Float64, N, N)
            fitnesses = get_all_payoffs(lattice, reputations, neighbor_pairs, game, verbose)
            return new(game, N, lattice, fitnesses, reputations, neighbors, neighbor_pairs, 1, reputation_payback, reputation_coop, verbose)
        end
    end

    function randomize_lattice!(
        pop::LatticePopulation,
        strategies::Array{Int64, 1}=collect(0:7)
        )
        pop.lattice = rand(strategies, size(pop.lattice))
        pop.fitnesses = get_all_payoffs(pop.lattice, pop.reputations, pop.neighbor_pairs, pop.game)
    end

    function decompose_strategy(
        strat::Int64
        )
        # returns an array of true and false for individual's strategy components
        # largest bit: defect/cooperate
        # middle bit: hoard/lend
        # smallest bit: shirk/payback
        return [Bool(strat ÷ 2), Bool(strat%2)]
    end

    function get_all_payoffs(
        lattice::Array{Int64, 2},
        reputation::Array{Int64, 2},
        neighbor_pairs::Array{Array{CartesianIndex{2},1},1},
        game::LatticeGame,
        verbose::Bool=false
        )
        fitnesses = zeros(Float64, size(lattice))
        for (pi, pair) in enumerate(neighbor_pairs)
            strategy = lattice[pair[1]]
            neighbor_strategy = lattice[pair[2]]
            rep = reputation[pair[1]]
            neighbor_rep = reputation[pair[2]]
            payoffs = pair_PGG_payoffs(strategy, neighbor_strategy, rep, neighbor_rep, game, verbose)
            fitnesses[pair[1]] += payoffs[1]
            fitnesses[pair[2]] += payoffs[2]
        end
        return fitnesses
    end

    function update_all_payoffs!(
        pop::LatticePopulation
        )
        pop.fitnesses = get_all_payoffs(pop.lattice, pop.reputations, pop.neighbor_pairs, pop.game, pop.verbose)
    end

    # function update_all_payoffs(
    #     lattice::Array{Int64, 2},
    #     neighbors::Dict{CartesianIndex,Array{CartesianIndex,1}},
    #     game::LatticeGame
    #     )
    #     fitnesses = zeros(Float64, size(lattice))
    #     for (i, indv) in enumerate(CartesianIndices(lattice))
    #         neighbor_strategies = Int64[lattice[neighb] for neighb in neighbors[indv]]
    #         fitnesses[indv] = get_payoffs(lattice[indv], neighbor_strategies, game)
    #     end
    #     return fitnesses
    # end
    #
    # function update_indv_payoffs!(
    #     pop::LatticePopulation,
    #     indv::CartesianIndex
    #     )
    #     pop.fitnesses[indv] = update_payoffs(pop, indv)
    #     for (ni, neighb) in enumerate(pop.neighbors[indv])
    #         pop.fitnesses[neighb] = update_payoffs(pop, neighb)
    #     end
    # end
    #
    # function update_payoffs(
    #     pop::LatticePopulation,
    #     indv::CartesianIndex
    #     )
    #     strategy = pop.lattice[indv]
    #     neighbor_strategies = Int64[pop.lattice[neighb] for neighb in pop.neighbors[indv]]
    #     return get_payoffs(strategy, neighbor_strategies, pop.game)
    # end

    function pair_PGG_payoffs(
        strategy_1::Int64,
        strategy_2::Int64,
        rep_1::Int64,
        rep_2::Int64,
        gp::LatticeGame,
        verbose::Bool=false
        )

        r, z, k, d, v = gp.r, gp.z, gp.k, gp.d, gp.v # get the game parameters

        if verbose
            println("r = $r, z = $z, k = $k, v = $v, d = $d")
            println("reputations are $rep_1 and $rep_2")
        end

        loan_1 = v*(rep_1==1 ? 1.0 : 1.0-d)
        loan_2 = v*(rep_2==1 ? 1.0 : 1.0-d)

        money_1 = 0.0 + loan_1
        money_2 = 0.0 + loan_2

        if verbose
            println("indv 1 has a " * (rep_1 == 1 ? "good" : "bad") * " reputation and is given $money_1")
            println("indv 2 has a " * (rep_2 == 1 ? "good" : "bad") * " reputation and is given $money_2")
        end

        strat_1 = decompose_strategy(strategy_1) # get individual strategy components
        strat_2 = decompose_strategy(strategy_2)

        if verbose
            println("strategies are $strat_1 and $strat_2")
        end

        coop_1, payback_1 = strat_1[1], strat_1[2]
        coop_2, payback_2 = strat_2[1], strat_2[2]

        if verbose
            println("indv 1: $coop_1, $payback_1")
            println("indv 2: $coop_2, $payback_2")
        end

        PGG_pot = r*(coop_1*money_1 + coop_2*money_2)

        if verbose
            println("pot size is $PGG_pot")
        end

        payoff_1 = PGG_pot/2 + (1.0-coop_1)*money_1
        payoff_2 = PGG_pot/2 + (1.0-coop_2)*money_2

        if verbose
            if payback_1
                println("focal indv pays back $(loan_1*z)")
            end
            if payback_2
                println("neighbor pays back $(loan_2*z)")
            end
        end


        payoff_1 -= loan_1*z*payback_1
        payoff_2 -= loan_2*z*payback_2

        if verbose
            println("payoffs are $payoff_1 and $payoff_2")
        end

        return [payoff_1, payoff_2]
    end

    function get_neighbors(
        lattice::Array{Int64, 2},
        wraparound::Bool=true
        )
        # obtain the neighbors for every individual in the lattice
        # this only needs to be run once, at the beginning
        # returns a Dict with a CartesianIndex() key corresponding to the individual
        # value is a list of CartesianIndices() corresponding to neighbors

        # unnecessary but allows for non-square lattices
        (width, height) = size(lattice)

        neighbors = Dict{CartesianIndex,Array{CartesianIndex,1}}()
        for (i, indv) in enumerate(CartesianIndices(lattice))
            # this loop just controls for the possibility that an individual
            # is along an edge
            # if it is, and wraparound is set, it handles that accordingly
            # otherwise, it grants fewer neighbors
            tmp_neighbors = CartesianIndex[]
            (x, y) = (indv[1], indv[2])
            if x != 1
                push!(tmp_neighbors, CartesianIndex(x-1, y))
            end
            if y != 1
                push!(tmp_neighbors, CartesianIndex(x, y-1))
            end
            if x != width
                push!(tmp_neighbors, CartesianIndex(x+1, y))
            end
            if y != height
                push!(tmp_neighbors, CartesianIndex(x, y+1))
            end

            # i am absolutely certain there is a less ugly way to do this
            if width > 2
                if x == 1 && wraparound
                    push!(tmp_neighbors, CartesianIndex(width, y))
                end
                if x == width && wraparound
                    push!(tmp_neighbors, CartesianIndex(1, y))
                end
            end
            if height > 2
                if y == 1 && wraparound
                    push!(tmp_neighbors, CartesianIndex(x, height))
                end
                if y == height && wraparound
                    push!(tmp_neighbors, CartesianIndex(x, 1))
                end
            end

            neighbors[indv] = tmp_neighbors
        end
        return neighbors
    end

    function get_neighbor_pairs(
        lattice::Array{Int64, 2},
        neighbors::Dict{CartesianIndex,Array{CartesianIndex,1}},
        remove_duplicates::Bool=true
        )
        doublets = Array{CartesianIndex{2}, 1}[]
        for (i, indv) in enumerate(CartesianIndices(lattice))
            for (ni, neighbor) in enumerate(neighbors[indv])
                doublet = [indv, neighbor]
                if remove_duplicates
                    sort!(doublet)
                end
                # println("$doublet, $(typeof(doublet)), $(typeof(doublets))")
                if doublet ∉ doublets
                    push!(doublets, doublet)
                end
            end
        end
        return doublets
    end

    function energy_function(
        pop::LatticePopulation,
        indv::CartesianIndex,
        neighbor::CartesianIndex
        )
        # return the "energy" difference between individuals
        # this sets the probability that one of them flips to the other's strategy
        payoff_diff = pop.fitnesses[indv] - pop.fitnesses[neighbor]
        energy = 1.0/(1.0+exp(payoff_diff/pop.game.κ))
        return energy
    end

    function evolve_sync!(
        pop::LatticePopulation
        )
        # evolves the lattice one generation
        # determines an ordering, then determines payoffs
        # (including transaction costs)
        # then updates the lattice as individuals switch strategies

        new_lattice = zeros(size(pop.lattice))
        update = rand(pop.N, pop.N)

        if pop.verbose
            println("lattice looks like:")
            display(pop.lattice)
            println("update array looks like:")
            display(update)
            println("payoff array looks like:")
            display(pop.fitnesses)
        end
        update_reputations!(pop)
        for (i, indv) in enumerate(CartesianIndices(pop.lattice))
            # pick a random neighbor for each individual and compare energies
            # the energy sets the probability of a state flip
            neighbor = rand(pop.neighbors[indv])
            energy = energy_function(pop, indv, neighbor)
            if pop.verbose
                println("individual $(indv.I) chose random neighbor $(neighbor.I)")
                println("individual $(indv.I) had payoff $(pop.fitnesses[indv])")
                println("neighbor $(neighbor.I) had payoff $(pop.fitnesses[neighbor])")
                println("the energy function is $energy")
                println("the update function is $(update[indv])")
            end
            if update[indv] < energy
                new_lattice[indv] = pop.lattice[neighbor]
                if pop.verbose println("so individual $(indv.I) updates their strategy") end
            else
                new_lattice[indv] = pop.lattice[indv]
                if pop.verbose println("so individual $(indv.I) keeps their strategy") end
            end
        end
        pop.generation += 1
        pop.lattice = new_lattice
        if pop.verbose
            println("lattice now looks like:")
            display(new_lattice)
        end
        update_all_payoffs!(pop)
        return pop
    end

    function strategy_string(
        strategy::Int64
        )
        decomposed_strat = decompose_strategy(strategy)
        strat = String[]
        push!(strat, decomposed_strat[1] ? "cooperate" : "defect")
        push!(strat, decomposed_strat[2] ? "lend" : "hoard")
        push!(strat, decomposed_strat[3] ? "payback" : "shirk")
        strat_string = "$(strat[1]) $(strat[2]) $(strat[3])"
        #println(strat_string)
        return strat_string
    end

    function update_reputations!(
        pop::LatticePopulation
        )
		new_reputations = zeros(size(pop.reputations))
        for (i, indv) in enumerate(CartesianIndices(pop.lattice))
            if pop.reputation_coop && pop.reputation_payback
                new_reputations[indv] = pop.lattice[indv] ÷ 3
            elseif pop.reputation_coop
                new_reputations[indv] = pop.lattice[indv] ÷ 2
            elseif pop.reputation_payback
                new_reputations[indv] = pop.lattice[indv] % 2
				#println("updated reputations")
            end
        end
		pop.reputations = new_reputations
    end


    # function evolve_async!(
    #     pop::LatticePopulation
    #     )
    #     # evolves the lattice one generation
    #     # determines an ordering, then determines payoffs
    #     # (including transaction costs)
    #     # then updates the lattice as individuals switch strategies
    #
    #     indv = rand(collect(CartesianIndices(pop.lattice)))
    #     update = rand()
    #
    #     if pop.verbose
    #         println("lattice looks like:")
    #         display(pop.lattice)
    #     end
    #
    #     neighbor = rand(pop.neighbors[indv])
    #     energy = energy_function(pop, indv, neighbor)
    #     if update < energy
    #         pop.lattice[indv] = pop.lattice[neighbor]
    #         update_indv_payoffs!(pop, indv)
    #         #if pop.verbose println("so individual $(indv.I) updates their strategy") end
    #     end
    #     pop.generation += 1
    #     #pop.lattice = new_lattice
    #     if pop.verbose
    #         println("lattice now looks like:")
    #         #display(new_lattice)
    #     end
    #     return pop
    # end

# final end statement to close the module
end
