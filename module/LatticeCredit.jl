#!/usr/bin/env julia

## LatticeCredit.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating loan-PGG-payback game on a lattice.

module LatticeCredit

    using Random

    export LatticeGame, LatticePopulation
    export evolve_sync!, evolve_async!, randomize_lattice!
    export decompose_strategy, strategy_string
    export get_payoffs, get_payoffs_fast, get_group_payoffs

    struct LatticeGame
        # new type for "game" parameters
        # this should minimize the amount of passing around that needs to be done

        r::Float64 # synergy parameter
        z::Float64 # interest on loans (from german "zins", meaning "interest")
        v::Float64 # extra value available to lend
        κ::Float64 # temperature
        k::Int64 # degree of network

        # constructor
        function LatticeGame(r::Float64, z::Float64)
            return new(r, z, 0.0, 0.1, 4)
        end
        function LatticeGame(r::Float64, z::Float64, v::Float64)
            return new(r, z, v, 0.1, 4)
        end
        function LatticeGame(r::Float64, z::Float64, v::Float64, κ::Float64, k::Int64)
            return new(r, z, v, κ, k)
        end

    end

    mutable struct LatticePopulation
        # new type for lattices
        # this should make it easier to pass around neighbors, etc.

        game::LatticeGame # game parameters
        N::Int64 # size of lattice
        lattice::Array{Int64, 2} # array of individuals' strategy strings
        fitnesses::Array{Float64, 2} # array of individuals' fitnesses
        # the above two are the only mutable attributes
        neighbors::Dict{CartesianIndex,Array{CartesianIndex,1}} # list of neighbors for each individual
        neighbor_pairs::Array{Array{CartesianIndex,1},1} # list of all pairs of neighbors
        generation::Int64
        verbose::Bool # turn this on for error tracking

        # constructor if game is already specified
        function LatticePopulation(
            game::LatticeGame,
            N::Int64,
            verbose::Bool=false
            )
            #lattice = zeros(Int64, N,N)
            # see above for explanation of next line
            #rand_init = CartesianIndices(lattice)[randperm(N^2)[1:floor(Int64,N^2*init_freq)]]
            #[lattice[x] = true for x in rand_init]
            lattice = zeros(Int64, N, N)
            neighbors = get_neighbors(lattice)
            neighbor_pairs = get_neighbor_pairs(lattice, neighbors)
            #fitnesses = zeros(Float64, N, N)
            fitnesses = get_all_payoffs(lattice, neighbors, game)
            return new(game, N, lattice, fitnesses, neighbors, neighbor_pairs, 1, verbose)
        end
    end

    function randomize_lattice!(
        pop::LatticePopulation,
        strategies::Array{Int64, 1}=collect(0:7)
        )
        pop.lattice = rand(strategies, size(pop.lattice))
        pop.fitnesses = get_all_payoffs(pop.lattice, pop.neighbors, pop.game)
    end

    function decompose_strategy(
        strat::Int64
        )
        # returns an array of true and false for individual's strategy components
        # largest bit: defect/cooperate
        # middle bit: hoard/lend
        # smallest bit: shirk/payback
        return [Bool(strat ÷ 4), Bool(strat%4 ÷ 2), Bool(strat%4%2)]
    end

    function get_all_payoffs(
        lattice::Array{Int64, 2},
        neighbors::Dict{CartesianIndex,Array{CartesianIndex,1}},
        game::LatticeGame
        )
        fitnesses = zeros(Float64, size(lattice))
        for (i, indv) in enumerate(CartesianIndices(lattice))
            strategy = lattice[indv]
            neighbor_strategies = Int64[lattice[neighb] for neighb in neighbors[indv]]
            payoffs = get_group_payoffs(strategy, neighbor_strategies, game)
            fitnesses[indv] += payoffs[1]
            for (ni, neighb) in enumerate(neighbors[indv])
                fitnesses[neighb] += payoffs[2][ni]
            end
        end
        return fitnesses
    end

    function update_all_payoffs!(
        pop::LatticePopulation
        )
        pop.fitnesses = get_all_payoffs(pop.lattice, pop.neighbors, pop.game)
    end

    function update_all_payoffs(
        lattice::Array{Int64, 2},
        neighbors::Dict{CartesianIndex,Array{CartesianIndex,1}},
        game::LatticeGame
        )
        fitnesses = zeros(Float64, size(lattice))
        for (i, indv) in enumerate(CartesianIndices(lattice))
            neighbor_strategies = Int64[lattice[neighb] for neighb in neighbors[indv]]
            fitnesses[indv] = get_payoffs(lattice[indv], neighbor_strategies, game)
        end
        return fitnesses
    end

    function update_indv_payoffs!(
        pop::LatticePopulation,
        indv::CartesianIndex
        )
        pop.fitnesses[indv] = update_payoffs(pop, indv)
        for (ni, neighb) in enumerate(pop.neighbors[indv])
            pop.fitnesses[neighb] = update_payoffs(pop, neighb)
        end
    end

    function update_payoffs(
        pop::LatticePopulation,
        indv::CartesianIndex
        )
        strategy = pop.lattice[indv]
        neighbor_strategies = Int64[pop.lattice[neighb] for neighb in pop.neighbors[indv]]
        return get_payoffs(strategy, neighbor_strategies, pop.game)
    end

    function get_group_payoffs(
        strategy::Int64,
        neighbor_strategies::Array{Int64, 1},
        gp::LatticeGame,
        verbose::Bool=false
        )
        # returns the payoffs for an individual after each step of lending and PGG.
        # strategy: individual strategy, ranges from 0 to 7
        # largest bit is defect/coop, then hoard/lend, then shirk/payback
        # for example, 4 would be cooperator who hoards and shirks
        # 3 would be a defector who lends and pays back
        # neighbor_strategies is the same
        if verbose
            println("strategies are $strategy and $neighbor_strategies")
        end
        r, z, k, v = gp.r, gp.z, gp.k, gp.v # get the game parameters

        focal_payoff = 0.0
        neighbor_payoffs = zeros(Float64, length(neighbor_strategies))
        focal_money = 1.0 + v
        neighbor_money = ones(Float64, k)

        strat = decompose_strategy(strategy) # get individual strategy components
        neighbor_strats = zeros(Bool, k, 3) # later we will do the same for neighbors
        for i in 1:k
            neighbor_strats[i,:] = decompose_strategy(neighbor_strategies[i])
        end
        neighbor_coop = neighbor_strats[:,1]
        neighbor_lend = neighbor_strats[:,2]
        neighbor_payback = neighbor_strats[:,3]

        coop, lend, payback = strat[1], strat[2], strat[3]

        loan_amount = v/k

        # issue loans and decompose neighbor strategies
        for i in 1:k
            focal_money -= lend*loan_amount
            neighbor_money[i] += lend*loan_amount
        end

        # play the PGG
        # cooperators invest all their money in the game
        # defectors invest none
        # we of course have to remember to subtract the focal individual's contribution
        # if they are a cooperator
        PGG_pot = r*((coop)*focal_money + sum(neighbor_coop .* neighbor_money))
        focal_payoff += PGG_pot/(k+1) - coop
        [neighbor_payoffs[i] += PGG_pot/(k+1) - neighbor_coop[i] for i in 1:k]

        # do the payback step
        for i in 1:k
            # if we lent the neighbor money, and they don't shirk,
            # then they pay us back with interest
            focal_payoff += lend*neighbor_payback[i]*loan_amount*z
            neighbor_payoffs[i] -= lend*neighbor_payback[i]*loan_amount*z
        end
        return [focal_payoff, neighbor_payoffs]
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
            display(payoffs)
        end

        for (i, indv) in enumerate(CartesianIndices(pop.lattice))
            # pick a random neighbor for each individual and compare energies
            # the energy sets the probability of a state flip
            neighbor = rand(pop.neighbors[indv])
            energy = energy_function(pop, indv, neighbor)
            if pop.verbose
                println("individual $(indv.I) chose random neighbor $(neighbor.I)")
                println("individual $(indv.I) had payoff $(payoffs[indv])")
                println("neighbor $(neighbor.I) had payoff $(payoffs[neighbor])")
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

    function evolve_async!(
        pop::LatticePopulation
        )
        # evolves the lattice one generation
        # determines an ordering, then determines payoffs
        # (including transaction costs)
        # then updates the lattice as individuals switch strategies

        indv = rand(collect(CartesianIndices(pop.lattice)))
        update = rand()

        if pop.verbose
            println("lattice looks like:")
            display(pop.lattice)
        end

        neighbor = rand(pop.neighbors[indv])
        energy = energy_function(pop, indv, neighbor)
        if update < energy
            pop.lattice[indv] = pop.lattice[neighbor]
            update_indv_payoffs!(pop, indv)
            #if pop.verbose println("so individual $(indv.I) updates their strategy") end
        end
        pop.generation += 1
        #pop.lattice = new_lattice
        if pop.verbose
            println("lattice now looks like:")
            #display(new_lattice)
        end
        return pop
    end

# final end statement to close the module
end
