#!/usr/bin/env julia

## LatticeCredit.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating loan-PGG-payback game on a lattice.

module LatticeCredit

    using Random

    export LatticeGame, LatticePopulation
    export evolve_sync!, evolve_async!
    export decompose_strategy, strategy_string, get_payoffs, get_payoffs_fast

    struct LatticeGame
        # new type for "game" parameters
        # this should minimize the amount of passing around that needs to be done

        r::Float64 # benefit to defecting
        z::Float64 # interest on loans (from german "zins", meaning "interest")
        l::Float64 # loan amount
        m::Float64 # initial money
        κ::Float64 # temperature
        k::Int64 # degree of network

        # constructor
        function LatticeGame(r::Float64, z::Float64)
            return new(r, z, 1.0, 4.0, 0.1, 4)
        end
        function LatticeGame(r::Float64, z::Float64, l::Float64, m::Float64, κ::Float64, k::Int64)
            return new(r, z, l, m, κ, k)
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

        # constructor
        # function LatticePopulation(
        #     r::Float64 # benefit to defecting
        #     z::Float64 # interest on loans (from german "zins", meaning "interest")
        #     l::Float64 # loan amount
        #     m::Float64 # initial money
        #     N::Int64,
        #     init_freq::Float64=0.5,
        #     verbose::Bool=false
        #     )
        #     lattice = zeros(Int64, N, N)
        #     # next line makes a list of the lattice's CartesianIndices,
        #     # randomly orders them,
        #     # takes the first init_freq (fraction) thereof,
        #     # then changes the corresponding lattice points to true
        #     rand_init = CartesianIndices(lattice)[randperm(N^2)[1:floor(Int64,N^2*init_freq)]]
        #     [lattice[x] = true for x in rand_init]
        #     neighbors = get_neighbors(lattice)
        #     neighbor_pairs = get_neighbor_pairs(lattice, neighbors)
        #     return new(LatticeGame(b, c, κ), N, lattice, neighbors, neighbor_pairs, 1, verbose)
        # end

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
            lattice = rand(0:7, N, N)
            neighbors = get_neighbors(lattice)
            neighbor_pairs = get_neighbor_pairs(lattice, neighbors)
            #fitnesses = zeros(Float64, N, N)
            fitnesses = update_all_payoffs(lattice, neighbors, game)
            return new(game, N, lattice, fitnesses, neighbors, neighbor_pairs, 1, verbose)
        end
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

    function who_pays(
        pop::LatticePopulation
        )
        # borrowed from LatticeCoop: function for determining who pays transaction costs
        # could easily be adapted if we decide to make lending binary
        # i.e., allow each individual to issue only one loan
        initiator = Dict{Array{CartesianIndex{2}, 1}, Bool}()
        for (npi, neighbor_pair) in enumerate(pop.neighbor_pairs)
            initiator[neighbor_pair] = rand(Bool)
        end
        return initiator
    end

    # function get_payoffs(
    #     pop::LatticePopulation,
    #     initiator::Dict{Array{CartesianIndex{2},1}, Bool}
    #     )
    #     # determine the payoff for every individual in the lattice
    #     payoffs = zeros(size(pop.lattice))
    #     for (i, indv) in enumerate(CartesianIndices(pop.lattice))
    #         # iterate over each individual
    #         tmp_payoffs = 0
    #         lattice_val = pop.lattice[indv]
    #         for (ni, neighbor) in enumerate(pop.neighbors[indv])
    #             # iterate over all their neighbors
    #             # get everyone's strategy
    #             neighbor_val = pop.lattice[neighbor]
    #             indv_strat = get_strategy(lattice_val)
    #             neighbor_strat = get_strategy(neighbor_val)
    #
    #             # decide who has to pay the transaction cost
    #             neighbor_pair = sort([indv, neighbor])
    #             pay_cost = xor(initiator[neighbor_pair], issorted([indv, neighbor]))
    #             tmp_payoff = transpose(indv_strat)*pop.game.A*neighbor_strat - pop.game.c*(pay_cost)
    #
    #             if pop.verbose
    #                 #println("individual $(indv.I) has strategy $(strategy_string(indv_strat)) and ordering $(ordering[indv])")
    #                 #println("neighbor $(neighbor.I) has strategy $(strategy_string(neighbor_strat)) and ordering $(ordering[neighbor])")
    #                 println("individual should pay transaction cost: $pay_cost")
    #                 println("payoff is $tmp_payoff")
    #             end
    #             tmp_payoffs += tmp_payoff
    #         end
    #         payoffs[indv] = tmp_payoffs
    #         if pop.verbose println("individual $(indv.I) had net payoff $tmp_payoffs") end
    #     end
    #     return payoffs
    # end

    function update_all_payoffs!(
        pop::LatticePopulation
        )
        for (i, indv) in enumerate(CartesianIndices(pop.lattice))
            pop.fitnesses[indv] = update_payoffs(indv, pop)
        end
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

    function get_payoffs_fast(
        strategy::Int64,
        neighbor_strategies::Array{Int64, 1},
        gp::LatticeGame,
        verbose::Bool=false
        )
        # a faster implementation of get_payoffs.
        # currently this is bugged.


        r, z, k, m, l = gp.r, gp.z, gp.k, gp.m, gp.l # get the game parameters

        payoff = 0.0
        money = m
        loan = l
        neighbor_money = m*ones(Float64, k)

        strat = decompose_strategy(strategy) # get individual strategy components
        neighbor_strats = zeros(Bool, k, 3) # later we will do the same for neighbors
        for i in 1:k
            neighbor_strats[i,:] = decompose_strategy(neighbor_strategies[i])
        end

        coop, lend, payback = strat[1], strat[2], strat[3]

        payoff += r/(k+1)*(coop*(m - k*lend + sum([neighbor_strats[x,1]*(m + lend - neighbor_strats[x,2]) for x in 1:k])))
        payoff -= coop*(m - k*lend + sum([neighbor_strats[x,1] * (m + lend - neighbor_strats[x,2]) for x in 1:k]))
        payoff += z*sum(lend*neighbor_strats[:,3] - payback*neighbor_strats[:,2])

        return payoff
    end

    function get_payoffs(
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
        r, z, k, m, l = gp.r, gp.z, gp.k, gp.m, gp.l # get the game parameters

        payoff = 0.0
        money = m
        loan = l
        neighbor_money = m*ones(Float64, k)

        strat = decompose_strategy(strategy) # get individual strategy components
        neighbor_strats = zeros(Bool, k, 3) # later we will do the same for neighbors
        for i in 1:k
            neighbor_strats[i,:] = decompose_strategy(neighbor_strategies[i])
        end

        coop, lend, payback = strat[1], strat[2], strat[3]

        # issue loans and decompose neighbor strategies
        for i in 1:k
            neighbor_lend = neighbor_strats[i,2] # second bit is lend/hoard
            money += loan*neighbor_lend # if neighbor is a lender, get some money
            neighbor_money[i] -= loan*neighbor_lend # they lose the same amount

            money -= loan*lend # if focal individual is a lender, lose some money
            neighbor_money[i] += loan*lend # neighbor gets the same amount
        end

        money = min(0.0, money)
        neighbor_money = min.(0.0, neighbor_money)
        # this step has to be done because people cannot contribute "negative" money to the PGG

        # play the PGG
        # cooperators invest all their money in the game
        # defectors invest none
        # we of course have to remember to subtract the focal individual's contribution
        # if they are a cooperator
        payoff += (r*(coop*money + sum(neighbor_money' * neighbor_strats[:,1]))/(k+1) - coop)
        if verbose
            println("payoff is $payoff")
        end

        # do the payback step
        for i in 1:k
            neighbor_lend = neighbor_strats[i,2]
            neighbor_payback = neighbor_strats[i,3]
            payoff -= neighbor_lend*payback*z # if the neighbor lent us money and we paid them back
            payoff += lend*neighbor_payback*z # if we lent the neighbor money and they pay us back
        end
        return payoff
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

        #ordering = assign_transactions(pop)
        #ordering = rand(pop.N, pop.N)
        #initiator = who_pays(pop)
        #payoffs = get_payoffs(pop, initiator)
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
        # strat[1] = (decomposed_strat[1] ? "defect" : "cooperate")
        # strat[2] = (decomposed_strat[2] ? "hoard" : "lend")
        # strat[3] = (decomposed_strat[3] ? "shirk" : "payback")
        println("$(strat[1]) $(strat[2]) $(strat[3])")
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
