#!/usr/bin/env julia

## compute_payoffs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Compute payoffs for a simple social credit example.

struct GameParams
    # a dummy type for easily storing/passing game variables
    r::Float64 # synergy factor for public goods game
    z::Float64 # interest on loans (from german "zins", meaning "interest")
    k::Int64 # number of neighbors per individual

    function GameParams(r::Float64, z::Float64, k::Int64)
        return new(r, z, k)
    end

    function GameParams(r::Float64, z::Float64)
        # constructor that assigns default value of k = 4
        return new(r, z, 4)
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

function get_payoffs(
    strategy::Int64,
    neighbor_strategies::Array{Int64, 1},
    gp::GameParams
    )
    # returns the payoffs for an individual after each step of lending and PGG.
    # strategy: individual strategy, ranges from 0 to 7
    # largest bit is defect/coop, then hoard/lend, then shirk/payback
    # for example, 4 would be cooperator who hoards and shirks
    # 3 would be a defector who lends and pays back
    # neighbor_strategies is the same

    r, z, k = gp.r, gp.z, gp.k # get the game parameters

    loan = 1.0 # the amount that is loaned

    initial_money = 1.0 # the amount of money everyone starts with
    # regrettably, this turns out to matter

    payoff = 0.0
    money = initial_money
    neighbor_money = initial_money*ones(k)

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

    println("money after loan is $money")
    println("neighbor money is $neighbor_money")

    money = min(0.0, money)
    neighbor_money = min.(0.0, neighbor_money)
    # this step has to be done because people cannot contribute "negative" money to the PGG

    # play the PGG
    # cooperators invest all their money in the game
    # defectors invest none
    # we of course have to remember to subtract the focal individual's contribution
    # if they are a cooperator
    payoff += (r*(coop*money + sum(neighbor_money' * neighbor_strats[:,1]))/k - coop)
    println("payoff is $payoff")

    # do the payback step
    for i in 1:k
        neighbor_lend = neighbor_strats[i,2]
        neighbor_payback = neighbor_strats[i,3]
        payoff -= neighbor_lend*payback*z # if the neighbor lent us money and we paid them back
        payoff += lend*neighbor_payback*z # if we lent the neighbor money and they pay us back
    end
    return payoff
end

# get the payoffs for a simple example: one neighbor.
# first argument: strategy
# second argument: array of neighbor strategies
# third: r, z, k (synergy factor, interest, number of neighbors)
get_payoffs(4, [4], GameParams(1.1, 1.2, 1))
