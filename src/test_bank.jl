#!/usr/bin/env julia

## cost.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test LatticeLending, a central-bank version of LatticeCredit with a pairwise PGG.

using Random
using Revise
using LatticeLending

using PyPlot

num_gens = 2 # number of generations to simulate
N = 3 # size of lattice

r, z, v, d, κ = 5.5, 1.2, 1.0, 1.0, 0.1

reputation_payback = true
reputation_coop = true

# initialize game and population
game = LatticeGame(r, z, v, d, κ, 4)
pop = LatticePopulation(game, N, reputation_payback, reputation_coop)
randomize_lattice!(pop, collect(0:3))

println(pop.lattice)
println(pop.reputations)
println(pop.fitnesses)



for j in 1:num_gens
    evolve_sync!(pop)
    println(pop.lattice)
    println(pop.reputations)
    println(pop.fitnesses)
    println(pop.generation)
end
