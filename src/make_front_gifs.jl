#!/usr/bin/env julia

## make_front_gifs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Make some GIFs of payback/shirk fronts in the high-cooperation regime.

using Random
using Revise
using LatticeLending

using Plots
pyplot()

N = 50 # size of lattice
num_gens = 400

v = 1.0
r = 3.5
z = 1.25
κ = 0.1 # temperature parameter
d = 0.8

reputation_payback = true
reputation_coop = false

histories = []
reputations = []

game = LatticeGame(r, z, v, d, κ, 4)
pop = LatticePopulation(game, N, reputation_payback, reputation_coop)
randomize_lattice!(pop, collect(0:3))

for i in 1:num_gens
    evolve_sync!(pop)
    if i > 200
        push!(histories, pop.lattice .- 2 + pop.reputations*2)
    end
end


anim = @animate for i=1:length(histories)
    heatmap(histories[i], aspect_ratio=1, clim=(0,3))
end

gif(anim, "figures/rep_fronts_z_$(z)_r_$(r)_d_$(d).gif", fps = 10)

# key:
# black (0): shirk, bad reputation
# purple (1): payback, bad reputation
# orange (2): shirk, good reputation
# cream (3): payback, good reputation
