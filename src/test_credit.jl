#!/usr/bin/env julia

## cost.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test LatticeCredit on a small lattice.
## Still to-do: implement more detailed bug checking.

using Random
using Revise
using LatticeCredit

using PyPlot

z = 1.3 # interest on loans
r = 1.5 # synergy parameter
l = 1.0 # loan amount
m = 4.0 # initial bankroll per individual
κ = 0.1 # temperature parameter

histories = []
overall_freqs = []

numgens = 1000
N = 5 # size of lattice

game = LatticeGame(r, z, 1.0, 4.0, κ, 4)
pop = LatticePopulation(game, N)

println("$(pop.lattice)")
#for i in 1:numgens
# while the lattice consists of a mix of types
while !(all(y->y==pop.lattice[1], pop.lattice))
    evolve_async!(pop)
    println("$(pop.lattice)")
end
println("$(pop.lattice)")

#clr = ["r", "b", "y", "g", "m"]

#c_list = [0.05, 0.25, 0.45, 0.65, 0.85]
#fig = plt.figure()
# for (ci, c) in enumerate(c_list)
#
#     println("evolving with transaction cost $c")
#
#     pop = LatticePopulation(b, c, κ, N)
#
#     history = []
#
#     push!(history, pop.lattice)
#     for i in 1:numgens
#         evolve!(pop)
#         push!(history, pop.lattice)
#     end
#
#     push!(histories, history)
#
#     coop_freqs = 1.0/N^2*[sum(x) for x in history]
#     push!(overall_freqs, coop_freqs)
#
#     plt.plot(collect(1:numgens+1), coop_freqs, label="c = $c",
#         c = clr[ci])
#     display(fig)
# end
#
# ax = fig.gca()
# ax.legend(loc=2)
# ax.set_xscale("log")
# ax.set_xlabel("time")
# ax.set_ylabel(L"\rho_c")
# ax.set_xlim([1,numgens])
# fig.tight_layout()
# display(fig)

# anim = @animate for i=1:length(history)
#     heatmap(history[i])
# end

#
# macro_lattices = []
# coop_freqs = []
#
# fig = figure()
#
# c_vals = [0.05, 0.25, 0.45, 0.65, 0.85]
# for c in c_vals
#     lattices = []
#     push!(lattices, lattice)
#     for i in 1:numgens-1
#         #println("$lattice")
#         println("$i")
#         global lattice = evolve(lattice, payoff_matrix, c, κ)
#         #imshow(lattice)
#         gcf()
#         push!(lattices, lattice)
#     end
#
#     coop_freq = [0.5]
#     for i in 1:numgens
#         push!(coop_freq, sum(lattices[i])/N^2)
#     end
#     plot(coop_freq, label="c = $c")
#     push!(macro_lattices, lattices)
#     push!(coop_freqs, coop_freq)
# end
# xscale("log")
# xlim([1,numgens])
# legend(loc=2)
# xlabel("time")
# ylabel("cooperator frequency")
# fig.savefig("test_coop_freqs.pdf")
#
# for (ci, c) in enumerate(c_vals)
#     figure()
#     anim = @animate for i=1:length(macro_lattices[ci])
#         heatmap(macro_lattices[ci][i])
#     end
#     gif(anim, "lattice_c_" * "$c" * "_fps15_take2.gif", fps = 15)
# end
