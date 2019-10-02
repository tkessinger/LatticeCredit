#!/usr/bin/env julia

## cost.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Iterate LatticeCredit over various values of r with no lending.

using Random
using Revise
using LatticeCredit

using PyPlot

numgens = 200 # number of generations to simulate
N = 20 # size of lattice

r_vals = collect(3:0.6:6) # values of synergy parameter r
v_vals = [0,1.0,4.0] # loanable amount v
# v_vals = [0] # uncomment this to remove lending
num_trials = 10
c_frac = zeros(Float64, length(r_vals), length(v_vals))
frequencies = zeros(Float64, length(r_vals), length(v_vals), 8)

z = 1.3 # interest on loans
κ = 0.1 # temperature parameter

for (ri, r) in enumerate(r_vals)
    for (vi, v) in enumerate(v_vals)
        println("simulating PGG with r = $r and v = $v")
        for i in 1:num_trials
            # initialize game and population
            game = LatticeGame(r, z, v, κ, 4)
            pop = LatticePopulation(game, N)
            randomize_lattice!(pop, [0,7])

            for j in 1:numgens
                evolve_sync!(pop)
            end
            # n.b.: this is integer division
            c_frac[ri,vi] += sum(x÷4 == 1 for x in pop.lattice)*1.0/num_trials/N^2
            [frequencies[ri,vi,j+1] += sum(x==j for x in pop.lattice)*1.0/num_trials/N^2 for j in 0:7]
        end
        println("ρ_c at r = $r and v = $v is $(c_frac[ri,vi])")
    end
end
fig = plt.figure()
for (vi, v) in enumerate(v_vals)
    plt.plot(r_vals, c_frac[:,vi],label="v = $v")
end
plt.xlabel(L"r")
plt.ylabel(L"\rho_c")
plt.legend(loc=2)
fig.tight_layout()
display(fig)
