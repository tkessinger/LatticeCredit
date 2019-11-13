#!/usr/bin/env julia

## cost.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test LatticeLending, a central-bank version of LatticeCredit
## with a pairwise PGG.
## Obtain frequency trajectories and plot them.

using Random
using Revise
using LatticeLending

using PyPlot

N = 40 # size of lattice

#r_vals = collect(1:0.3:4) # values of synergy parameter r
#v_vals = [1.0,2.0,4.0] # loanable amount v
#v_vals = [2.0]
# v_vals = [0] # uncomment this to remove lending
num_gens = 200

v = 2.0
z_vals = [1.2, 1.5, 1.8]#, 1.6, 2.0] # interest on loans
r_vals = [3.5, 3.8, 4.0]
#z = 1.5
#r = 1.5
κ = 0.1 # temperature parameter
d = 0.8

function tracking_label(reputations)
    if reputations[1] && reputations[2]
        return "both"
    elseif reputations[1]
        return "payback"
    elseif reputations[2]
        return "coop"
    else
        return "neither"
    end
end
fig, axs = plt.subplots(2*length(r_vals), 2*length(z_vals),
    figsize=(12, 12), sharey="row", sharex="col")
for reputation_payback in [false, true]
    for reputation_coop in [false, true]
        for (zi, z) in enumerate(z_vals)
            for (ri, r) in enumerate(r_vals)
                frequencies = zeros(Float64, num_gens, 4)

                game = LatticeGame(r, z, v, d, κ, 4)
                pop = LatticePopulation(game, N, reputation_payback, reputation_coop)
                #println(pop)
                randomize_lattice!(pop, collect(0:3))
                for i in 1:num_gens
                    evolve_sync!(pop)
                    [frequencies[i, j+1] += sum(x==j for x in pop.lattice)*1.0/N^2 for j in 0:3]
                end

                ax = axs[ri+length(r_vals)*reputation_coop, zi+length(z_vals)*reputation_payback]

                [ax.plot(collect(1:num_gens), frequencies[:,j+1], label="type $j") for j in 0:3]

                if ri == length(r_vals) && reputation_coop == true
                    ax.set_xlabel("time")
                end
                if zi == 1 && reputation_payback == false
                    ax.set_ylabel(L"\rho")
                end
                ax.set_ylim([0,1])
                plt.legend(loc=2)
                ax.set_title("$(tracking_label([reputation_payback, reputation_coop])), r = $r, z = $z")
            end
        end
    end
end
fig.tight_layout()
fig.suptitle("d = $d, N = $N")
plt.subplots_adjust(top=0.94)
display(fig)
