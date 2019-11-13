#!/usr/bin/env julia

## test_bank_compare_tracking.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test LatticeLending, a central-bank version of LatticeCredit with a pairwise PGG.
## Compare different reputation tracking methods across different values of r and z.

using Random
using Revise
using LatticeLending

using PyPlot

N = 20 # size of lattice

r_vals = collect(1:0.2:5) # values of synergy parameter r
#v_vals = [1.0,2.0,4.0] # loanable amount v
v_vals = [2.0]
# v_vals = [0] # uncomment this to remove lending
num_trials = 25
num_gens = 200

v = 2.0
z_vals = [1.0, 1.2, 1.5, 2.0, 3.0] # interest on loans
κ = 0.1 # temperature parameter
d = 0.9

function sim_and_get_frequencies(
    r_vals::Array{Float64, 1},
    z_vals::Array{Float64, 1},
    N::Int64,
    v::Float64,
    d::Float64,
    κ::Float64,
    num_trials::Int64,
    num_gens::Int64,
    reputation_tracking::Array{Bool, 1}
    )

    c_frac = zeros(Float64, length(r_vals), length(z_vals))
    payback_frac = zeros(Float64, length(r_vals), length(z_vals))
    frequencies = zeros(Float64, length(r_vals), length(z_vals), 4)

    reputation_payback, reputation_coop = reputation_tracking[1], reputation_tracking[2]

    for (ri, r) in enumerate(r_vals)
        for (zi, z) in enumerate(z_vals)
            println("simulating PGG with r = $r and z = $z")
            println("reputations: $reputation_tracking")

            for i in 1:num_trials
                # initialize game and population
                game = LatticeGame(r, z, v, d, κ, 4)
                pop = LatticePopulation(game, N, reputation_payback, reputation_coop)
                #println(pop)
                randomize_lattice!(pop, collect(0:3))

                for j in 1:num_gens
                    evolve_sync!(pop)
                end
                # n.b.: this is integer division
                c_frac[ri,zi] += sum(x÷2 == 1 for x in pop.lattice)*1.0/num_trials/N^2
                payback_frac[ri, zi] += sum(x%2 == 1 for x in pop.lattice)*1.0/num_trials/N^2
                [frequencies[ri,zi,j+1] += sum(x==j for x in pop.lattice)*1.0/num_trials/N^2 for j in 0:3]
            end
        println("ρ_c at r = $r and z = $z is $(c_frac[ri,zi])")
        end
    end
    return (c_frac, payback_frac, frequencies)
end

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

function plot_frequencies(r_vals, z_vals, N, v, d, κ, num_trials, num_gens)

    for (zi, z) in enumerate(z_vals)
        fig = plt.figure()
        for tracking_params in collect(Base.product([false, true], [false, true]))
            tracking_params = collect(tracking_params)
            c_frac, payback_frac, frequencies = sim_and_get_frequencies(r_vals, [z], N, v, d, κ,
                num_trials, num_gens, tracking_params)
            plt.plot(r_vals, c_frac[:,1], label=tracking_label(tracking_params))
        end
        plt.xlabel(L"r")
        plt.ylabel(L"\rho_{\mathrm{cooperate}}")
        plt.ylim([0,1])
        plt.legend(loc=2)
        plt.title("z = $z")
        fig.tight_layout()
        display(fig)
    end



    fig = plt.figure()
    for (zi, z) in enumerate(z_vals)
        fig = plt.figure()

        for tracking_params in collect(Base.product([false, true], [false, true]))
            tracking_params = collect(tracking_params)
            c_frac, payback_frac, frequencies = sim_and_get_frequencies(r_vals, [z], N, v, d, κ,
                num_trials, num_gens, tracking_params)
            plt.plot(r_vals, payback_frac[:,1], label=tracking_label(tracking_params))
        end
        plt.xlabel(L"r")
        plt.ylabel(L"\rho_{\mathrm{payback}}")
        plt.ylim([0,1])
        plt.legend(loc=2)
        plt.title("z = $z")
        fig.tight_layout()
        display(fig)
    end
end

function subplot_frequencies(r_vals, z_vals, N, v, d, κ, num_trials, num_gens)

    fig, axs = plt.subplots(2, length(z_vals),
        figsize=(14, 7), sharey="row", sharex="col")

    for (zi, z) in enumerate(z_vals)
        for tracking_params in collect(Base.product([false, true], [false, true]))
            tracking_params = collect(tracking_params)
            c_frac, payback_frac, frequencies = sim_and_get_frequencies(r_vals, [z], N, v, d, κ,
                num_trials, num_gens, tracking_params)
            ax = axs[1, zi]
            if zi == 1
                ax.set_ylabel(L"\rho_{\mathrm{cooperate}}")
            end
            ax.plot(r_vals, c_frac[:,1], label=tracking_label(tracking_params))
            ax.set_ylim([0,1])
            ax.set_title("z = $z")


            ax = axs[2, zi]
            if zi == 1
                ax.set_ylabel(L"\rho_{\mathrm{payback}}")
            end
            ax.plot(r_vals, payback_frac[:, 1], label=tracking_label(tracking_params))
            ax.set_xlabel(L"r")
            ax.set_ylim([0,1])
        end
    end
    plt.legend(loc=2)
    fig.tight_layout()
    plt.suptitle("d = $d, N = $N")
    plt.subplots_adjust(top=0.92)
    display(fig)
    return fig
end

for d in [0.9, 0.5, 0.1, 0.0]
    fig = subplot_frequencies(r_vals, z_vals, N, v, d, κ, num_trials, num_gens)
    plt.savefig("figures/final_freqs_N" * "$N" * "_d0" * "$(floor(Int64, d*10))" * ".pdf")
end
