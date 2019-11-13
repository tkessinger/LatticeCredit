#!/usr/bin/env julia

## test_bank_variable_r.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test LatticeLending, a central-bank version of LatticeCredit with a pairwise PGG.

using Random
using Revise
using LatticeLending

using PyPlot

N = 10 # size of lattice

r_vals = collect(1:0.3:4) # values of synergy parameter r
#v_vals = [1.0,2.0,4.0] # loanable amount v
v_vals = [2.0]
# v_vals = [0] # uncomment this to remove lending
num_trials = 50
num_gens = 400
c_frac = zeros(Float64, length(r_vals), length(v_vals))
payback_frac = zeros(Float64, length(r_vals), length(v_vals))
frequencies = zeros(Float64, length(r_vals), length(v_vals), 4)

v = 2.0
z_vals = [1.2, 1.6, 2.0] # interest on loans
κ = 0.1 # temperature parameter
d = 0.5

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

function plot_frequencies(r_vals, z_vals, N, v, d, κ, num_trials, num_gens,
    reputation_payback, reputation_coop)

    c_frac, payback_frac, frequencies = sim_and_get_frequencies(r_vals, z_vals, N, v, d, κ,
        num_trials, num_gens, [reputation_payback, reputation_coop])

end

reputation_payback = true
reputation_coop = false

c_frac, payback_frac, frequencies = sim_and_get_frequencies(r_vals, z_vals, N, v, d, κ,
    num_trials, num_gens, [true, false])

fig = plt.figure()
for (zi, z) in enumerate(z_vals)
    plt.plot(r_vals, c_frac[:,zi],label="z = $z")
end
plt.xlabel(L"r")
plt.ylabel(L"\rho_{\mathrm{cooperate}}")
plt.legend(loc=2)
plt.title("reputation tracked on payback")
fig.tight_layout()
display(fig)


fig = plt.figure()
for (zi, z) in enumerate(z_vals)
    plt.plot(r_vals, payback_frac[:,zi],label="z = $z")
end
plt.xlabel(L"r")
plt.ylabel(L"\rho_{\mathrm{payback}}")
plt.legend(loc=2)
plt.title("reputation tracked on payback")
fig.tight_layout()
display(fig)

c_frac, payback_frac, frequencies = sim_and_get_frequencies(r_vals, z_vals, N, v, d, κ,
    num_trials, num_gens, [false, true])

fig = plt.figure()
for (zi, z) in enumerate(z_vals)
    plt.plot(r_vals, c_frac[:,zi],label="z = $z")
end
plt.xlabel(L"r")
plt.ylabel(L"\rho_{\mathrm{cooperate}}")
plt.legend(loc=2)
plt.title("reputation tracked on cooperation")
fig.tight_layout()
display(fig)


fig = plt.figure()
for (zi, z) in enumerate(z_vals)
    plt.plot(r_vals, payback_frac[:,zi],label="z = $z")
end
plt.xlabel(L"r")
plt.ylabel(L"\rho_{\mathrm{payback}}")
plt.legend(loc=2)
plt.title("reputation tracked on cooperation")
fig.tight_layout()
display(fig)

c_frac, payback_frac, frequencies = sim_and_get_frequencies(r_vals, z_vals, N, v, d, κ,
    num_trials, num_gens, [true, true])

fig = plt.figure()
for (zi, z) in enumerate(z_vals)
    plt.plot(r_vals, c_frac[:,zi],label="z = $z")
end
plt.xlabel(L"r")
plt.ylabel(L"\rho_{\mathrm{cooperate}}")
plt.legend(loc=2)
plt.title("reputation tracked on both")
fig.tight_layout()
display(fig)


fig = plt.figure()
for (zi, z) in enumerate(z_vals)
    plt.plot(r_vals, payback_frac[:,zi],label="z = $z")
end
plt.xlabel(L"r")
plt.ylabel(L"\rho_{\mathrm{payback}}")
plt.legend(loc=2)
plt.title("reputation tracked on both")
fig.tight_layout()
display(fig)
