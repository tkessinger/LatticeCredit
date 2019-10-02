#!/usr/bin/env julia

## plt_credit_freqs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from LatticeCredit simulations.

using CSV, PyPlot, Statistics

# output info about the fixation probabilities?
verbose = false

# load simulation output as a dataframe
runs = CSV.read("output/test_redo_20_3.csv")

# dicts to store fixation probabilities
type_freqs = Dict{Tuple{Float64, Float64, Float64, Int64},Array{Float64, 1}}()

c_freqs = Dict{Tuple{Float64, Float64, Float64},Float64}()

N = sort(unique(runs[:N]))[1]

# get unique values from the runs dataframe
z_vals = sort(unique(runs[:z]))
v_vals = sort(unique(runs[:v]))
r_vals = sort(unique(runs[:r]))
focal_types = sort(unique(runs[:focal_type]))

param_combs = collect(Base.product(z_vals, v_vals, r_vals, focal_types))

# iterate over all update types, ks, and ws

for (pi, param_comb) in enumerate(param_combs)
    z, v, r, focal_type = param_comb
    c_param_comb = (z, v, r)
    #println("$param_comb")
    type_freqs[param_comb] = zeros(Float64, 8)
    c_freqs[c_param_comb] = 0
    tmp_runs = runs[(runs[:z] .== z) .& (runs[:v] .== v) .& (runs[:r] .== r) .& (runs[:focal_type] .== focal_type), :]
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
        freqs = parse.(Float64,String.(split(run[:mean_freqs], r";|,|\[|\]")[2:end-1]))
        type_freqs[param_comb] += freqs/size(tmp_runs, 1)
        c_freqs[c_param_comb] += sum(freqs[5:end])/size(tmp_runs, 1)
    end
end

#[println("$keyt, $(mean(p_c_trajectories[(keyt)])[2])") for keyt in Base.product(k_vals, w_vals, update_types)]

color_list = ["tab:blue", "tab:orange", "tab:green", "tab:red",
    "tab:purple", "tab:gray", "tab:pink"]

# fig, axs = plt.subplots(length(z_vals), length(focal_types),
#     figsize=(14,7))#, sharey="row", sharex="col")
fig = plt.figure()
for (zi, z) in enumerate(z_vals)
    #for (fi, focal_type) in enumerate(focal_types)
        #ax=axs[zi, fi]
        for (vi, v) in enumerate(v_vals)
            plot(r_vals, [c_freqs[z, v, r] for r in r_vals], label="v = $v")
        end
    #end
end
plt.ylabel(L"\rho_c")
plt.xlabel(L"r")
plt.legend(loc=2)
plt.title("all strategies")
plt.tight_layout()
display(fig)
