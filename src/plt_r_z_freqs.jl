#!/usr/bin/env julia

## plt_credit_freqs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from LatticeLending simulations.
## Look at dependence of type frequencies on r and z.

using CSV, PyPlot, Statistics

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

# output info about the fixation probabilities?
verbose = false

# load simulation output as a dataframe
runs = CSV.read("output/lending_r_z_3.csv")

# dicts to store fixation probabilities
type_freqs = Dict{Tuple{Float64, Float64, Float64, Bool, Bool},Array{Float64, 1}}()
c_freqs = Dict{Tuple{Float64, Float64, Float64, Bool, Bool},Float64}()
p_freqs = Dict{Tuple{Float64, Float64, Float64, Bool, Bool},Float64}()


N = sort(unique(runs[:N]))[1]

# get unique values from the runs dataframe
z_vals = sort(unique(runs[:z]))
d_vals = sort(unique(runs[:d]))
r_vals = sort(unique(runs[:r]))

reputation_coop = [false, true]
reputation_payback = [false, true]

param_combs = collect(Base.product(z_vals, d_vals, r_vals, reputation_payback, reputation_coop))

# iterate over all update types, ks, and ws

for (pi, param_comb) in enumerate(param_combs)
    z, d, r, payback, coop = param_comb
    #println("$param_comb")
    type_freqs[param_comb] = zeros(Float64, 4)
    c_freqs[param_comb] = 0.0
    p_freqs[param_comb] = 0.0
    tmp_runs = runs[(runs[:z] .== z) .& (runs[:d] .== d) .& (runs[:r] .== r) .& (runs[:reputation_payback] .== payback) .& (runs[:reputation_coop] .== coop), :]
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
        freqs = parse.(Float64,String.(split(run[:mean_freqs], r";|,| |\[|\]")[2:end-1]))
        type_freqs[param_comb] += freqs/size(tmp_runs, 1)
        c_freqs[param_comb] += sum([freqs[3], freqs[4]])/size(tmp_runs, 1)
        p_freqs[param_comb] += sum([freqs[2], freqs[4]])/size(tmp_runs, 1)
    end
end

zmin, zmax, rmin, rmax = z_vals[1], z_vals[end], r_vals[1], r_vals[end]
rskip = 1

for coop in reputation_coop
    for payback in reputation_payback
        fig, axs = plt.subplots(2, length(d_vals)-1, figsize=(12,3))
        for (di, d) in enumerate(d_vals[2:end])
            c_freq_im = zeros(ceil(Int64, length(r_vals)/rskip), length(z_vals))
            p_freq_im = zeros(ceil(Int64, length(r_vals)/rskip), length(z_vals))
            for (zi, z) in enumerate(z_vals)
                for (ri, r) in enumerate(r_vals[1:rskip:end])
                    c_freq_im[ri, zi] = c_freqs[z, d, r, payback, coop]
                    p_freq_im[ri, zi] = p_freqs[z, d, r, payback, coop]
                end
            end
            # ax = axs[1, di]
            # im = ax.imshow(c_freq_im, vmin=0, vmax=1, extent = [zmin, zmax, rmin, rmax])
            # ax.set_xlabel("z")
            # ax.set_ylabel("r")
            # ax.set_title("d = $d, coop = $coop, payback = $payback")
            # #fig.colorbar(im)
            # if di == length(d_vals)
            #     cax = divider.append_axes("right", size="5%", pad=0.05)
            #     fig.colorbar(im, cax = cax)
            # end
            ax = axs[1,di]
            im = ax.imshow(p_freq_im, origin = "lower",
                aspect=length(z_vals)*rskip/length(r_vals),
                vmin=0, vmax=1,
                extent = [zmin, zmax, rmin, rmax])
            ax.set_xlabel("z")
            ax.set_ylabel("r")
            ax.set_title("d = $d")

            ax = axs[2,di]
            im = ax.imshow(c_freq_im, origin = "lower",
                aspect=length(z_vals)*rskip/length(r_vals),
                vmin=0, vmax=1,
                extent = [zmin, zmax, rmin, rmax])
            ax.set_xlabel("z")
            ax.set_ylabel("r")
            ax.set_title("d = $d")
            #fig.colorbar(im)
        #    cbar = ax.cax.colorbar(im)
        #    cbar = grid.cbar_axes[0].colorbar(im)
        end
        fig.suptitle(tracking_label([payback, coop]))
        #fig.tight_layout()
        plt.subplots_adjust(top=0.85)
        display(fig)
    end
end

#[println("$keyt, $(mean(p_c_trajectories[(keyt)])[2])") for keyt in Base.product(k_vals, w_vals, update_types)]

color_list = ["tab:blue", "tab:orange", "tab:green", "tab:red",
    "tab:purple", "tab:gray", "tab:pink"]

# fig, axs = plt.subplots(length(z_vals), length(focal_types),
#     figsize=(14,7))#, sharey="row", sharex="col")
# fig = plt.figure()
# for (zi, z) in enumerate(z_vals)
#     #for (fi, focal_type) in enumerate(focal_types)
#         #ax=axs[zi, fi]
#         for (vi, v) in enumerate(v_vals)
#             plot(r_vals, [c_freqs[z, v, r] for r in r_vals], label="v = $v")
#         end
#     #end
# end
# plt.ylabel(L"\rho_c")
# plt.xlabel(L"r")
# plt.legend(loc=2)
# plt.title("all strategies")
# plt.tight_layout()
# display(fig)

# fig, axs = plt.subplots(1, length(d_vals)-1, figsize = (20, 8))
# for (di, d) in enumerate(d_vals[2:end])
#     ax = axs[di]
#     ax.plot([zmin, zmax], [rmin, rmax]*d)
#     ax.set_xlim([zmin, zmax])
#     ax.set_xlabel("z")
#     ax.set_ylim([rmin, rmax])
#     ax.set_ylabel("r")
# end
# display(fig)
