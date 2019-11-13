#!/usr/bin/env julia

## test_payoffs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test LatticeCredit payoff function.

using Revise
using LatticeCredit

r = 2.0
z = 0.0
v = 0.0
game = LatticeGame(r, z, v)
println("no lending")
println("payoffs for a cooperator with four cooperator neighbors")
println("there should be 5 in the pot, or 10 after synergy")
println("total payoff of 1 each")
payoffs = (get_group_payoffs(4, [4,4,4,4], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")

println("payoffs for a defector with four defector neighbors")
println("there should be 0 in the pot")
println("total payoff of 0 each")
payoffs = (get_group_payoffs(0, [0,0,0,0], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")

r = 2.0
z = 1.1
v = 0.0
game = LatticeGame(r, z, v)
println("still no lending")
println("payoffs for a cooperator with four cooperator neighbors")
println("there should be 5 in the pot, or 10 after synergy")
println("total payoff of 1 each")
payoffs = (get_group_payoffs(4, [4,4,4,4], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")

println("payoffs for a cooperator with four defector neighbors")
println("there should be 1 in the pot, or 2 after synergy")
println("total payoff of -0.6 for focal indv and 0.4 for neighbors")
payoffs = (get_group_payoffs(4, [0,0,0,0], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")

println("payoffs for a cooperator lender with four cooperator payback neighbors")
println("here interest is charged, but the focal indv has nothing to lend")
println("there should be 5 in the pot, or 10 after synergy")
println("total payoff of 1 each")
payoffs = (get_group_payoffs(7, [7,7,7,7], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")

r = 2.0
z = 2.0
v = 1.0
game = LatticeGame(r, z, v)
println("now we'll give the focal indv an extra 1 to invest or lend as they see fit")

println("payoffs for a cooperator lender with four cooperator payback neighbors")
println("focal indv lends 0.25 to each neighbor")
println("they all invest")
println("there should be 6 in the pot, or 12 after synergy")
println("everyone gets 2.4 back, or 1.4 after cooperating")
println("focal indv should get 2.9, everyone else 0.9")
payoffs = (get_group_payoffs(7, [7,7,7,7], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")

println("payoffs for a cooperator lender with four cooperator shirker neighbors")
println("focal indv lends 0.25 to each neighbor")
println("they all invest")
println("there should be 6 in the pot, or 12 after synergy")
println("everyone gets 2.4 back, or 1.4 after cooperating")
payoffs = (get_group_payoffs(7, [6,6,6,6], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")

println("payoffs for a cooperator hoarder with four defector neighbors")
println("there should be 2 in the pot, or 4 after synergy")
println("focal indv should get -0.2, everyone else 0.8")
payoffs = (get_group_payoffs(4, [0,0,0,0], game))
println("focal payoff: $(payoffs[1])")
println("neighbor payoffs: $(payoffs[2])")
