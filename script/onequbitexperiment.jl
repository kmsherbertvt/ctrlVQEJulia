#= Identify the number of Trotter steps required for each Trotterized method
    to reach the accuracy obtained by numerical integration.

This tests a single square pulse applied to a two-level system.
=#

##########################################################################################
#                                       RUN DATA

# include("../src/utils.jl")
# include("../src/device.jl")
# include("../src/pulse.jl")
# include("../src/evolve.jl")
# include("../analysis/random.jl")
# include("../analysis/stocksolutions.jl")
# include("../analysis/experiment.jl")
# include("../analysis/experiments/onequbit.jl")
#
# import .Experiments
#
# # 2-LEVEL SYSTEM
# import .OneQubitSquarePulseExperiment
# expmt = OneQubitSquarePulseExperiment.Control(10.0, 6.0, .02, 4.4)
# setup = Experiments.initialize(expmt)
# index = OneQubitSquarePulseExperiment.createindex(seed_=0:9, k_=0:15)
# open("dat/onequbitsquarepulse.csv", "w") do io
#     Experiments.collectdata(io, expmt, setup, index, writeheader=true)
# end

# # 3-LEVEL SYSTEM
# import .OneQutritSquarePulseExperiment
# expmt = OneQutritSquarePulseExperiment.Control(10.0, 6.0, .02, 4.4, 0.3)
# setup = Experiments.initialize(expmt)
# index = OneQutritSquarePulseExperiment.createindex(seed_=0:9, k_=0:15)
# open("dat/onequtritsquarepulse.csv", "w") do io
#     Experiments.collectdata(io, expmt, setup, index, writeheader=true)
# end

##########################################################################################
#                                       PLOT DATA

# filehandle = "onequbitsquarepulse"
filehandle = "onequtritsquarepulse"

using DataFrames
using DelimitedFiles: readdlm
using Plots

# READ IN DATA
data, header = readdlm("dat/$filehandle.csv", '\t', header=true)
df = DataFrame(data, vec(header))


# PLOT ALL POINTS
function Plots.scatter!(field::String, label::String)
    x = df.numsteps[df[!,field] .> 0]
    y = df[!,field][df[!,field] .> 0]
    scatter!(x, y, label=label, markershape=:auto, markeralpha=0.7)
end

plt = scatter(
    xlabel = "Time Steps",
    ylabel = "Infidelity with Exact Solution",
    xscale=:log, yscale=:log
)
scatter!("FI", "Initial State")
scatter!("F_ODE", "ODE")
scatter!("F_Direct", "Direct")
scatter!("F_Lanczos", "Lanczos")
scatter!("F_Rotate", "Rotate")
scatter!("F_Prediag_1", "Prediag 1")
scatter!("F_Prediag_2", "Prediag 2")


png(plt, "fig/$filehandle.png")
gui(plt)
