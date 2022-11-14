#= Identify the number of Trotter steps required for each Trotterized method
    to reach the accuracy obtained by numerical integration.

This tests ctrlq evolve vs mine.
=#

# ##########################################################################################
# #                                       RUN DATA
#
# include("../src/utils.jl")
# include("../src/device.jl")
# include("../src/pulse.jl")
# include("../src/evolve.jl")
# include("../analysis/random.jl")
# include("../analysis/stocksolutions.jl")
# include("../analysis/experiment.jl")
# include("../analysis/experiments/ctrlqaccuracy.jl")
#
# import .Experiments
#
# import .CtrlQExperiment
# expmt = CtrlQExperiment.Control(10.0, 5.0, 0.02, hash("device"))
#
# setup = Experiments.initialize(expmt)
#
# # VERIFY THAT ctrlq CORRECTLY DIAGONALIZES THE DEVICE
# import LinearAlgebra: diag
# @assert diag(setup.myham.dsham.toarray()) ≈ setup.ΛD
#
# # nonzeros = findall(x -> abs(x)>1e-10, ctrlqH)
# # ctr = 0
# # for ix in nonzeros
# #     (i,j) = Tuple(ix)
# #     if i == j; continue; end
# #     println("$i $j $(ctrlqH[i,j])")
# #     global ctr += 1
# # end
# # println("$ctr of-diagonal elements")
# #
# # import LinearAlgebra: diag, norm
# # ctrlΛ = diag(ctrlqH)
# # err = norm(ctrlΛ-setup.ΛD)
# # println("$err")
#
#
# open("dat/ctrlqaccuracy.csv", "a") do io
#     Experiments.collectdata(io, expmt, setup, 0:15, writeheader=true)
# end


##########################################################################################
#                                       PLOT DATA

using DataFrames
using DelimitedFiles: readdlm
using Plots

# READ IN DATA
data, header = readdlm("dat/ctrlqaccuracy.csv", '\t', header=true)
df = DataFrame(data, vec(header))

# PLOT CURVES
function Plots.scatter!(field::String, label::String)
    # x = df.numsteps[df[!,field] .> 0 .&& df.Ω .== 0.02]
    # y = df[!,field][df[!,field] .> 0 .&& df.Ω .== 0.02]
    x = df.numsteps[df.Ω .== 0.02]
    y = df[!,field][df.Ω .== 0.02]
    scatter!(x, y, label=label)
end

plt = scatter(
    title = "Time Evolution via Direct Exponentiation",
    xlabel = "Time Steps",
    ylabel = "Infidelity with ODE Solution",
    xscale=:log, yscale=:log,
    ylims=[1e-16,1], ytick=10.0 .^ (-1:-2:-16),
    legend=:bottomleft,
)
scatter!("FI", "Initial State")
scatter!("F_ctrlq", "CtrlQ")
scatter!("F_julia", "New Code")


png(plt, "fig/ctrlqaccuracy.png")
gui(plt)
