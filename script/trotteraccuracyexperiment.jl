#= Construct Trotter accuracy plots for increasing system sizes. =#

mn = [
#    m n     N
    (2,1), # 2
    (3,1), # 3
    (2,2), # 4
    (2,3), # 8
    (3,2), # 9
    (2,4), # 16
    (3,3), # 27
    (2,5), # 32
    (2,6), # 64
    (3,4), # 81
]

rundata  = false
plotdata = true


if rundata
    include("../src/utils.jl")
    include("../src/device.jl")
    include("../src/pulse.jl")
    include("../src/evolve.jl")
    include("../analysis/random.jl")
    include("../analysis/stocksolutions.jl")
    include("../analysis/experiment.jl")
    include("../analysis/experiments/randomsquarepulse.jl")
end


for (m, n) in mn

##########################################################################################
#                                       RUN DATA
if rundata
    filename = "trotteraccuracy.m=$m.n=$n"

    import .Experiments
    import .TrotterAccuracyExperiment

    expmt = TrotterAccuracyExperiment.Control(n=n, m=m)
    setup = Experiments.initialize(expmt)

    # display(preface.pulses)
    # display(preface.device)

    index = TrotterAccuracyExperiment.createindex(seed_=0:9, k_=0:15)
    open("dat/$filename.csv", "w") do io
        Experiments.collectdata(io, expmt, setup, index, writeheader=true)
    end
end

##########################################################################################
#                                       PLOT DATA
if plotdata
    filename = "trotteraccuracy.m=$m.n=$n"

    using DataFrames
    using DelimitedFiles: readdlm
    using Plots

    function Plots.scatter!(df::DataFrame, field::String, label::String)
        x = df.numsteps[df[!,field] .> 0]
        y = df[!,field][df[!,field] .> 0]
        println("$label: $(length(x))/$(size(df,1)) positive infidelities")
        scatter!(x, y, label=label, markershape=:auto, markeralpha=0.7)
    end

    # READ IN DATA
    data, header = readdlm("dat/$filename.csv", '\t', header=true)
    df = DataFrame(data, vec(header))

    # PLOT ALL POINTS
    plt = scatter(
        title  = "File: $filename",
        xlabel = "Time Steps", ylabel = "Infidelity with Exact Solution",
        xscale=:log, yscale=:log,
        legend=:bottomleft,
    )
    scatter!(df, "FI", "Initial State")
    scatter!(df, "F_Direct", "Direct")
    scatter!(df, "F_Lanczos", "Lanczos")
    scatter!(df, "F_Rotate_K", "Rotate (K)")
    scatter!(df, "F_Prediag_1K", "Prediag (1K)")
    scatter!(df, "F_Prediag_2K", "Prediag (2K)")
    scatter!(df, "F_Rotate_T", "Rotate (T)")
    scatter!(df, "F_Prediag_1T", "Prediag (1T)")
    scatter!(df, "F_Prediag_2T", "Prediag (2T)")

    # png(plt, "fig/$filename.png")
    gui(plt)
end

end # END mn LOOP
