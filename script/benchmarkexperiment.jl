#= Construct Trotter accuracy plots for increasing system sizes. =#

rundata  = true
plotdata = true

filename = "benchmark.11_18_2022.external"
modes = [
    0,  # ODE
    1,  # DIRECT
    2,  # LANCZOS
    3,  # ROTATE (Kronec)
    # 4,  # ROTATE (Tensor)
    # 5,  # PREDIAG (Kronec, order=1)
    # 7,  # PREDIAG (Kronec, order=2)
    # 6,  # PREDIAG (Tensor, order=1)
    # 8,  # PREDIAG (Tensor, order=2)
]

m_ = [2,3]

##########################################################################################
#                                       RUN DATA
if rundata
    include("../src/utils.jl")
    include("../src/device.jl")
    include("../src/pulse.jl")
    include("../src/evolve.jl")
    include("../analysis/random.jl")
    include("../analysis/stocksolutions.jl")
    include("../analysis/experiment.jl")
    include("../analysis/experiments/randomsquarepulse.jl")

    import .Experiments
    import .BenchmarkExperiment

    expmt = BenchmarkExperiment.Control(numsteps=30)

    for mode in modes; for m in m_
        result = nothing
        n = 1
        while n <= 24 ÷ m
            # RUN SINGLE EXPERIMENT
            index = BenchmarkExperiment.createindex(m, n, mode)
            open("dat/$filename.csv", "a") do io
                result = Experiments.collectdata(io, expmt, index)
            end

            # TERMINATE AFTER FIRST TRIAL THAT TAKES OVER A SECOND
            if minimum(result.benchmark).time > 1.0e9; break; end

            # UPDATE QUBIT COUNT
            n += 1
        end
    end; end
end

##########################################################################################
#                                       PLOT DATA
if plotdata
    using DataFrames
    using DelimitedFiles: readdlm
    using Plots

    modelabel = Dict(
        0 => "ODE",
        1 => "Direct",
        2 => "Lanczos",
        3 => "Rotate (K)",
        4 => "Rotate (T)",
        5 => "Prediag (1K)",
        6 => "Prediag (1T)",
        7 => "Prediag (2K)",
        8 => "Prediag (2T)",
    )

    function Plots.scatter!(df::DataFrame, m::Integer, mode::Integer, label::String)
        mask = (df.m .== m .&& df.mode .== mode)
        x = df.n[mask]
        y = df.time[mask] / 1.0e9   # CONVERT ns -> s
        scatter!(x, y, label=label, markershape=:auto, markeralpha=0.7)
    end

    # READ IN DATA
    data, header = readdlm("dat/$filename.csv", '\t', header=true)
    df = DataFrame(data, vec(header))

    # PLOT ALL POINTS
    plt = scatter(
        title  = "Benchmarking",
        xlabel = "Number of Qubits", ylabel = "Time (s)",
        yscale=:log, ylims=[1e-4, 1e2],
        legend=:topleft,
    )
    for mode in modes; for m in m_
        scatter!(df, m, mode, "$(modelabel[mode]) m=$m")
    end; end

    # png(plt, "fig/$filename.png")
    gui(plt)
end
