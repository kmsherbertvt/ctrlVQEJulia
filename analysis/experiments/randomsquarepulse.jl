#= Experiments using randomized initial state/square pulse/device
    to assess accuracy and time of evolution methods, scaled for increasing system size.
=#


""" Calculate accuracy of Trotterized evolution methods as number of time steps increases.

Accuracy is measured by overlap with numerically-integrated solution.
Note that, at some point, at least some Trotterized methods
    are likely to be *more* accurate than the numerically-integrated solution,
    so be careful how you interpret results!

See `Experiments` module for general documentation.

"""
module TrotterAccuracyExperiment
    import ..Utils
    import ..Devices
    import ..Pulses
    import ..Evolutions

    import ..Experiments
    import ..StockSolutions

    import Random: MersenneTwister
    import ..RandomConstructs

    struct Control <: Experiments.Control
        n::Integer              # NUMBER OF QUBITS
        m::Integer              # NUMBER OF LEVELS PER QUBIT
        # PULSE PARAMETERS
        pulseseed::Integer      # RANDOM SEED USED TO GENERATE PULSE
        T::Real                 # PULSE DURATION            (*NOT* RANDOMIZED)
        ν::Real                 # MEAN PULSE FREQUENCY
        Ω::Real                 # MAX PULSE AMPLITUDE
        W::Integer              # MEAN # OF WINDOWS
        σν::Union{Real,Nothing} # √(σ²) OF PULSE FREQUENCY
        σΩ::Union{Real,Nothing} # RANDOMIZE AMPLITUDE?
        σW::Union{Real,Nothing} # RANDOMIZE # OF WINDOWS?
        # DEVICE PARAMETERS
        deviceseed::Integer     # RANDOM SEED USED TO GENERATE DEVICE
        ω::Real                 # MEAN RESONANCE FREQUENCY
        δ::Real                 # MEAN ANHARMONICITY
        g::Real                 # MEAN COUPLING STRENGTH
        σC::Union{Real,Nothing} # PAIR-WISE CHANCE OF COUPLING (beyond linear)
        σω::Union{Real,Nothing} # √(σ²) OF RESONANT FREQUENCIES
        σδ::Union{Real,Nothing} # √(σ²) OF ANHARMONICITY
        σg::Union{Real,Nothing} # √(σ²) OF COUPLING STRENGTH
    end

    """
        Control()

    Default constructor: pass any control variables as keyword arguments

    Mean pulse/device values default to typical IBM-like numbers;
        variance values default to 1/9 of the mean
        (resulting in Gamma distribution with shape α=3).
    """
    Control(; n=1, m=2,
        pulseseed=hash("pulse"), T=10.0, ν=5.0, Ω=0.02, W=1, σν=ν/9, σΩ=true, σW=nothing,
        deviceseed=hash("device"), ω=5.0, δ=0.3, g=0.02, σC=nothing, σω=ω/9, σδ=δ/9, σg=g/9,
    ) = Control(n, m,
        pulseseed, T, ν, Ω, W, σν, σΩ, σW,
        deviceseed, ω, δ, g, σC, σω, σδ, σg,
    )

    struct Setup <: Experiments.Setup
        pulses::AbstractVector{<:Pulses.PulseTemplate}  # SINGLE-QUBIT PULSE
        device::Devices.Device                          # SINGLE-QUBIT DEVICE
    end

    function Experiments.initialize(expmt::Control)
        pulserng = MersenneTwister(expmt.pulseseed)
        pulses = [
            RandomConstructs.squarepulse(pulserng,
                expmt.T, expmt.ν, expmt.Ω, expmt.W;
                σν=expmt.σν, σΩ=expmt.σΩ, σW=expmt.σW
            ) for q in 1:expmt.n
        ]
        device = RandomConstructs.transmondevice(MersenneTwister(expmt.deviceseed),
            expmt.n, expmt.ω, expmt.δ, expmt.g;
            σC=expmt.σC, σω=expmt.σω, σδ=expmt.σδ, σg=expmt.σg
        )
        return Setup(pulses, device)
    end

    struct Independent <: Experiments.Independent
        seed::Integer                   # RANDOM SEED USED TO GENERATE INITIAL STATE
        numsteps::Integer               # NUMBER OF TROTTER STEPS
    end

    function Experiments.mapindex(expmt::Control, i::Integer)
        seed, k = divrem(i, 100)
        return Independent(seed, 2^k)
    end

    struct Result <: Experiments.Result
        ψI::AbstractVector{<:Number}        # INITIAL STATE
        ψ0::AbstractVector{<:Number}        # ANALYTICAL SOLUTION
        ψ_Direct::AbstractVector{<:Number}  # NUMERICALLY CALCULATED SOLUTIONS
        ψ_Lanczos::AbstractVector{<:Number}             # ⋮
        ψ_Rotate_K::AbstractVector{<:Number}            # ⋮
        ψ_Prediag_1K::AbstractVector{<:Number}          # ⋮
        ψ_Prediag_2K::AbstractVector{<:Number}          # ⋮
        ψ_Rotate_T::AbstractVector{<:Number}            # ⋮
        ψ_Prediag_1T::AbstractVector{<:Number}          # ⋮
        ψ_Prediag_2T::AbstractVector{<:Number}          # ⋮
    end

    function Experiments.runtrial(expmt::Control, setup::Setup, xvars::Independent)
        ψI = RandomConstructs.statevector(MersenneTwister(xvars.seed), expmt.m^expmt.n)
        return Result(ψI,
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.ODE),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Direct;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Lanczos;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Rotate;
                    numsteps=xvars.numsteps,
                    qubitapplymode=Evolutions.Kronec()),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=1,
                    qubitapplymode=Evolutions.Kronec()),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=2,
                    qubitapplymode=Evolutions.Kronec()),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Rotate;
                    numsteps=xvars.numsteps,
                    qubitapplymode=Evolutions.Tensor()),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=1,
                    qubitapplymode=Evolutions.Tensor()),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=2,
                    qubitapplymode=Evolutions.Tensor()),
        )
    end

    struct Dependent <: Experiments.Dependent
        FI::Real                # FIDELITY   BETWEEN BLACK-BOX SOLUTION AND INITIAL STATE
        F_Direct::Real          # FIDELITIES BETWEEN BLACK-BOX AND SIMULATED SOLUTIONS
        F_Lanczos::Real                             #  ⋮
        F_Rotate_K::Real                            #  ⋮
        F_Prediag_1K::Real                          #  ⋮
        F_Prediag_2K::Real                          #  ⋮
        F_Rotate_T::Real                            #  ⋮
        F_Prediag_1T::Real                          #  ⋮
        F_Prediag_2T::Real                          #  ⋮
    end

    function Experiments.synthesize(::Control, ::Setup, xvars::Independent, result::Result)
        ψ0 = result.ψ0
        return Dependent(
            Utils.infidelity(ψ0, result.ψI),
            Utils.infidelity(ψ0, result.ψ_Direct),
            Utils.infidelity(ψ0, result.ψ_Lanczos),
            Utils.infidelity(ψ0, result.ψ_Rotate_K),
            Utils.infidelity(ψ0, result.ψ_Prediag_1K),
            Utils.infidelity(ψ0, result.ψ_Prediag_2K),
            Utils.infidelity(ψ0, result.ψ_Rotate_T),
            Utils.infidelity(ψ0, result.ψ_Prediag_1T),
            Utils.infidelity(ψ0, result.ψ_Prediag_2T),
        )
    end

    """
    Convenience function to produce a list of indices,
        given specified lists of independent variables.
    """
    function createindex(;seed_=nothing, k_=nothing)
        return (100*seed + k for (seed, k) in Iterators.product(seed_, k_))
    end

end # END MODULE



""" Benchmark performance (time and memory) of each evolution method.

See `Experiments` module for general documentation.

"""
module BenchmarkExperiment
    using BenchmarkTools

    import ..Utils
    import ..Devices
    import ..Pulses
    import ..Evolutions

    import ..Experiments
    import ..StockSolutions

    import Random: MersenneTwister
    import ..RandomConstructs

    struct Control <: Experiments.Control
        numsteps::Integer       # NUMBER OF TROTTER STEPS FOR TIME EVOLUTION METHODS
        # STATE PARAMETERS
        stateseed::Integer      # RANDOM SEED USED TO GENERATE INITIAL STATE
        # PULSE PARAMETERS
        pulseseed::Integer      # RANDOM SEED USED TO GENERATE PULSE
        T::Real                 # PULSE DURATION            (*NOT* RANDOMIZED)
        ν::Real                 # MEAN PULSE FREQUENCY
        Ω::Real                 # MAX PULSE AMPLITUDE
        W::Integer              # MEAN # OF WINDOWS
        σν::Union{Real,Nothing} # √(σ²) OF PULSE FREQUENCY
        σΩ::Union{Real,Nothing} # RANDOMIZE AMPLITUDE?
        σW::Union{Real,Nothing} # RANDOMIZE # OF WINDOWS?
        # DEVICE PARAMETERS
        deviceseed::Integer     # RANDOM SEED USED TO GENERATE DEVICE
        ω::Real                 # MEAN RESONANCE FREQUENCY
        δ::Real                 # MEAN ANHARMONICITY
        g::Real                 # MEAN COUPLING STRENGTH
        σC::Union{Real,Nothing} # PAIR-WISE CHANCE OF COUPLING (beyond linear)
        σω::Union{Real,Nothing} # √(σ²) OF RESONANT FREQUENCIES
        σδ::Union{Real,Nothing} # √(σ²) OF ANHARMONICITY
        σg::Union{Real,Nothing} # √(σ²) OF COUPLING STRENGTH
    end

    """
        Control()

    Default constructor: pass any control variables as keyword arguments

    Mean pulse/device values default to typical IBM-like numbers;
        variance values default to 1/9 of the mean
        (resulting in Gamma distribution with shape α=3).
    """
    Control(; numsteps=100, stateseed=hash("state"),
        pulseseed=hash("pulse"), T=10.0, ν=5.0, Ω=0.02, W=1, σν=ν/9, σΩ=true, σW=nothing,
        deviceseed=hash("device"), ω=5.0, δ=0.3, g=0.02, σC=nothing, σω=ω/9, σδ=δ/9, σg=g/9,
    ) = Control(numsteps, stateseed,
        pulseseed, T, ν, Ω, W, σν, σΩ, σW,
        deviceseed, ω, δ, g, σC, σω, σδ, σg,
    )

    struct Setup <: Experiments.Setup
    end

    function Experiments.initialize(expmt::Control)
        return Setup()
    end

    struct Independent <: Experiments.Independent
        m::Integer                          # NUMBER OF STATES PER QUBIT
        n::Integer                          # NUMBER OF QUBITS
        mode::Integer                       # EVOLUTION MODE
    end

    function Experiments.mapindex(expmt::Control, i::Integer)
        i, m = divrem(i, 10)
        mode, n = divrem(i, 100)
        return Independent(m, n, mode)
    end

    struct Result <: Experiments.Result
        pulses::AbstractVector{<:Pulses.PulseTemplate}  # SINGLE-QUBIT PULSE
        device::Devices.Device                          # SINGLE-QUBIT DEVICE
        ψI::AbstractVector{<:Number}                    # INITIAL STATE
        benchmark::BenchmarkTools.Trial                 # BENCHMARK TRIAL
    end

    function Experiments.runtrial(expmt::Control, setup::Setup, xvars::Independent)
        m, n, mode = xvars.m, xvars.n, xvars.mode

        # SET UP CONTROL OBJECTS
        pulserng = MersenneTwister(expmt.pulseseed)
        pulses = [
            RandomConstructs.squarepulse(pulserng,
                expmt.T, expmt.ν, expmt.Ω, expmt.W;
                σν=expmt.σν, σΩ=expmt.σΩ, σW=expmt.σW
            ) for q in 1:n
        ]
        device = RandomConstructs.transmondevice(MersenneTwister(expmt.deviceseed),
            n, expmt.ω, expmt.δ, expmt.g;
            σC=expmt.σC, σω=expmt.σω, σδ=expmt.σδ, σg=expmt.σg
        )
        ψI = RandomConstructs.statevector(MersenneTwister(expmt.stateseed), m^n)


        # PRE-DIGAONALIZE DEVICE HAMILTONIAN
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS


        # DIFFERENT EVOLUTION METHOD FOR EACH MODE
        if xvars.mode == 0      # ODE
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.ODE;
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 1  # DIRECT
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Direct; numsteps=$expmt.numsteps,
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 2  # LANCZOS
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Lanczos; numsteps=$expmt.numsteps,
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 3  # ROTATE: Kronec
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Rotate; numsteps=$expmt.numsteps,
                    qubitapplymode=Evolutions.Kronec(),
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 4  # ROTATE: Tensor
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Rotate; numsteps=$expmt.numsteps,
                    qubitapplymode=Evolutions.Tensor(),
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 5  # PREDIAG: Kronec, order=1
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Prediag; numsteps=$expmt.numsteps,
                    qubitapplymode=Evolutions.Kronec(), suzukiorder=1,
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 6  # PREDIAG: Tensor, order=1
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Prediag; numsteps=$expmt.numsteps,
                    qubitapplymode=Evolutions.Tensor(), suzukiorder=1,
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 7  # PREDIAG: Kronec, order=2
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Prediag; numsteps=$expmt.numsteps,
                    qubitapplymode=Evolutions.Kronec(), suzukiorder=1,
                    ΛD=$ΛD, UD=$UD)
        elseif xvars.mode == 8  # PREDIAG: Tensor, order=2
            benchmark = @benchmark Evolutions.evolve($ψI, $pulses, $device,
                    Evolutions.Prediag; numsteps=$expmt.numsteps,
                    qubitapplymode=Evolutions.Tensor(), suzukiorder=2,
                    ΛD=$ΛD, UD=$UD)
        else
            errror("Unsupported mode $mode")
        end

        return Result(pulses, device, ψI, benchmark)
    end

    struct Dependent <: Experiments.Dependent
        time::Real          # MINIMUM EXECUTION TIME (ns)
        gctime::Real        # MINIMUM COMPILATION TIME (ns)
        memory::Integer     # MINIMUM MEMORY CONSUMPTION (bytes)
        allocs::Integer     # MINIMUM NUMBER OF ALLOCATIONS
    end

    function Experiments.synthesize(::Control, ::Setup, xvars::Independent, result::Result)
        mintrial = minimum(result.benchmark)
        return Dependent(
            mintrial.time,
            mintrial.gctime,
            mintrial.memory,
            mintrial.allocs,
        )
    end

    """ Convenience function to produce an index from a set of independent variables. """
    function createindex(m::Integer, n::Integer, mode::Integer)
        return m + 10*n + 1000*mode
    end

end # END MODULE
