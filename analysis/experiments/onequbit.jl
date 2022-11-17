#= Identify the number of Trotter steps required for each Trotterized method
    to reach the accuracy obtained by numerical integration.

This file contains experiments for a single-body system that can be solved exactly.
=#

module OneQubitSquarePulseExperiment
    import ..Utils
    import ..Devices
    import ..Pulses
    import ..Evolutions

    import ..Experiments
    import ..StockSolutions

    import Random: MersenneTwister
    import ..RandomConstructs

    struct Control <: Experiments.Control
        T::Float64                  # PULSE DURATION
        ν::Float64                  # PULSE FREQUENCY
        Ω::Float64                  # PULSE AMPLITUDE
        ω::Float64                  # QUBIT RESONANCE FREQUENCY
    end

    struct Setup <: Experiments.Setup
        pulses::AbstractVector{<:Pulses.PulseTemplate}  # SINGLE-QUBIT PULSE
        device::Devices.Device                          # SINGLE-QUBIT DEVICE
    end

    function Experiments.initialize(expmt::Control)
        return Setup(
            [Pulses.BasicSquarePulse(expmt.T, expmt.ν, expmt.Ω)],
            Devices.Transmon([expmt.ω], [0]),
        )
    end

    struct Independent <: Experiments.Independent
        seed::Integer                   # RANDOM SEED USED TO GENERATE INITIAL STATE
        numsteps::Int                   # NUMBER OF TROTTER STEPS
    end

    function Experiments.mapindex(expmt::Control, i::Integer)
        seed, k = divrem(i, 100)
        return Independent(seed, 2^k)
    end

    struct Result <: Experiments.Result
        ψI::Vector{ComplexF64}          # INITIAL STATE
        ψ0::Vector{ComplexF64}          # ANALYTICAL SOLUTION
        ψ_ODE::Vector{ComplexF64}       # NUMERICALLY CALCULATED SOLUTIONS
        ψ_Direct::Vector{ComplexF64}            # ⋮
        ψ_Lanczos::Vector{ComplexF64}           # ⋮
        ψ_Rotate::Vector{ComplexF64}            # ⋮
        ψ_Prediag_1::Vector{ComplexF64}         # ⋮
        ψ_Prediag_2::Vector{ComplexF64}         # ⋮
    end

    function Experiments.runtrial(expmt::Control, setup::Setup, xvars::Independent)
        rng = MersenneTwister(xvars.seed)
        ψI = RandomConstructs.statevector(rng, 2)
        return Result(ψI,
            StockSolutions.onequbitsquarepulse(ψI, expmt.T, expmt.ν, expmt.Ω, expmt.ω),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.ODE),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Direct;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Lanczos;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Rotate;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=1),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=2),
        )
    end

    struct Dependent <: Experiments.Dependent
        FI::Float64             # FIDELITY   BETWEEN EXACT SOLUTION AND INITIAL STATE
        F_ODE::Float64          # FIDELITIES BETWEEN EXACT AND SIMULATED SOLUTIONS
        F_Direct::Float64                           #  ⋮
        F_Lanczos::Float64                          #  ⋮
        F_Rotate::Float64                           #  ⋮
        F_Prediag_1::Float64                        #  ⋮
        F_Prediag_2::Float64                        #  ⋮
    end

    function Experiments.synthesize(::Control, ::Setup, xvars::Independent, result::Result)
        ψ0 = result.ψ0
        return Dependent(
            Utils.infidelity(ψ0, result.ψI),
            Utils.infidelity(ψ0, result.ψ_ODE),
            Utils.infidelity(ψ0, result.ψ_Direct),
            Utils.infidelity(ψ0, result.ψ_Lanczos),
            Utils.infidelity(ψ0, result.ψ_Rotate),
            Utils.infidelity(ψ0, result.ψ_Prediag_1),
            Utils.infidelity(ψ0, result.ψ_Prediag_2),
        )
    end

    """ Convenience function to produce a list of indices,
            given specified lists of independent variables.
    """
    function createindex(;seed_=nothing, k_=nothing)
        return (100*seed + k for (seed, k) in Iterators.product(seed_, k_))
    end

end # END MODULE



module OneQutritSquarePulseExperiment
    import ..Utils
    import ..Devices
    import ..Pulses
    import ..Evolutions

    import ..Experiments
    import ..StockSolutions

    import Random: MersenneTwister
    import ..RandomConstructs

    struct Control <: Experiments.Control
        T::Float64                  # PULSE DURATION
        ν::Float64                  # PULSE FREQUENCY
        Ω::Float64                  # PULSE AMPLITUDE
        ω::Float64                  # QUBIT RESONANCE FREQUENCY
        δ::Float64                  # QUBIT ANHARMONICITY
    end

    struct Setup <: Experiments.Setup
        pulses::AbstractVector{<:Pulses.PulseTemplate}  # SINGLE-QUBIT PULSE
        device::Devices.Device                          # SINGLE-QUBIT DEVICE
    end

    function Experiments.initialize(expmt::Control)
        return Setup(
            [Pulses.BasicSquarePulse(expmt.T, expmt.ν, expmt.Ω)],
            Devices.Transmon([expmt.ω], [expmt.δ]),
        )
    end

    struct Independent <: Experiments.Independent
        seed::Integer                   # RANDOM SEED USED TO GENERATE INITIAL STATE
        numsteps::Int                   # NUMBER OF TROTTER STEPS
    end

    function Experiments.mapindex(expmt::Control, i::Integer)
        seed, k = divrem(i, 100)
        return Independent(seed, 2^k)
    end

    struct Result <: Experiments.Result
        ψI::Vector{ComplexF64}          # INITIAL STATE
        ψ0::Vector{ComplexF64}          # ANALYTICAL SOLUTION
        ψ_ODE::Vector{ComplexF64}       # NUMERICALLY CALCULATED SOLUTIONS
        ψ_Direct::Vector{ComplexF64}            # ⋮
        ψ_Lanczos::Vector{ComplexF64}           # ⋮
        ψ_Rotate::Vector{ComplexF64}            # ⋮
        ψ_Prediag_1::Vector{ComplexF64}         # ⋮
        ψ_Prediag_2::Vector{ComplexF64}         # ⋮
    end

    function Experiments.runtrial(expmt::Control, setup::Setup, xvars::Independent)
        rng = MersenneTwister(xvars.seed)
        ψI = RandomConstructs.statevector(rng, 3)
        return Result(ψI,
            StockSolutions.onequtritsquarepulse(ψI,
                    expmt.T, expmt.ν, expmt.Ω, expmt.ω, expmt.δ),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.ODE),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Direct;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Lanczos;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Rotate;
                    numsteps=xvars.numsteps),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=1),
            Evolutions.evolve(ψI, setup.pulses, setup.device, Evolutions.Prediag;
                    numsteps=xvars.numsteps, suzukiorder=2),
        )
    end

    struct Dependent <: Experiments.Dependent
        FI::Float64             # FIDELITY   BETWEEN EXACT SOLUTION AND INITIAL STATE
        F_ODE::Float64          # FIDELITIES BETWEEN EXACT AND SIMULATED SOLUTIONS
        F_Direct::Float64                           #  ⋮
        F_Lanczos::Float64                          #  ⋮
        F_Rotate::Float64                           #  ⋮
        F_Prediag_1::Float64                        #  ⋮
        F_Prediag_2::Float64                        #  ⋮
    end

    function Experiments.synthesize(::Control, ::Setup, xvars::Independent, result::Result)
        ψ0 = result.ψ0
        return Dependent(
            Utils.infidelity(ψ0, result.ψI),
            Utils.infidelity(ψ0, result.ψ_ODE),
            Utils.infidelity(ψ0, result.ψ_Direct),
            Utils.infidelity(ψ0, result.ψ_Lanczos),
            Utils.infidelity(ψ0, result.ψ_Rotate),
            Utils.infidelity(ψ0, result.ψ_Prediag_1),
            Utils.infidelity(ψ0, result.ψ_Prediag_2),
        )
    end

    """ Convenience function to produce a list of indices,
            given specified lists of independent variables.
    """
    function createindex(;seed_=nothing, k_=nothing)
        return (100*seed + k for (seed, k) in Iterators.product(seed_, k_))
    end

end # END MODULE
