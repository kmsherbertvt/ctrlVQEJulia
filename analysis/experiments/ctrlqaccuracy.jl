#= Trotter experiment contrasting DIRECT method from that implemented in ctrlq. =#


module CtrlQExperiment
    import PyCall: pyimport

    import ..Utils
    import ..Devices
    import ..Pulses
    import ..Evolutions

    import ..Experiments
    import ..StockSolutions

    import Random: MersenneTwister
    import ..RandomConstructs

    struct Control <: Experiments.Control
        # PULSE PARAMETERS
        T::Float64              # PULSE DURATION
        ν::Float64              # PULSE FREQUENCY (ALL QUBITS)
        Ω::Float64              # PULSE AMPLITUDE (ALL QUBITS)
        # DEVICE PARAMETERS
        deviceseed::Integer     # RANDOM SEED USED TO GENERATE DEVICE
    end

    struct Setup <: Experiments.Setup
        ψI::Vector{ComplexF64}                          # INITIAL STATE
        pulses::AbstractVector{<:Pulses.PulseTemplate}  # JULIA PULSES
        device::Devices.Device                          # JULIA DEVICE
        HD                                              # DEVICE HAMILTONIAN
        ΛD                                              # DEVICE EIGENVALUES
        UD                                              # DEVICE EIGENBASIS
        a_                                              # QUBIT OPERATORS
        pobjc                                           # ctrlq PULSES
        myham                                           # ctrlq DEVICE
        cvqe_evolve                                     # ctrlq EVOLVE METHOD INTERFACE
    end

    function Experiments.initialize(expmt::Control)
        n = 4
        m = 3
        pulses = [Pulses.BasicSquarePulse(expmt.T, expmt.ν, expmt.Ω) for q in 1:n]
        device = RandomConstructs.transmondevice(MersenneTwister(expmt.deviceseed),
            n, 5.0, 0.3, 0.02,
        ) # PARAMETERS SELECTED TO QUALITATIVELY MATCH ctrlq DEFAULTS

        # PRE-DIAGONALIZE STATIC HAMILTONIAN
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
        a_ = Utils.algebra(n, m, basis=UD)          # LIST OF ROTATED ANNIHILATION OPERATORS

        # CONSTRUCT ctrlq PULSE OBJECT
        ctrlq_solve = pyimport("ctrlq.lib.solve")
        pobjc = ctrlq_solve.pulsec(
            [pulse.amplitudes for pulse ∈ pulses],
            [pulse.steptimes  for pulse ∈ pulses],
            [pulse.frequency  for pulse ∈ pulses],
            expmt.T, length(pulses), length(pulses[1].amplitudes),
        )

        # CONSTUCT ctrlq TRANSMON OBJECT
        cvqe = pyimport("ctrlq.cvqe")
        myham = cvqe.transmon(
            nqubit=length(device),
            nstate=m,
            mham=zeros(1,1),
            Hstatic=Matrix(HD),
        )

        # GRAB THE INITIAL STATE FROM ctrlq DEFAULT
        ψI = myham.initial_state[:,1]

        # PYOBJECT TO CALL ctrlq EVOLVE
        cvqe_evolve = pyimport("ctrlq.cvqe.evolve")
        return Setup(ψI, pulses, device, HD, ΛD, UD, a_, pobjc, myham, cvqe_evolve)
    end

    struct Independent <: Experiments.Independent
        numsteps::Integer               # NUMBER OF TROTTER STEPS
    end

    function Experiments.mapindex(expmt::Control, i::Integer)
        return Independent(2^i)
    end

    struct Result <: Experiments.Result
        ψ0::Vector{ComplexF64}              # SOLUTION OBTAINED FROM NUMERICAL INTEGRATION
        ψ_ctrlq::Vector{ComplexF64}         # SOLUTION OBTAINED FROM CtrlQ code
        ψ_julia::Vector{ComplexF64}         # SOLUTION OBTAINED FROM Julia code
    end

    function Experiments.runtrial(expmt::Control, setup::Setup, xvars::Independent)
        ψ0 = Evolutions.evolve(setup.ψI, setup.pulses, setup.device, Evolutions.ODE)

        ψ_ctrlq = setup.cvqe_evolve.evolve(
            setup.ψI,
            setup.pobjc,
            setup.myham,
            solver="trotter",
            nstep=xvars.numsteps,
        )[:,1]

        ψ_julia = Evolutions.evolve(setup.ψI, setup.pulses, setup.device, Evolutions.Direct;
                    numsteps=xvars.numsteps)

        return Result(ψ0, ψ_ctrlq, ψ_julia)
    end

    struct Dependent <: Experiments.Dependent
        FI::Float64             # FIDELITY BETWEEN NUMERICAL INTEGRATION AND INITIAL STATE
        F_ctrlq::Float64        # CtrlQ FIDELITY
        F_julia::Float64        # Julia FIDELITY
    end

    function Experiments.synthesize(::Control,
            setup::Setup, xvars::Independent, result::Result)
        ψ0 = result.ψ0
        return Dependent(
            Utils.infidelity(ψ0, setup.ψI),
            Utils.infidelity(ψ0, result.ψ_ctrlq),
            Utils.infidelity(ψ0, result.ψ_julia),
        )
    end

end # END MODULE
