#= Code to evaluate the gradient-per-time-step of pulses in time. =#

module Gradients

import LinearAlgebra: eigen, Hermitian, Diagonal, norm, lmul!, rmul!, mul!
import TensorOperations: ncon

import ..Utils
import ..Pulses
import ..Devices
import ..Evolutions: IOBasisMode, QubitApplyMode, Kronec, Tensor, DeviceBasis, QubitBasis

import ..Evolutions


"""
    gradientsignal(
        ψI::Vector{ComplexF64},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        observable::Hermitian;
        iobasis::IOBasisMode = DeviceBasis(),
        r::Integer = 2000,
        qubitapplymode::QubitApplyMode = Kronec(),

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ψ),                      # SIZE OF STATEVECTOR
        n = length(device),                 # NUMBER OF QUBITS
        m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
        T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
        t_= range(0,T,r+1),                 # TIME GRID
        Δt= T / r,                          # DURATION OF EACH TIME STEP

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        ΛD = nothing,                       # EIGENVALUES OF STATIC HAMILTONIAN
        UD = nothing,                       # CORRESPONDING EIGENVECTORS
        V  = nothing,                       # REPEATED DEVICE ACTION
        a  = nothing,                       # SINGLE-QUBIT ANNIHILATION OPERATOR

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        tmpV = Vector{ComplexF64}(undef, N),    # FOR MATRIX-VECTOR MULTIPLICATION
        tmpM_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n],   # QUBIT-WISE DRIVE OPERATORS
        tmpK_ = nothing,                        # FOR APPLYING OPERATORS
                                                    # (default depends on `qubitapplymode`)
        ∂Ω = Vector{ComplexF64}(undef, r+1),    # STORES GRADIENT FUNCTION
    )

TODO: Formally document.

This method calculates the gradient of the expectation value of an observable Ô,
    with respect to the pulse amplitudes at each Trotterized time step.
One may use the chain rule to use this gradient to calculate the gradient
    with respect to arbitrary parameters that form an overall pulse shape.

The obserable Ô should be given as an N×N Hermitian matrix,
    ie. a matrix cast to the `Hermitian` type by `LinearAlgebra.Hermitian(matrix)`.

All other parameters are as found in `Evolutions.evolve!` using the `Rotate` method.

"""
function gradientsignal(
    ψI::Vector{ComplexF64},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    observable::Hermitian;
    iobasis::IOBasisMode = DeviceBasis(),
    r::Integer = 2000,
    qubitapplymode::QubitApplyMode = Kronec(),

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψI),                     # SIZE OF STATEVECTOR
    n = length(device),                 # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,r+1),                 # TIME GRID
    τ = T / r,                          # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    ΛD = nothing,                       # EIGENVALUES OF STATIC HAMILTONIAN
    UD = nothing,                       # CORRESPONDING EIGENVECTORS
    V  = nothing,                       # REPEATED DEVICE ACTION
    a  = nothing,                       # SINGLE-QUBIT ANNIHILATION OPERATOR

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    tmpV = Vector{ComplexF64}(undef, N),    # FOR MATRIX-VECTOR MULTIPLICATION
    tmpM_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n],   # QUBIT-WISE DRIVE OPERATORS
    tmpK_ = nothing,                        # FOR APPLYING OPERATORS
                                                # (default depends on `qubitapplymode`)
    ∂Ω = Matrix{Float64}(undef, r+1, n), # STORES GRADIENT FUNCTION
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS
    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end; if V === nothing
        V = UD* Diagonal(exp.((-im*τ) * ΛD)) *UD'  # REPEATED DEVICE ACTION
    end; if a === nothing
        a = Utils.a_matrix(m)                       # SINGLE-QUBIT ANNIHILATION OPERATOR
    end; if tmpK_ === nothing
        tmpK_ = (
            qubitapplymode isa Kronec ?
                [Matrix{ComplexF64}(undef, m^q, m^q) for q ∈ 1:n] :
            qubitapplymode isa Tensor ? [
                Dims(m for _ in 1:n),                   # RESHAPING DIMENSIONS
                [[[-q, q] for q in 1:n]..., n:-1:1],    # TENSOR INDICES
                -n:-1,                                  # OUTPUT PERMUTATION
                zeros(Bool, n+1),                       # ADJOINT FLAG
            ] : error("Invalid `QubitApplyMode` object. (How did you manage that???)")
        )
    end

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa DeviceBasis; ψI = mul!(tmpV, UD, ψI); end

    ######################################################################################
    #                                 INITIALIZATION

    ##########################
    # INITIALIZE ψ AS |ψI⟩
    ##########################
    ψ = copy(ψI)

    ##########################
    # INITIALIZE σ AS U†HU|ψI⟩
    ##########################
    σ = copy(ψI)

    # TIME-EVOLUTION UNDER TRAPEZOIDAL RULE
    σ .= _step(σ, t_[1], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)
    Utils.transform!(σ, V, tmpV)
    for i ∈ 2:r
        σ .= _step(σ, t_[i], τ, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)
        Utils.transform!(σ, V, tmpV)
    end
    σ .= _step(σ, t_[end], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)

    Utils.transform!(σ, UD', tmpV)      # ROTATE INTO DEVICE BASIS
    σ .*= exp.( (im*T) * ΛD)            # ROTATE PHASES FOR FINAL exp(𝒊 HD T)
    Utils.transform!(σ, UD, tmpV)       # ROTATE OUT OF DEVICE BASIS
    #= TODO: Eventually H may be in device or qubit basis, interaction or lab frame
                but for now I'm assuming it's in qubit basis, interaction frame. =#

    Utils.transform!(σ, observable, tmpV)       # CALCULATE Ô|ψ⟩

    # REVERSE TIME-EVOLUTION UNDER TRAPEZOIDAL RULE
    Utils.transform!(σ, UD', tmpV)      # ROTATE INTO DEVICE BASIS
    σ .*= exp.( (-im*T) * ΛD)           # ROTATE PHASES FOR FIRST exp(-𝒊 HD T)
    Utils.transform!(σ, UD, tmpV)       # ROTATE OUT OF DEVICE BASIS

    σ .= _step(σ, t_[end], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_, true)
    Utils.transform!(σ, V', tmpV)
    for i ∈ reverse(2:r)
        σ .= _step(σ, t_[i], τ, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_, true)
        Utils.transform!(σ, V', tmpV)
    end
    σ .= _step(σ, t_[1], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_, true)

    ######################################################################################
    #                               GRADIENT CALCULATION
    _fillgradients!(∂Ω, 1, σ,ψ, t_, τ/2, pulses, qubitapplymode, n, a, UD, ΛD, tmpV,tmpM_,tmpK_)

    σ .= _step(σ, t_[1], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)
    ψ .= _step(ψ, t_[1], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)

    Utils.transform!(σ, V, tmpV)
    Utils.transform!(ψ, V, tmpV)
    for i ∈ 2:r
        _fillgradients!(∂Ω, i, σ,ψ, t_, τ, pulses, qubitapplymode, n, a, UD, ΛD, tmpV,tmpM_,tmpK_)

        σ .= _step(σ, t_[i], τ, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)
        ψ .= _step(ψ, t_[i], τ, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)

        Utils.transform!(σ, V, tmpV)
        Utils.transform!(ψ, V, tmpV)
    end

    _fillgradients!(∂Ω, r+1, σ,ψ, t_, τ/2, pulses, qubitapplymode, n, a, UD, ΛD, tmpV,tmpM_,tmpK_)
    # NO REASON TO PERFORM THE LAST TIME-EVOLUTION STEP...

    # TEMP: Do it anyways..?
    σ .= _step(σ, t_[end], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)
    ψ .= _step(ψ, t_[end], τ/2, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)

    ######################################################################################

    return ∂Ω
end




""" Auxiliary function to evolve a single step in time. """
function _step(ψ, t, τ, pulses,
        qubitapplymode, n, a, tmpV, tmpM_, tmpK_, adjoint=false)
    ######################################################################################
    #                                 SINGLE TIME STEP

    # PREPARE QUBIT DRIVES
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)

        # CONSTRUCT AND EXPONENTIATE MATRIX
        tmpM_[q] .= z .* a  # ADD za TERM
                            # THE z' a' TERM IS ACCOUNTED FOR BY THE `Hermitian` VIEW

        tmpM_[q] .= exp(( ((-1)^adjoint)* -im*τ) .* Hermitian(tmpM_[q]))
            # THIS LAST STEP SHOULD BE THE ONLY ONE REQUIRING ANY ALLOCATIONS
    end

    # APPLY QUBIT DRIVES
    if qubitapplymode isa Kronec
        # KRONECKER MODE: CONSTRUCT FULL-BODY OPERATOR
        O = Utils.kron_concat(tmpM_, tmpK_)
        return mul!(tmpV, O, ψ)
    elseif qubitapplymode isa Tensor
        # TENSOR MODE: RESHAPE AND CONTRACT
        ψ_ = reshape(ψ, tmpK_[1])   # *NOT* A COPY; MUTATIONS APPLY TO BOTH
        ψ_ .= ncon(
            [tmpM_..., ψ_],                         # LIST OF TENSORS
            tmpK_[2],    # LIST OF INDICES ON EACH TENSOR
            tmpK_[4], :cache,                       # ENABLE CACHING
            output=tmpK_[3],                        # FINAL PERMUTATION
        )
        # ψ HAS ALREADY BEEN UPDATED, IN MUTATIONS OF ψ_
        return ψ
    else
        error("Invalid `QubitApplyMode` object. (How did you manage that???)")
    end
end


""" Auxiliary function to calculate gradients for a single step in time. """
function _fillgradients!(∂Ω, i, σ, ψ, t_, τ, pulses,
        qubitapplymode, n, a, UD, ΛD, tmpV, tmpM_, tmpK_)
    # INITIALIZE ALL QUBIT OPERATORS TO I
    for q ∈ 1:n
        tmpM_[q] .= one(a)
    end
    #= TODO: Right now we're assigning all qubit operators for each gradient value,
            so that we can re-use all the code we already have for time evolution.
        But actually here we can do at least the tensor contraction much much faster,
            since we only need to do a single contraction for each gradient value.
            (ie. only need to permute and un-permute dimensions once).
        In fact, I think we can do the dimension permutation in a cyclic fashion
            to calculate *each* gradient value successively,
            to get a total runtime as low as a single time-step.
        But that can wait 'til we've got our own manual tensor contraction.

        In thinking about the kron method, I realized the `on` method should have
            a pre-allocation option available similar to kron_concat,
            but it wouldn't be any faster, I think, than just doing kron_concat
            with a bunch of identity matrices, so we needn't mess with it.
    =#

    # PREPARE QUBIT DRIVES
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        ν = Pulses.frequency(pulses[q], t_[i])
        z = exp(im*ν*t_[i])

        # CONSTRUCT MATRIX
        # a = zeros(eltype(a), size(a))
        # a[1,2] = 1
        tmpM_[q] .= z .* a                  # ADD za TERM
        tmpM_[q] .= Hermitian(tmpM_[q])     # Hermitian ADDS a'

        # APPLY AS OPERATOR
        if qubitapplymode isa Kronec
            # KRONECKER MODE: CONSTRUCT FULL-BODY OPERATOR
            O = Utils.kron_concat(tmpM_, tmpK_)
            mul!(tmpV, O, ψ)
        elseif qubitapplymode isa Tensor
            tmpV .= ψ       # BE CAREFUL NOT TO MUTATE ψ BEFORE WE MEAN TO!
            # TENSOR MODE: RESHAPE AND CONTRACT
            ψ_ = reshape(tmpV, tmpK_[1])    # *NOT* A COPY; MUTATIONS APPLY TO BOTH
            ψ_ .= ncon(
                [tmpM_[q], ψ_],                         # LIST OF TENSORS
                tmpK_[2],    # LIST OF INDICES ON EACH TENSOR
                tmpK_[4], :cache,                       # ENABLE CACHING
                output=tmpK_[3],                        # FINAL PERMUTATION
            )
            # tmpV HAS ALREADY BEEN UPDATED, IN MUTATIONS OF ψ_
        else
            error("Invalid `QubitApplyMode` object. (How did you manage that???)")
        end

        # CALCULATE GRADIENT
        σAψ = -im*τ * (σ' * tmpV)       # THE BRAKET ⟨σ|A|ψ⟩

        # # TEMP: GLOBAL PHASE IS NOT GLOBAL THROUGH PROJECTION
        # σAψ *= exp(-im*angle(tmpV[1]))
        # σAψ *= exp(+im*angle(   σ[1]))


        ∂Ω[i,q] = σAψ + σAψ'            # THE GRADIENT ⟨σ|A|ψ⟩ + ⟨ψ|A|σ⟩

        # REPLACE QUBIT OPERATOR WITH IDENTITY
        tmpM_[q] .= one(a)
    end
end




end # END MODULE
