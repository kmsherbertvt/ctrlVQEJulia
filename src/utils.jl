#= Provide matrix representations for common Hamiltonians. =#
module Utils
import LinearAlgebra: kron, I, eigen, Eigen


"""
    a_matrix(nstates::Integer=2)

Matrix representation of the bosonic annihilation operator on a single qubit.

Note that the matrix representation must be truncated to `nstates` levels.

"""
function a_matrix(nstates::Integer=2)
    a = zeros((nstates,nstates))
    for i ∈ 1:nstates-1
        a[i,i+1] = √i               # BOSONIC ANNIHILATION OPERATOR
    end
    return a
end

"""
    on(op::AbstractMatrix{<:Number}, q::Integer, n::Integer)

Expand the single-qubit matrix operator `op` to act on qubit `q` of `n` qubits.

In other words, apply Kronecker-products such that identity `I` acts on the other qubits.

"""
function on(op::AbstractMatrix{<:Number}, q::Integer, n::Integer)
    A = ones(1,1)                   # A 1x1 IDENTITY MATRIX
    I = one(op)                     # AN IDENTITY MATRIX MATCHING THE DIMENSIONS OF `op`
    for i ∈ 1:n
        A = kron(A, i == q ? op : I)    # `op` ACTS ON QUBIT `q`, `I` ACTS ON ALL OTHERS
    end
    return A
end

"""
    kron_concat(ops::AbstractVector{<:AbstractMatrix{<:Number}})

Concatenate a sequence of operators with the Kronecker product.

"""
function kron_concat(ops::AbstractVector{<:AbstractMatrix{<:Number}})
    O = Matrix(I,1,1)
    for q ∈ eachindex(ops)
        O = kron(O, ops[q])
    end
    return O
    #= TODO: Pre-allocate O, ...and each intermediate. :?
        There *must* be an easy way to kron a sequence of matrices together...
        Maybe reshaping O into a tensor..? Ironic...
    =#
end

"""
    kron_concat(ops::AbstractMatrix{<:Number}, n::Integer)

Concatenate a repeated string of an operator with the Kronecker product.

"""
function kron_concat(op::AbstractMatrix{<:Number}, n::Integer)
    O = Matrix(I,1,1)
    for q ∈ 1:n
        O = kron(O, op)
    end
    return O
    #= TODO: Pre-allocate O, ...and each intermediate. :?
        There *must* be an easy way to kron a sequence of matrices together...
        Maybe reshaping O into a tensor..? Ironic...
    =#
end



"""
    algebra(
        n::Integer,
        m::Integer=2;
        basis::Union{AbstractMatrix{<:Number},Nothing}=nothing,
    )

Construct a vector of annihilation operators acting on each of n m-level systems.

Optionally, rotate these operators into the provided basis.

These matrices, in conjunction with their adjoints,
    form a complete algebra for the space of n m-level systems.

"""
function algebra(
    n::Integer,
    m::Integer=2;
    basis::Union{AbstractMatrix{<:Number},Nothing}=nothing,
)
    a_ = a_matrix(m)                        # SINGLE-QUBIT ANNIHILATION OPERATOR
    a = [on(a_, q, n) for q in 1:n]         # EACH OPERATOR, ACTING ON FULL HILBERT SPACE
    if !(basis === nothing)
        for q ∈ 1:n
            a[q] .= basis' * a[q] * basis   # CONJUGATE WITH BASIS
        end
    end
    return a
end


"""
    dressedbasis(H::AbstractMatrix{<:Number})

Calculate the eigenvalues and eigenvectors for a Hermitian matrix `H`.

This method returns an `Eigen` object, the same type as LinearAlgebra.eigen(H).

Dressing does three things:
    1. Re-order eigenpairs such that the ith eigenvector is the one with
        the largest magnitude in the ith element.
       (This is distinct from the LinearAlgebra.eigen order,
            which uses ascending order of eigenvalues.)
    2. Multiply each eigenvector by a global phase such that the ith element
        is real and positive.
    3. Zero out any elements with magnitude less than machine precision.

"""
function dressedbasis(H::AbstractMatrix{<:Number})
    N = size(H)[1]
    Λ, U = eigen(H)
    #= TODO: I don't know if this is sophisticated enough to use Hermitian algorithm when H is Hermitian. We should do some testing to find out... =#

    # IMPOSE PERMUTATION
    σ = Vector{Int}(undef,N)
    for i in 1:N
        perm = sortperm(abs.(U[i,:]), rev=true) # STABLE SORT BY "ith" COMPONENT
        perm_ix = 1                             # CAREFULLY HANDLE TIES
        while perm[perm_ix] ∈ σ[1:i-1]
            perm_ix += 1
        end
        σ[i] = perm[perm_ix]
    end
    Λ .= Λ[  σ]
    U .= U[:,σ]

    #= TODO: This strategy for selecting the permutation might not be the most efficient.
        It's not *bad*: should be no worse than O(N² log N).
        But a less careful handling of ties would be only O(N²). I'm not sure. :?
    =#

    # IMPOSE PHASE
    for i in 1:N
        U[:,i] .*= U[i,i] < 0 ? -1 : 1           # ASSUMES REAL TYPE
        # U[:,i] .*= exp(-im*angle(U[i,i]))        # ASSUMES COMPLEX TYPE
    end

    # IMPOSE ZEROS
    Λ[abs.(Λ) .< eps(eltype(Λ))] .= zero(eltype(Λ))
    U[abs.(U) .< eps(eltype(U))] .= zero(eltype(U))
    return Eigen(Λ,U)

    #= TODO:
        Each of these impositions should properly be their own function.
        The phase and zeros impositions should have two methods each, Real vs Complex types.

            real phase should just multiply U[:,i] by -1 if U[i,i] < 0.
            complex phase does as written
                except use .*= operator once you're guaranteed types needn't change.

            complex zeros should reinterpret as (larger) Float array
                and zero out all small components.
            I don't know exactly how to do the reinterpret with the correct typing,
                short of writing a different method for each possible complex type.
            I absolutely do *not* want to assume Float64.
                It's *probably* true but I don't trust it!

            phase need only ever take matrices, but zeros takes anything
            Unfortunately, knowing the NUMBER type probably means
                having to dispatch separately on the array type.
            I think we can do that just fine with a Union rather than separate methods,
                but you can see how this whole thing gets WAY more complicated
                than it should be when you want to zero out individual components.
            abs turns out to be WAY easier...
    =#
end


end
