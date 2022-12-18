#= Provide an interface to PySCF calculations.

Input: geometry, or molecule+parameters
Output: h[p,q], h[p,q,r,s], E_HF, E_FCI

The interface with PySCF is intentionally Klingon-ish and restrictive,
    to motivate a pure Julianic solution sooner than later!
For example, I'm assuming sto-3g basis always, and no use for MO-AO transformation.
Personally I think atomic orbitals should be the starting point of quantum algorithms,
    since they are sort of the last "physically intuitive" object.
But I'm just interested in catching the code up to ctrlq right now, so this will do.

=#

module Molecules

import LinearAlgebra: I, Hermitian, norm

import TensorOperations: @tensor
import PyCall: pyimport

import ..Utils

GeometryType = Vector{Pair{String, Tuple{Float64,Float64,Float64}}}

function H2_geometry(r::Real)::GeometryType
    return [
        "H" => (zero(typeof(r)), zero(typeof(r)), zero(typeof(r))),
        "H" => (zero(typeof(r)), zero(typeof(r)),             r  ),
    ]
end

struct Molecule
    geometry::GeometryType      # ATOMIC LABELS PAIRED WITH COORDINATES
    N::Int                      # SIZE OF BASIS SET, OR NUMBER OF *SPATIAL* ORBITALS
    n::Int                      # TOTAL NUMBER OF *SPIN* ORBITALS
    Ne::Int                     # TOTAL NUMBER OF ELECTRONS
    Nα::Int                     # NUMBER OF SPIN-UP   ELECTRONS
    Nβ::Int                     # NUMBER OF SPIN-DOWN ELECTRONS
    E_HF::Float64               # HARTREE-FOCK ENERGY
    E_FCI::Float64              # FULL CONFIGURATION INTERACTION ENERGY
    h0::Float64                 # NUCLEAR REPULSION ENERGY
    h1::Array{Float64,2}        # SINGLE-BODY TENSOR
    h2::Array{Float64,4}        # TWO-BODY TENSOR
    # H = h0 + ∑ h1[p,q] a†[p] a[q] + ∑ h2[p,q,r,s] a†[p] a†[q] a[r] a[s]
end

function Molecule(geometry::GeometryType; skip_FCI=false)
    # CREATE MOLECULE
    ######################################################################################
    mol = pyimport("pyscf.gto").Mole()
    mol.atom = geometry     # TODO: Verify pair converts to tuple, or add comprehension.
    mol.basis = "sto-3g"
    mol.build()
        # NOTE: No options for ions or multiplicities. Include in the Julianic solution!
    N = mol.nao_nr().tolist() # nao_nr returns scalar numpy array sometimes??
    n = 2N
    Nα, Nβ = mol.nelec
    Ne = Nα + Nβ

    # RUN CALCULATIONS
    ######################################################################################
    mf = pyimport("pyscf.scf").RHF(mol).run()
    E_HF = mf.e_tot

    if !skip_FCI
        fci = pyimport("pyscf.fci").FCI(mf).run()
        E_FCI = fci.e_tot
    end

    # PROCESS INTEGRALS INTO HAMILTONIAN TENSORS
    ######################################################################################

    # CALCULATE ATOMIC-ORBITAL INTEGRALS
    T = mol.intor("int1e_kin_sph")      # KINETIC ENERGY
    U = mol.intor("int1e_nuc_sph")      # NUCLEAR REPULSION
    h = T + U                           # SINGLE-BODY INTEGRALS
    v = mol.intor("int2e_sph")          # TWO-BODY INTEGRALS

    # TRANSFORM INTO MOLECULAR-ORBITAL BASIS
    C = mf.mo_coeff                    # AO -> MO TRANSFORMATION MATRIX
    @tensor H[i,j] := C'[i,p] * h[p,q] * C[q,j]
    @tensor V[i,j,k,l] := v[p,q,r,s] * C[p,i] * C[q,j] * C[r,k] * C[s,l]

    # CONSTRUCT TENSORS
    h0 = mol.energy_nuc()               # NUCLEAR REPULSION

    h1 = zeros(Float64, n, n)           # ONE-BODY TENSOR
    for i in 1:N; for j in 1:N
        h1[2i-1,2j-1]     = H[i,j]      # SPIN-UP   TERMS
        h1[2i  ,2j  ] = H[i,j]          # SPIN-DOWN TERMS
    end; end

    h2 = zeros(Float64, n, n, n, n)     # TWO-BODY TENSOR
    for i in 1:N; for j in 1:N; for k in 1:N; for l in 1:N
        h2[2i-1, 2j-1, 2l-1, 2k-1] = V[i,k,j,l] / 2
        h2[2i-1, 2j  , 2l  , 2k-1] = V[i,k,j,l] / 2
        h2[2i  , 2j-1, 2l-1, 2k  ] = V[i,k,j,l] / 2
        h2[2i  , 2j  , 2l  , 2k  ] = V[i,k,j,l] / 2
    end; end; end; end
        #= Following code in ADAPT repo but I still think the index swap is a bit odd.

                Swapping k and j on right switches
                    from chemist's a†a a†a to physicist's normal-ordering.
                Why swap k and l on the left, though?
                    Apparently it's some convention to make adjoints easier to think about,
                        but I don't get it yet...
        =#

        #= NOTE: This is a highly wasteful tensor entirely ignorant of fermionic symmetries.

            Eg. [0,1,1,0] represents the same as -[1,0,1,0], which is also [1,0,0,1]
                Even more poignantly, [0,0,0,0] is the same as... 0. ^_^
                But, this'll do for a brute-force solution.

                DO NOT REPLICATE IN THE PURE JULIANIC VERSION TO COME!!! ^_^
        =#

    return Molecule(geometry, N, n, Ne, Nα, Nβ, E_HF, E_FCI, h0, h1, h2)
end

function fermi_a(q,n)
    Iq = [1.0  0.0; 0.0  1.0]       # IDENTITY
    Zq = [1.0  0.0; 0.0 -1.0]       # PARITY
    aq = [0.0  1.0; 0.0  0.0]       # ANNIHILATOR

    ops = [ (p < q) ? Iq : (p > q) ? Zq : aq for p in 1:n]
    return Utils.kron_concat(ops)
end

function molecular_hamiltonian(mol::Molecule; m=nothing)
    #= m is number of levels to extend 2-level fermionic operator to.
        Only use this for calculating gradient,
        where we need to apply H to transmon state-vector in the middle of time evolution.
    =#
    # NOTE: This code is as brute-force as it gets. DO NOT REPLICATE!!!

    n = mol.n           # NUMBER OF QUBITS == 2 * NUMBER OF SPATIAL ORBITALS
    N = 2^n             # SIZE OF HILBERT SPACE (NOTE: *NOT* THE SAME AS mol.N)

    H = mol.h0 * Matrix(I,N,N)                              # h0

    for i in 1:n; for j in 1:n                              # h1
        if abs(mol.h1[i,j]) < eps(Float64); continue; end
        H += mol.h1[i,j] * fermi_a(i,n)' * fermi_a(j,n)
    end; end

    for i in 1:n; for j in 1:n; for k in 1:n; for l in 1:n  # h2
        if abs(mol.h2[i,j,k,l]) < eps(Float64); continue; end
        H += mol.h2[i,j,k,l] * fermi_a(i,n)' * fermi_a(j,n)' * fermi_a(k,n) * fermi_a(l,n)
    end; end; end; end

    if m !== nothing
        H = Utils.extendoperator(H, n, m)
    end

    return Hermitian(H)
end

function measure_energy(H::Hermitian, ψ::Vector{ComplexF64})
    N = size(H,1)
    # IF NEEDED, PROJECT Ψ ONTO QUBIT SPACE
    if length(ψ) > N
        n = Int(ceil(log2(N)))              # NUMBER OF QUBITS
        ψ = Utils.qubitspace(ψ, n)          # PROJECTED STATEVECTOR
        ψ ./= norm(ψ)                       # RE-NORMALIZED
    end
    return real(Utils.expectation(H,ψ))     # ENERGY CALCULATION
end






end # END MODULE
