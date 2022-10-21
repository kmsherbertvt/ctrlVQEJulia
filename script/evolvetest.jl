#= We want to confirm 'trotter' method duplicated from ctrlq code is functional.



In particular, evolving 2 qubits with 3 levels each starting in the state |11⟩
    [0, 0, 0, 0, 1, 0, 0, 0, 0]
for 10ns under a pulse pattern
    29.0 GHz [ 0.5 |2.3ns| 0.4 |5.6ns| 0.3 |6.1ns| 0.2 ]
    31.0 GHz [ 0.3 |1.0ns| 0.1 |3.4ns| 0.5 |5.6ns| 0.6 ]
on a device with
    (ω = 2π 4.808049015463495, δ = 2π 0.3101773613134229)
                =| g = 2π 0.018312874435769682 |=
    (ω = 2π 4.833254817254613, δ = 2π 0.2916170385725456)
with 2000 Trotter steps should produce a statevector a bit like:
    [-0.15351954+0.23975274j  0.10354723-0.1407912j   0.05421186+0.09792655j
      0.63627648-0.11943777j -0.5198033 -0.34088801j -0.05686797-0.03488353j
      0.1904338 -0.02687666j  0.15815612+0.00389777j  0.05597458-0.02764928j]
up to a global phase.

=#

include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")

module TrotterValidation
import ..Devices
import ..Pulses
import ..Evolutions

##########################################################################################
#                               EXPECTED FINAL RESULT
ψ_ctrlq = Dict(   # MAPS DURATION OF PULSE TO EXPECTED STATEVECTOR
     0.0 => [0, 0, 0, 0, 1.0, 0, 0, 0, 0],
     0.1 => [
        -1.67480756e-03+6.42397602e-05im, -2.01210100e-03-6.34942126e-02im,
        -1.66531657e-03+3.83130079e-04im,  1.99008074e-03-6.11714361e-03im,
         9.94245767e-01-1.68308318e-04im, -5.12187158e-03-7.58555789e-02im,
        -1.95361042e-03+1.36588065e-05im,  6.54234910e-05-3.99805000e-02im,
        -3.13291742e-03+5.25145054e-04im,
     ],
     0.5 => [
        -0.03814677+0.00732511im, -0.04807000-0.28225712im, -0.01820459+0.03649146im,
         0.04278763-0.03019730im,  0.86703808-0.01916422im, -0.10843451-0.31243879im,
        -0.04627545+0.00167817im, -0.01151319-0.20789797im, -0.04855472+0.05245380im,
    ],
     1.0 => [
        -0.11473402+0.04389236im, -0.16671459-0.39486411im,  0.06837837+0.11447633im,
         0.10187140-0.05832751im,  0.57461647-0.11692289im, -0.25736235-0.35018002im,
        -0.15832829+0.01367583im, -0.13897147-0.38795092im,  0.00763658+0.20442401im,
    ],
     5.0 => [
        -0.01497412-0.35080042im,  0.06454326-0.53487738im,  0.09211357-0.08776135im,
        -0.00901190-0.09606471im, -0.20740631-0.45950085im, -0.13465213+0.38943825im,
         0.24062247+0.01153718im,  0.12970638-0.24059826im,  0.00231078-0.06525294im,
    ],
    10.0 => [
        -0.15351954+0.23975274im,  0.10354723-0.14079120im,  0.05421186+0.09792655im,
         0.63627648-0.11943777im, -0.51980330-0.34088801im, -0.05686797-0.03488353im,
         0.19043380-0.02687666im,  0.15815612+0.00389777im,  0.05597458-0.02764928im,
    ],
    20.0 => [
        -0.07793078-0.17825898im, -0.12078268-0.01304937im, -0.00625089-0.04476663im,
         0.49426001-0.30091994im, -0.47424059-0.50686593im,  0.08929867-0.25512956im,
         0.05882941-0.18257537im,  0.13132036+0.03093926im,  0.02495000+0.00095248im,
    ]
)


fidelity(ψ,φ) = 1 - abs2(ψ'*φ)

# CONSTRUCT INITIAL STATE
ψI::Vector{Number} = [0, 0, 0, 0, 1, 0, 0, 0, 0]


# DESIGN PULSE PATTERNS
T = 10.0
pulses = [
    Pulses.BasicSquarePulse(T, 29.0, [0.5, 0.4, 0.3, 0.2], [2.3, 5.6, 6.1]),
    Pulses.BasicSquarePulse(T, 31.0, [0.3, 0.1, 0.5, 0.6], [1.0, 3.4, 5.6]),
]


# # JUST FOR CONTEXT
# display(ψ_ctrlq[T])
# println("Distance should be much much lower than $(fidelity(ψI, ψ_ctrlq[T]))")



# DESIGN DEVICE
device = Devices.Transmon(
    2π .* [4.8080490154634950, 4.8332548172546130],
    2π .* [0.3101773613134229, 0.2916170385725456],
    Dict(
        Devices.QubitCouple(1,2) => 2π*0.018312874435769682,
    ),
)








##########################################################################################
#                           ROTATION INTO/OUT OF ctrlq'S BASIS

# CONSTRUCT AND DIAGONALIZE THE DEVICE HAMILTONIAN
import LinearAlgebra: eigen
nstates = 3; n=2; N = nstates^n
HD = Devices.static_hamiltonian(device, nstates)
ΛD, UD = eigen(HD)

# REORDER BASES BY OVERLAP WITH COMPUTATIONAL BASIS
#   ALSO SELECT SIGN SO THAT DIAGAONAL IS NON-NEGATIVE
# TEMP: To match ctrlq, let's re-order basis states by overlap with computational basis. (..why would you do this? Is there some computational advantage to it, or did Oinam just want to impose some sense onto the eigenbasis? Did he know numpy already sorts them by eigenvalue? I'd guess so, which would imply he didn't like that ordering... must find out why.)
order = [i for (i,j) ∈ Tuple.(argmax(abs.(UD), dims=1))]
new_UD = Matrix{Number}(undef,N,N)
for i ∈ 1:N
    new_UD[:,order[i]] = UD[:,i] * sign(UD[order[i],i])
end
UD = new_UD

# RUNNING THE ctrlq CODE, THE INITIAL STATE WAS IMPLICITLY IN ctrlq BASIS
#   MY CODE EXPECTS COMPUTATIONAL BASIS, SO WE MUST ROTATE *OUT* OF ctrlq BASIS
ψI = UD * ψI

##########################################################################################
#                           PERFORM THE ACTUAL EVOLUTION
ψx = Evolutions.evolve(ψI, pulses, device, Evolutions.Trotter; numsteps=2000)

# MY CODE OUTPUTS IN THE COMPUTATIONAL BASIS, SO WE MUST ROTATE *INTO* ctrlq BASIS
ψx = UD' * ψx




##########################################################################################
#                           CHECK FIDELITY WITH EXPECTED RESULT

# display(ψx)
println("Distance from ψI: $(fidelity(ψx, ψI))")
println("Distance from ctrlq: $(fidelity(ψx, ψ_ctrlq[T]))")


end
