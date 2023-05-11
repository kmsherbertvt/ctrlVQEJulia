# Energy Functions
##### 5/11/2023

In another note (`timeevolution.md`), we have discussed algorithms to simulate time evolution of a quantum-computational state under some controllable drive Hamiltonian.
The next step in ctrl-VQE is to measure the expectation value $E$ of an observable, typically a molecular Hamiltonian.
The ultimate goal is to refine our control signals so that we prepare a state whose expectation value $E$ is as small as possible.
But before we get to that, we need to address some technical issues with measuring $E$ correctly, related to the distinction between logical vs physical states.

## Logical States and Physical States
Let's say we have a Hermitian operator $\hat O$ which acts on a _logical state_ $|Ψ⟩$.
Our goal is to measure the expectation value:

$$ E ≡ ⟨Ψ|\hat O|Ψ⟩ $$

For quantum chemistry applications,
    the operator $\hat O$ is often the second-quantized fermionic Hamiltonian
    associated with a particular molecule in a particular geometry,
    mapped onto qubit operators in a way that preserves fermionic anti-symmetry (eg. with the Jordan-Wigner transformation).
Thus, the logical state lives in the Hilbert space spanned by some number $n$ logical qubits.

On the other hand, we have prepared (after time evolution) a _physical state_ $|ψ⟩$.
While a _logical_ qubit has only two potential states (labeled $|0⟩$, $|1⟩$),
    the "physical qubit" used to _implement_ a qubit will usually have a multitude of additional potential states.
The physical state $|ψ⟩$ lives in the much larger Hilbert space spanned by each physical system.
There are as many ways of encoding logical states into physical states as there are quantum computers,
    but in this note, we'll label the states for each physical qubit as $|0⟩$, $|1⟩$, $|2⟩$, $|3⟩$, ...,
    and we'll define our logical space as the span of the first two ($|0⟩$ and $|1⟩$) from each qubit.
In order to measure $E$, we need to define some way of mapping any given physical state $|ψ⟩$ onto a logical state $|Ψ⟩$,
    which will involve some form of projecting out unwanted states.
There are multiple strategies for doing so.

The other discrepancy between logical and physical states relates to the _reference frame_.
Consider this: if we have applied some control fields to prepare some logical state,
    and then we _stop_ applying control fields,
    the logical state should not change from then on.
On the other hand, the physical qubits continue to evolve under the static device Hamiltonian $H_0$ featured in the `timeevolution.md` note.
In this sense, the physical state _does_ change, and we must account for this in our mapping from $|ψ⟩$ to $|Ψ⟩$.
There are, once again, multiple strategies for doing so.

### Projection Strategies

Neglecting for the time being the latter problem of reference frames, let us survey strategies for projecting out unwanted states.

#### Projected

While all these strategies are projection strategies,
    I consider this first approach to be the "default" - the most intuitive and also the most desirable when it can be realized experimentally.
Therefore, I've christened it the "Projected" strategy.

Let the single-qubit projection operator $π≡|0⟩⟨0|+|1⟩⟨1|$, and the $n$-qubit projection operator $Π≡π^{⊗n}$.
Now let the "projected" logical state $|Ψ_P⟩$ be the physical state $|ψ⟩$ restricted to the logical Hilbert space:

$$ |Ψ_Π⟩ ≡ Π|ψ⟩ $$

If the physical state $|ψ⟩$ has any support on states outside the logical Hilbert space,
    this projection results in an _unnormalized_ logical state $|Ψ_Π⟩$,
    meaning the expectation value $⟨Ψ_Π|\hat O|Ψ_Π⟩$ is not a _physically meaningful_ energy.

_However_, that is actually okay for typical VQE experiments.
The unnormalization of $|Ψ⟩$ always serves to bring $E$ closer to zero than it otherwise would be.
Thus, as long as the ground state energy is negative (as it always is for bound chemical systems),
    the energy is only ever minimized by normalized logical states.
In other words, the optimization procedure in VQE will naturally penalize "leakage" to states outside the logical Hilbert space.
This is quite desirable, and as such, this strategy is typically the one we'd _like_ to use when possible.

#### Normalized

Consider the same projection operator from above.
Now let the normalized logical state be $|Ψ_N⟩$ be:

$$ |Ψ_N⟩ ≡ Π|ψ⟩ / \sqrt{⟨ψ|Π|ψ⟩}$$

The denominator ensures that $|Ψ_N⟩$ is normalized.
This is easily shown: calculate $⟨Ψ_N|Ψ_N⟩$, remembering that projection operators are Hermitian and idempotent (ie. $Π^† Π=Π$).

This projection protocol ensures every energy $E$ is always a physical one,
    but in doing so, it "hides" any leakage from the optimizer.
Thus, in principle, optimizing an energy derived from this projection may well lead to very ill-behaved physical states.

#### Biased

The "Projected" and "Normalized" strategies are well-suited to numerical simulations,
    but they may or may not be experimentally realizable, depending on the measurement protocol.

Let's model the measurement protocol as a function $M$ which takes some physical state
    and outputs a vector of "results" for each qubit.
For example, for a single qubit $n=1$, we may find $M(|ψ⟩) = 0$, or $M(|ψ⟩) = 1$.
If $|ψ⟩$ is a superposition of multiple basis vectors,
    $M(|ψ⟩)$ is probabilistic, and it must be performed many times to gather good statistics
    (re-preparing $|ψ⟩$ each time, of course. Quantum experiments are hard!).

Now, what happens if $|ψ⟩$ has support on states outside the logical Hilbert space?
There are at least three possibilities, depending on the sophistication of one's experimental apparatus:
1. The protocol can discriminate every potential physical state: $M(|ψ⟩)$ could be 0, 1, 2, 3, etc.
2. The protocol discriminates logical states from non-logical states: $M(|ψ⟩)$ could be 0, 1, or FAILURE.
3. The protocol artificially interprets non-logical states as logical states: $M(|ψ⟩)$ could be 0 or 1 only.

In either of the first two possibilities, the "Projected" and "Normalized" strategies are both viable.
Unfortunately, the third possibility seems to be the more common one.
A proper treatment would require detailed knowledge of the measurement protocol,
    and the relative probabilities that each physical basis state is measured as each binary vector,
    but a somewhat justifiable heuristic is to assume that any non-logical states are measured as 1.
Thus, any leakage artificially biases the expectation value toward more "occupied" orbitals.

Using this heuristic, we can define the biased single-qubit projection operator $π_B≡|0⟩⟨0| + \sum_{i=1}^{∞}|1⟩⟨i|$,
    and the biased $n$-qubit projection operator as $Π_B≡π_B^{⊗n}$.
Then the biased logical state $|Ψ_B⟩$ is:

$$ |Ψ_B⟩ = Π_B|ψ⟩$$

Note that $|Ψ_B⟩$ is automatically normalized, so leakage is apparently invisible to the optimizer.
However, this particular choice of bias always serves to increase particle number,
    so in typical chemistry applications we can still penalize it by including additional terms like $⟨Ψ_B|(\hat N-N_e)|Ψ_B⟩$ to our cost function,
    where $\hat N$ is the particle number operator and $N_e$ is the number of electrons we expect to see in our molecule.
Indeed, we'd probably like to include such terms _anyway_, to guarantee we find states satisfying expected symmetries.



### Frame Rotation Strategies

Let us now survey strategies to account for the reference frame.

We _want_ the logical state to be "at rest" with respect to the device's static Hamiltonian.
We will call this the "rotating" frame, indicated by a subscript $_R$,
    while the one "at rest" with nothing in particular we will call the "lab frame", indicated by a lack of subscript.
Rotating from one frame to another at time $T$ is done via the interaction picture formalism, applying rotations of $e^{±iTH_0}$.
For example the physical state represented in each frame may be written in terms of the other as follows:

$$ |ψ⟩ → e^{-iT H_0} |ψ_R⟩ $$

If the logical state is to be at rest with respect to the rotating frame,
    we should apply our projection strategy
    to the physical state in the rotating frame, ie. $|ψ_R⟩$ rather than $|ψ⟩$:

$$ |Ψ⟩ = Π |ψ_R⟩ $$

(I'm writing just $Π$ here for the "Projected" strategy,
    but you may include a normalization factor or substitute $Π→Π_B$,
    according to your preferred strategy.)

This is perfectly fine for simulation -
    indeed, the `timeevolution.md` notes use the interaction picture throughout,
    and its conclusions give $|ψ_R⟩$ rather than $|ψ⟩$!
But is it tractable for experiment?

Honestly, we're not quite sure.
But it seems like measurement happens _to_ the device, not _within_ the device,
    and the _thing which is measured_ is $|ψ⟩$, not $|ψ_R⟩$!
Thus, we should write our logical state out explicitly as:

$$ |Ψ⟩ = Π e^{iT H_0} |ψ⟩ $$

We need to explicitly include a frame rotation _prior_ to applying one of the projection strategies discussed above.
So how do we handle this?

#### Reverse Evolution

Perhaps the most obvious option is to simply prepare the state $e^{iT H_0} |ψ⟩$ instead of $|ψ⟩$.
_In silico_, this is trivial: we simply perform the projection step on $|ψ_R⟩$.
Experimentaly, this strategy requires finding a control signal which effects the unitary evolution $e^{iT H_0}$.

Note that the control signal effecting the _inverse_ evolution $e^{-iT H_0}$ is trivial:
    you simply do _nothing_ for a time $T$, allowing the device to evolve under its static Hamiltonian.
So the goal is simply to "reverse time" on your device.
I feel like there may be a clever way to just _do_ that -
    like, reverse the polarity or something. I dunno.
But if not, you need to find the pulses which explicitly implement the target unitary $e^{iT H_0}$.
This is a traditional problem in quantum control,
    and might be left to the hardware people,
    since the same pulses should work for any choice of variational parameters,
    as long as the pulse duration remains fixed.
But eliminating that quantum control problem is one of the primary motivations for ctrl-VQE,
    and it is rather annoying to have it crop up at the very end...

#### Virtual Phase Rotations

Let us consider a special case, using the model transmon Hamiltonian from `timeevolution.md`:

$$ H_0 = \sum_q ( ω_q \hat a_q^† \hat a_q - \frac{δ_q}{2} \hat a_q^† \hat a_q^† \hat a_q \hat a_q ) + \sum_{⟨pq⟩} g_{pq} (\hat a_p^† \hat a_q + \hat a_q^† \hat a_p ) $$

Recall that $\hat a_q$ is the bosonic annihilation operator acting on transmon $q$,
    $ω_q$ and $δ_q$ are the resonance frequency and anharmonicity of transmon $q$,
    and $g_{pq}$ is the coupling strength between transmons $p$ and $q$.
If we neglect the coupling terms (ie. set $g_{pq}≃0$),
    the unitary $e^{iT H_0}$ reduces to applying a diagonal phase rotation operator $u_q$ on each qubit:

$$ u_q = e^{i ω_q T \hat a_q^† \hat a_q} e^{-i δ_q T \hat a_q^† \hat a_q^† \hat a_q \hat a_q} $$

When restricted to the (two-level) logical Hilbert space and the occupation-number basis,
    this rotation is equivalent to the standard "phase gate" $P_q(ϕ)$:

$$ P_q(ϕ) = e^{iϕ(I-Z_q)/2} $$

I is the identity operator; $Z_q$ is the Pauli-Z operator acting on qubit $q$.
In our case, $ϕ=ω_q T$.

_This_ unitary is _especially_ easy to implement in most frameworks,
    because phase gates are usually implemented "virtually",
    meaning no actual quantum control has to take place.
Instead, they only modify the phases with which subsequent pulses are applied.
In our case, the only subsequent pulses are those used in measurement.

> From what I gather,
>   this is the solution that IBM uses to implement their traditional gate-based VQEs.
> I don't yet have a complete picture, though.
> First of all, coupling strengths $g_{pq}$ are _not_ normally zero.
> This method should certainly work to zeroth order when $g_{pq}, δ_q \ll ω_q$,
>     but a more thorough analysis is required to account for first-order error.
> Second, I don't understand yet what effect a virtual phase gate has on states outside the logical Hilbert space -
>     does $P_q$ actually approximate $u_q$, or does this method incur major leakage error?

#### Observable Rotation

Up to now, we have had to perform the frame rotation in the physical Hilbert space,
    and then project onto the logical Hilbert space.
Let us now restrict ourselves to the "Projected" strategy,
    and consider the special case where $[Π, H_0]=0$.
(Note that the model transmon Hamiltonian from above, with $g_{pq}≃0$, satisfies this condition.)
Then, defining $H_Π ≡ Π H_0 Π$ as the restriction of $H_0$ to the logical Hilbert space,
    the following relation can be shown:

$$ Π e^{iT H_0} = e^{iT H_Π} Π $$

This means we can project first, then perform the frame rotation in the _logical_ Hilbert space.
Writing out our energy function:

$$ E = ⟨ψ| Π^\dagger e^{-iT H_Π} \hat O e^{iT H_Π} Π |ψ⟩ $$

We can write this in a more familiar way as follows:

$$ \begin{align}
E &= ⟨Ψ|\tilde O_T|Ψ⟩ \\
|Ψ⟩ &≡ Π|ψ⟩ \\
\tilde O_T &≡ e^{-iT H_Π} \hat O e^{iT H_Π}
\end{align} $$

Now we are taking our logical state to be the physical state in the _lab frame_,
    projected to the logical Hilbert space,
    and we have applied the frame rotation to the _observable_.
Technically, we are considering the operator $\hat O$ to be "at rest" in the rotating frame,
    and the rotated operator $\tilde O_T$ is the same operator transformed into the lab frame.

> Calculating $\tilde O_T$ is done "offline", on the classical side.
> It is only experimentally tractable when $H_Π$ contains no coupling terms.
> On the one hand, the condition $[Π, H_0]=0$ is _probably_ only satisfied when there are indeed no coupling terms.
> On the other, there usually _are_ coupling terms.
> A more thorough analysis is required to understand the error accumulated from neglecting them.

NOTE: This strategy should work just as well for the "Normalized" projection strategy,
    as long as $[Π, H_0]=0$.
But the "Biased" projection strategy uses a non-Hermitian projector $Π_B$,
    meaning the operator $H_{Π_B} = Π_B H_0 Π_B$ is _not_ a restriction to the logical Hilbert space.
