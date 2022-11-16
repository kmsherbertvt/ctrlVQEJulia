# Time Evolution
##### 11/15/2022

The idea of ctrl-VQE is to directly optimize the experimental parameters used when operating a quantum computer.
Typically (at least for a transmon-based device) these parameters correspond to a time-dependent EM wave (eg. a microwave laser).

In order to study ctrl-VQE, we can simulate this time-dependent system _in silico_.
The key objective, then, is to solve for the quantum state of the quantum computer after being subjected to some time-dependent fields,
and then ultimately understand how to vary those time-dependent fields to achieve a quantum state best satisfying some property
  (eg. minimizing the energy of a molecular Hamiltonian).
This notebook focuses on the former problem, for transmon qubits.

### Glossary of Symbols
- $n=$ the number of qubits in a device
- $m=$ the number of modes (or states) considered for each qubit
- $N=$ the size of the Hilbert space considered, $N≡m^n$
- $ω=$ qubit resonance frequencies
- $δ=$ qubit anharmonicities
- $g=$ qubit coupling constants
- $Ω=$ EM modulating wave amplitude
- $ν=$ EM carrier wave frequencies
- $T=$ total duration of EM wave
- $r=$ number of steps in a discretized time evolution
- $Δt=$ length of step in a discretized time evolution

Throughout, assume fundamental constants like $\hbar$ are unity, so that energy and frequency are interchangeable.

### The Model
Our problem is one of a time-dependent perturbation:
  the structure of the quantum computer itself confers a "static" Hamiltonian $\hat H_0$,
  while the controllable EM wave confers a time-dependent "drive" Hamiltonian $\hat V(t)$.
  
$$ \hat H = \hat H_0 + \hat V(t) $$

#### Time-Independent Hamiltonian
Associate with each of $n$ transmons a set of vibrational modes $|0⟩$, $|1⟩$, $|2⟩$, etc.
If we approximate the vibrational modes of transmon $q$ as residing in a harmonic potential, we can say they are spaced apart by $ω_q$ energy units.
If we go a step further and assume an ever-so-slightly _anharmonic_ potential,
  we can say that the spacing of each energy level is reduced by a small amount $δ_q$ at each level.
Finally, we will say that some transmons are coupled,
  so that the energies of transmons $p$ and $q$ will be ever-so-slightly modified by an amount $g_{pq}$.

We write the (time-independent) Hamiltonian of the quantum computer as:

$$ \hat H_0 = \sum_q ( ω_q \hat a_q^† \hat a_q - \frac{δ_q}{2} \hat a_q^† \hat a_q^† \hat a_q \hat a_q ) + \sum_{⟨pq⟩} g_{pq} (\hat a_p^† \hat a_q + \hat a_q^† \hat a_p )  $$

The operators $\hat a^†$, $\hat a$ are bosonic creation and annihilation operators, (note that transmons $p$ and $q$ are _distinguishable_).
The notation $⟨pq⟩$ denotes a sum over unique pairs.

#### Time-Dependent Hamiltonian
Consider the controllable EM wave as a pure sine wave (eg. a laser) with tunable frequency and amplitude.
The energy perturbation induced by this EM wave will depend on the field's polarization and the dipole moment of each transmon,
  but we can abstract away those details and say that each transmon $q$ experiences an effective energy amplitude $Ω_q$.
Further, we'll say that multiple EM waves can be applied to individually address each transmon.
Therefore, the overall drive Hamiltonian (drawing on results from quantum field theory) is:

$$ \hat V(t) = \sum_q Ω_q (e^{iν_q t} \hat a_q + e^{-iν_q t} \hat a_q^†) $$

The parameters $Ω_q$ and $ν_q$ are the parameters that can be controlled.
They are generally time-dependent, and in particular we will assume $Ω_q(t<0)=Ω_q(t>T)=0$ for some "pulse duration" $T$.

### Time Evolution Basics
Represent as $|ψ⟩$ the quantum state of our quantum computer, which exists in the Hilbert space spanned by the product of all individual transmon spaces.
In the Schrödinger picture, the time evolution of $|ψ⟩$ is given by:

$$ i \frac{∂}{∂t} |ψ⟩ = \hat H |ψ⟩ $$

By truncating each transmon to just $m$ modes, we can represent $|ψ⟩$ and $\hat H$ at any given time
  with a finite $N$-dimensional vector and a finite $N×N$ matrix, where $N≡m^n$.
Thus, Schrödinger's equation is a linear first-order system of (many) differential equations,
  and in principle it can be integrated numerically.
  
Alternatively, we can solve this system symbolically, by momentarily pretending that $|ψ⟩$ and $\hat H$ are scalar...

$$ |ψ(t)⟩ = e^{ -i \int_0^t \hat H dt } |ψ(0)⟩ $$

When $H$ is time-independent, solving the time-evolution of $|ψ⟩$ just reduces to characterizing the matrix-exponential operator $e^{-it \hat H}$.
When $H$ is time-dependent, a robust analytical treatment of the strange exponentiated-matrix-integral object
  requires the rather arcane "time-ordering operator".
Its effect on a numerical treatment is somewhat more intuitive:
  by discretizing time, the integral becomes a sum, and the exponential-of-sums becomes a time-ordered product-of-exponentials:

$$\begin{array}{rl}
  |ψ(t)⟩ &= \exp \left\[ -i Δt \displaystyle\sum_{i=1}^r \hat H_i \right\] |ψ(0)⟩ \\
         &= e^{ -i Δt \hat H_r } … e^{ -i Δt \hat H_2 } · e^{ -i Δt \hat H_1 } |ψ(0)⟩
\end{array}$$

The first step, turning the integral into a sum, can be done in a number of different ways which will be discussed a little later on.

The last step, turning an exponential-of-sums into a product-of-exponentials, is known as Trotterization,
  and it is only exact when all the terms commute.
This will not generally be the case!
There will be some error which tends to scale by $O(Δt^2)$, so the error vanishes as you increase $r$.
In fact there are more complicated ways to Trotterize which have better error scaling,
  but I don't actually know how they interact with the time-ordering operator,
  so they're best omitted from this discussion.

#### Interaction Picture
For time-dependent perturbations, it is helpful to adopt instead the "interaction picture", which focuses on a slightly-modified equation:

$$ i \frac{∂}{∂t} |ψ_I⟩ = \hat V_I |ψ_I⟩ $$

The point of this equation is to abstract away the time-independent part of $\hat H = \hat H_0 + \hat V(t)$.
Since $H_0$ is time-independent, it makes sense to work with the static time-evolution operator $e^{-it \hat H_0}$.
Now,

$$\begin{array}{rl}
   |ψ_I(t)⟩ &≡ e^{it \hat H_0} |ψ(t)⟩ \\
\hat V_I(t) &≡ e^{it \hat H_0} \hat V(t) e^{-it \hat H_0}
\end{array}$$

Vaguely speaking, $|ψ_I⟩$ represents the state $|ψ⟩$ as if it were not subject to the time-independent $\hat H_0$.
More precisely, and more usefully, you can solve for $|ψ⟩$ simply by evolving $|ψ_I⟩$ in time subject _only_ to the time-independent $\hat H_0$.

We do this because, at least for oscillating field problems, the system of differential equations tends to look simpler in the interaction picture.
Therefore, we will synthesize the results of this section and the previous to form our two main time-evolution methods:

##### Time Evolution by Solving a System of Linear Differential Equations

1. Truncate each transmon up to $m$ modes.
2. Represent the quantum state $|ψ_I⟩$ with a finite vector $\vec ψ$ and the operator $\hat V_I$ with a finite matrix $\mathbf{V}$.
3. Solve the following system of differential equations:

$$ \dot{\vec ψ} = -i \mathbf{V}(t) · \vec ψ $$

##### Time Evolution by Trotterization

1. Truncate each transmon up to $m$ modes.
2. Represent the quantum state $|ψ_I⟩$ with a finite vector $\vec ψ$ and the operator $\hat V_I$ with a finite matrix $\mathbf{V}$.
3. Evaluate the product:

$$ \vec ψ(t) = e^{ -i Δt \mathbf{V_r} } … e^{ -i Δt \mathbf{V_2} } · e^{ -i Δt \mathbf{V_1} } \vec ψ(0) $$

#### Trapezoidal Integration
Recall that the Trotterization method required discretizing an integral into a sum, and doing so introduced subscripts into our operators:

$$ \int_0^T \hat H(t) dt → \sum_{i=1}^r \hat H_i Δt $$

We need to decide what $\hat H_i$ means!

Just think of the good-old Riemann sums you used when you were first learning what an "integral" even is.
We had _left_-handed sums, _right_-handed sums, and then a more symmetric _trapezoidal rule_.

The left-handed sum would evaluate a function at the "left-handed" time point, ie. the earlier time for each step.
In that case, $\hat H_i ≡ \hat H((i-1)Δt)$.
The right-handed sum uses the "right-handed" time point, producing the somewhat-more-elegant definition $\hat H_i ≡ \hat H(iΔt)$.
Finally, the trapezoidal rule _averages_ the function at both time points:

$$ \hat H_i ≡ \frac{1}{2} \left\[ \hat H((i-1)Δt) + \hat H(iΔt) \right] $$

This one tends to be quite a bit more accurate, so it's the one we're going to use.
(Like in Trotterization, there are higher-order formulae that peform significantly better with fewer $r$,
  but I don't know how they interact with the time-ordering operator in the next step,
  so I've elected to not pursue them.)

Now let's see what effect this has on our Trotterization method:

$$\begin{array}{rl}
  \vec ψ(t) &= e^{ -i Δt \mathbf{V_r} } … e^{ -i Δt \mathbf{V_2} } · e^{ -i Δt \mathbf{V_1} } \vec ψ(0) \\
            &= e^{ -i \frac{Δt}{2} \left\[\mathbf{V}(T)+\mathbf{V}(T-Δt)\right\] } … e^{ -i \frac{Δt}{2} \left\[\mathbf{V}(2Δt)+\mathbf{V}(Δt)\right\] } · e^{ -i \frac{Δt}{2} \left\[\mathbf{V}(Δt)+\mathbf{V}(0)\right\] } \vec ψ(0) \\
            &≈ \left( e^{ -i \frac{Δt}{2} \mathbf{V}(T) } · e^{ -i \frac{Δt}{2} \mathbf{V}(T-Δt) } \right) … \left( e^{ -i \frac{Δt}{2} \mathbf{V}(2Δt) } · e^{ -i \frac{Δt}{2} \mathbf{V}(Δt) } \right) · \left( e^{ -i \frac{Δt}{2} \mathbf{V}(Δt) } · e^{ -i \frac{Δt}{2} \mathbf{V}(0) } \right) \vec ψ(0) \\
            &= e^{ -i \frac{Δt}{2} \mathbf{V}(T) } · e^{ -i Δt \mathbf{V}(T-Δt) } … e^{ -i Δt \mathbf{V}(2Δt) } · e^{ -i Δt \mathbf{V}(Δt) } · e^{ -i \frac{Δt}{2} \mathbf{V}(0) } \vec ψ(0) \\
\end{array}$$

Unpacking the new terms in the exponentials via Trotterization, and being careful to keep all terms time-ordered,
  we discover that the trapezoidal rule asks us to perform a time evolution at all $r+1$ time points in our discrete grid,
  but to treat the very first and last points as though the time-spacing $Δt$ were halved.
Note that, while this derivation seems to incur some extra error from factoring the left- and right-handed contributions at each time-step,
  it does not affect the error-_scaling_ already implicit in the original Trotterization formula.

Please note that **all Trotter evolution methods implemented in code use the trapezoidal rule**.
However, for the remainder of these notes, we will omit explicit reference to it, to keep notation a little simpler.


### Evolution Methods

The code includes a number of distinct evolution methods; this section attempts to digest each one mathematically.

#### `ODE`

`ODE` is short for "Ordinary Differential Equation" and at present,
  this method is the only performing time evolution by solving a system of linear differential equations.
  
$$ \dot{\vec ψ} = -i \mathbf{V}(t) · \vec ψ $$

There are in fact _countless_ ways to numerically solve a system of differential equations.
I don't know very much about them, which is why _this_ method gives Julia the problem and lets _her_ decide how to solve it.
In other words, you should consider this method a black box.
An unsophisticated benchmarking analysis indicates that Julia selects a fairly low-overhead $O(N^3)$ algorithm,
  but there is presently no particular guarantee that scaling behavior is consistent for larger systems.

#### `Direct`

`Direct` indicates "Direct Exponentiation" at each time step.
It takes the Trotter time-evolution formula very literally.

$$ \vec ψ(t) = e^{ -i Δt \mathbf{V_r} } … e^{ -i Δt \mathbf{V_2} } · e^{ -i Δt \mathbf{V_1} } \vec ψ(0) $$

Each of these matrix exponentials is individually calculated, using Julia's default matrix exponentiation.
_Assuming_ Julia is implemented the way _I_ would have implemented it, it is equivalent to the following:

- For each time point $t$, diagonalize $\mathbf{V_i} = \mathbf{U_{V_i}} · \mathbf{Λ_{V_i}} · \mathbf{U_{V_i}^†}$,
    where $\mathbf{U_{V_i}}$ is a unitary matrix of eigenvectors and $\mathbf{Λ_{V_i}}$ is a diagonal matrix of eigenvalues.
  This step is the bottleneck, costing $O(N^3)$ runtime.
- Next, calculate $e^{-iΔt\mathbf{Λ_{V_i}}}$.
  If you like, this is the desired matrix exponential _in the eigenbasis_.
  Because $\mathbf{Λ_{V_i}}$ is diagonal, this step costs only $O(N)$ runtime.
- Finally, return to the original basis by conjugating with $\mathbf{U_{V_i}}$, which involves an $O(N^3)$ matrix multiplication:

$$ \mathbf{E_i} = \mathbf{U_{V_i}} · e^{-iΔt\mathbf{Λ_{V_i}}} · \mathbf{U_{V_i}^†} $$

We have so far accumulated a total runtime of $O(rN^3)$ and the following formula:

$$ \vec ψ(t) = \mathbf{E_r} … \mathbf{E_2} · \mathbf{E_1} \vec ψ(0) $$

We contract from right to left, performing a sequence of matrix-vector calculations (each $O(N^2)$).
Thus, the total runtime remains $O(rN^3)$.

There is one more detail to mention, which is how $\mathbf{V_i}$ itself is computed prior to Julia's matrix exponentiation.
Recall that $\mathbf{V_i}$ is a matrix representation of the interaction-picture Hamiltonian:

$$ \hat V_I(t) ≡ e^{it \hat H_0} \hat V(t) e^{-it \hat H_0} $$

When working in the _device basis_, where each amplitude in $\vec ψ$ corresponds to an eigenvector of $\hat H_0$,
  the conjugating factor $e^{it \hat H_0}$ is diagonal; it's really just all the eigenvalues of $\hat H_0$ mulitiplied by $it$.
Therefore, once you have a matrix representation of $\hat V(t)$,
  obtaining $\mathbf{V_i}$ is very efficient: you just mulitply each element by an easily-determined phase factor.
_However_, if we are in the device basis, constructing the matrix representation of $\hat V(t)$ is not so easy.
The most efficient strategy seems to be to retain in memory an $N×N$ matrix representation $\mathbf{a_q}$ of each $\hat a_q$ in the device-basis.
Then the matrix representation of $\hat V(t)$ is $\sum_q \mathbf{V^{(q)}}$, where

$$ \mathbf{V^{(q)}} = Ω_q (e^{iν_q t} \mathbf{a_q} + e^{-iν_q t} \mathbf{a_q^†}) $$

is just a simple matrix sum over $N×N$ operators.
Thus, obtaining each $\mathbf{V_i}$ is an $O(nN^2)$ operation.
While there may be a more efficient construction of $\hat V(t)$ in the qubit basis (see the `Rotate` method!),
  the conjugation would then require at least one $O(N^3)$ matrix multiplication step,
  so this is probably the best implementation for this particular evolution method.

Each $\mathbf{E_i}$ can be treated successively so memory requirements can in principle be bounded by $O(N^2)$.
However, maintaining each $\hat a_q$ in memory increases this to $O(nN^2)$.


#### `Lanczos`

This method is almost entirely identical to `Direct`, except that it avoids calculating the matrix exponential $\mathbf{E_i}$ explicitly.
Instead, a "Lanczos" scheme is used to construct a "Krylov" sub-space of the generating matrix $\mathbf{V_i}$,
  and this Krylov subspace is used to compute the action of the matrix exponential $\mathbf{E_i}$ on $\vec ψ$ directly.

I'm not very familiar with the Lanczos scheme _or_ the Krylov subspace, but from what I understand,
  the runtime of a proper implementation _should_ be bounded by $O(rN^2)$, but with a very steep overhead.
Memory requirements are probably also a very steep overhead,
  but should scale asymptotically as $O(nN^2)$, just as in `Direct`.
The studious reader can research the method used [here](https://jutho.github.io/KrylovKit.jl/stable/man/matfun/#KrylovKit.exponentiate).

#### `Rotate`

`Rotate` indicates "Rotating back and forth between the device and drive bases at each time step".
It is a Trotterized method, but I don't want to work in matrix representation just yet:

$$ |ψ(t)⟩ = e^{ -i Δt \hat V_I(t_r) } … e^{ -i Δt \hat V_I(t_2) } · e^{ -i Δt \hat V_I(t_1) } |ψ(0)⟩ $$

but we will break each time step into device and drive components.
It will turn out that the device components can be treated entirely _independent_ of time, and the drive components can be treated _qubit-wise_,
  significantly reducing the number and size of computations.

Recall that the interaction picture Hamiltonian is

$$ \hat V_I(t_i) ≡ e^{it_i \hat H_0} \hat V(t_i) e^{-it_i \hat H_0} $$

Conjugation by the unitary $e^{it_i \hat H_0}$ is in essence representing the drive Hamiltonian $\hat V(t_i)$ in an alternate, time-dependent basis.
Substituting this into our Trotterized time evolution, the conjugating terms (because they are unitary) can pass _out_ of the enveloping exponential:

$$\begin{array}{rl}
  |ψ(t)⟩ &= \exp\left\[ e^{it_r \hat H_0} (-iΔt\hat V(t_r)) e^{-it_r \hat H_0}\right\] … \exp\left\[ e^{it_2 \hat H_0} (-iΔt\hat V(t_2)) e^{-it_2 \hat H_0}\right\] · \exp\left\[ e^{it_1 \hat H_0} (-iΔt\hat V(t_1)) e^{-it_1 \hat H_0}\right\] |ψ(0)⟩ \\
         &= e^{it_r \hat H_0} · \exp\left\[ -i Δt \hat V(t_r) \right\] · \left( e^{-it_r \hat H_0} · e^{it_{r-1} \hat H_0} \right) · \exp\left\[ -i Δt \hat V(t_{r-1}) \right\] … \exp\left\[ -i Δt \hat V(t_2) \right\] · \left( e^{-it_2 \hat H_0} · e^{it_1 \hat H_0} \right) · \exp\left\[ -i Δt \hat V(t_1) \right\] · e^{-it_1 \hat H_0} |ψ(0)⟩ \\
         &= e^{it_r \hat H_0} · \exp\left\[ -i Δt \hat V(t_r) \right\] · e^{-i(t_r - t_{r-1}) \hat H_0} · \exp\left\[ -i Δt \hat V(t_{r-1}) \right\] … \exp\left\[ -i Δt \hat V(t_2) \right\] · e^{-i(t_2 - t_1) \hat H_0} · \exp\left\[ -i Δt \hat V(t_1) \right\] · e^{-it_1 \hat H_0} |ψ(0)⟩ \\
         &= e^{it_r \hat H_0} · \exp\left\[ -i Δt \hat V(t_r) \right\] · e^{-iΔt \hat H_0} · \exp\left\[ -i Δt \hat V(t_{r-1}) \right\] … \exp\left\[ -i Δt \hat V(t_2) \right\] · e^{-iΔt \hat H_0} · \exp\left\[ -i Δt \hat V(t_1) \right\] · e^{-it_1 \hat H_0} |ψ(0)⟩
\end{array}$$

Aside from the very first and last time steps, the device Hamiltonian $H_0$ is always treated with the same "ligand" operator $e^{-iΔt \hat H_0}$.
Let $\mathbf{H_0}$ be the matrix representation of our device Hamiltonian in the _qubit basis_.
The matrix exponential, as described in the `Direct` method, is computed via the following steps:
- Diagonalize $\mathbf{H_0} = \mathbf{U_0} · \mathbf{Λ_0} · \mathbf{U_0^†}$,
    where $\mathbf{U_0}$ is a unitary matrix of eigenvectors and $\mathbf{Λ_0}$ is a diagonal matrix of eigenvalues.
  This step is the bottleneck, costing $O(N^3)$ runtime.
- Next, calculate $e^{-iΔt\mathbf{Λ_0}}$.
  This is in fact the ligand operator computed in the in the _device basis_.
  Because $\mathbf{Λ_0}$ is diagonal, this step costs only $O(N)$ runtime.
- Finally, rotate back to the qubit basis by conjugating with $\mathbf{U_0}$.

$$ \mathbf{L_0} = \mathbf{U_0} · e^{-iΔt\mathbf{Λ_0}} · \mathbf{U_0^†} $$

Like in the `Direct` method, calculating this matrix exponential has a runtime of $O(N^3)$.
_Unlike_ the `Direct` method, this is the _only_ time we need to perform an $O(N^3)$ operation!

We also need to evolve $H_0$ for the first and last points.
But in the trapezoidal rule, the first time point is actually when $t=0$, so that operator is actually just identity.
For the last time point ( $t=T$ ), we can re-use the same factorization $\mathbf{H_0} = \mathbf{U_0} · \mathbf{Λ_0} · \mathbf{U_0}^†$.
Therefore, after applying the last drive, we apply the matrix $\mathbf{U_0}^†$ to rotate back to the device basis ( $O(N^2)$ ),
  then we apply the final exponential $e^{iT \mathbf{Λ_0}}$, which has only $O(N)$ runtime.
Now we already _have_ our final statevector in the device basis, but we can choose to apply $\mathbf{U_0}$ to rotate back to the qubit basis.

Now let's look more closely at the drive Hamiltonian.

$$ \hat V(t_i) = \sum_q Ω_q (e^{iν_q t_i} \hat a_q + e^{-iν_q t_i} \hat a_q^†) $$

Since we're in the qubit basis, it makes sense to represent each annihilation operator $\hat a_q$ with just the one $m×m$ operator $\mathbf{a}$.
Then our qubit-wise drive Hamiltonian is represented by:

$$ \mathbf{V^{(q)}}(t_i) = Ω_q (e^{iν_q t} \mathbf{a} + e^{-iν_q t} \mathbf{a^†}) $$

Let $\mathbf{D^{(q)}}(t_i)$ be the matrix exponential $e^{-iΔt \mathbf{V^{(q)}}(t_i)}$.
The full-body operator $e^{-i Δt \hat V(t_i)}$ is therefore representable as the matrix product operator $\mathbf{D}(t_i)$:

$$ e^{-i Δt \hat V} → \mathbf{D} ≡ \mathbf{D^{(q)}} ⊗ … \mathbf{D^{(2)}} ⊗ \mathbf{D^{(1)}} $$

Now we are left with two questions:
1. How do we apply the matrix product operator $\mathbf{D}(t_i)$ to the statevector $\vec ψ$?
2. How do we calculate each $m×m$ matrix exponential $\mathbf{D^{(q)}}$?

Each of these questions is answered in multiple ways in the code, and is therefore worthy of its own subsection.

##### Applying the matrix product operator $\mathbf{D}(t)$ to the statevector $\vec ψ$
The conceptually-easiest approach is to apply the Kronecker product directly to calculate an $N×N$ matrix,
  which has a runtime of $O(nN^2)$,
  and then apply this matrix directly to $\vec ψ$, which has a runtime of $O(N^2)$.
So far we've discussed applying the ligand matrix $\mathbf{L}$, which requires a runtime of $O(N^2)$ at each time step,
  or $O(rN^2)$ overall assuming $\mathbf{L}$ is pre-computed.
Applying $\mathbf{D}$ with this method increases the asymptotic complexity to $O(rnN^2)$.
This is a significant improvement over the `Direct` method but it could be even better.

The trickier-but-more-efficient approach is to reinterpret $\vec ψ$ as an $n$-body tensor $ψ_{i_n,…i_2,i_1}$,
  then apply $\mathbf{D}$ with a tensor contraction:

$$ \mathbf{D}ψ_{i_n,…i_2,i_1} = \mathbf{D_{i_n,j_n}^{(q)}} … \mathbf{D_{i_2,j_2}^{(2)}} \mathbf{D_{i_1,j_1}^{(1)}} ψ_{j_n,…j_2,j_1} $$

In principle, this contraction can be computed with a runtime of $O(nN)$.
As long as it is implemented well, this results in a runtime bounded by applying the ligand operator: $O(rN^2)$.
This is, I think, as good as it gets, unless sparsity in $\mathbf{L}$ can be exploited.

In the code, these two methods can be selected via the keyword argument `qubitapplymode`.

##### Calculating the matrix exponential $\mathbf{D^{(q)}}(t)$
The most straightforward approach, and the one taken by the `Rotate` method, is to directly exponentiate.
As described in the `Direct` method, this _presumably_ involves diagonalizing $\mathbf{V^{(q)}}$, which has a runtime of $O(m^3)$.
Repeating for each $q$ incurs a runtime of $O(nm^3)$ at each timestep,
  negligible next to the $O(N^2)=O(m^{2n})$ runtime to apply the ligand operator $\mathbf{L}$.

This is perfectly fine and it (spoilers) don't actually get better.
But, committing to the idea of performing all diagonalizations just once during the whole algorithm, we could _try_ to improve it;
  this is what the `Prediag` method does.

Consider the diagonalization of $\mathbf{V^{(q)}}$.

$$ \mathbf{V^{(q)}}(t) = \mathbf{U_q}(t) \mathbf{Λ_q}(t) \mathbf{U_q^†}(t) $$

We'd really _like_ to be able to compute $\mathbf{D^{(q)}} = e^{-iΔt \mathbf{V^{(q_)}}}$ as $\mathbf{U_q} e^{-iΔt \mathbf{Λ_q}(t)} \mathbf{U_q^†}$,
  where $\mathbf{U_q}$ is independent of time.
The diagonal $\mathbf{Λ_q}(t)$ would presumably be calculable in $O(m)$ runtime, as would its matrix exponential,
  and then conjugation with $\mathbf{U_q}$ would give $\mathbf{D^{(q)}}$ in just $O(m^2)$ time.
The problem is that $\mathbf{U_q}$ is _not_ independent of time, and it takes $O(m^3)$ runtime to calculate it.
(Again, it's not _actually_ a problem, but this section is trying to solve it anyways...)

What if we diagonalized not $\mathbf{V^{q}}$, but each individual _term_ in $\mathbf{V^{q}}$?
Then we can "Trotterize" to evaluate the matrix exponential.
Perhaps this will remove time-dependence from the eigenvectors?

Yes, it will work, but careful!
We do not want to evaluate the matrix exponential of a matrix like $e^{-iz Δt \mathbf{a}}$, because $\mathbf{a}$ _is not Hermitian_.
Moreover, $\mathbf{a}$ and $\mathbf{a^†}$ do not commute, so the $\mathbf{D^{(q)}}$ resulting from Trotterization will _not be unitary_:
  our time evolution will not preserve the norm of the statevector, resulting in numerical instability.
Therefore, we will re-cast our problem into a different pair of operators:

$$\begin{array}{rl}
  \hat Q &≡ \hat a + \hat a^† \\
  \hat P &≡ i(\hat a - \hat a^†)
\end{array}$$

These operators are the so-called "canonical" operators,
  which have a rich history in the development of both classical and quantum mechanics
  ( $\hat Q$ is known as a "canonical coordinate" and $\hat P$ as its "canonical momentum").
They are useful to us here and now because they are Hermitian.
We will write their matrix representations (truncated to $m$ modes) as $\mathbf{Q}$ and $\mathbf{P}$,
  their eigenvector matrices as $\mathbf{U_Q}$ and $\mathbf{U_P}$,
  and the diagonal matrix of eigenvalues as $\mathbf{Λ}$ (they share the same eigenvalues!).

Decomposing our qubit-wise drive Hamiltonian into a sum over _these_ operators:

$$\begin{array}{rl}
  \mathbf{V^{(q)}}(t) &= x_q(t) \mathbf{Q} + y_q(t) \mathbf{P} \\
               x_q(t) &≡ \Re ( Ω_q e^{iν_q t} ) \\
               y_q(t) &≡ \Im ( Ω_q e^{iν_q t} ) \\
\mathbf{D^{(q)}}(t) &≡ e^{-iΔt \mathbf{V^{(q)}} } \\
                    &= e^{-iΔt (x_q(t) \mathbf{Q} + y_q(t) \mathbf{P}) } \\
                    &≈ e^{-i x_q(t) Δt \mathbf{Q}} · e^{-i y_q(t) Δt \mathbf{P}} \\
                    &= \mathbf{U_Q} · e^{-i x_q(t) Δt \mathbf{Λ}} · \mathbf{U_Q^†} · \mathbf{U_P} · e^{-i y_q(t) Δt \mathbf{Λ}} · \mathbf{U_P^†}
\end{array}$$

Here is a formula which permits us to calculate $\mathbf{D^{(q)}}(t)$ using pre-diagonalized matrices!
The _only_ trouble is that it seems to require at least two matrix-matrix multiplications at each time step, which have... $O(m^3)$ runtime.
(I say at least two because $\mathbf{U_Q^†} · \mathbf{U_P}$ can be pre-computed,
  and the diagonal factors are efficiently folded into an adjacent matrix.)
Yeah, diagonalizing a matrix has the same asymptotic complexity as matrix multiplication.
The direct exponentiation has a little overhead from requiring a matrix multiplication after solving the eigenvalue problem,
  but since this formula itself requires _two_ matrix multiplications, that balances out...
  
Oh, but can't we do matrix-vector multiplications instead..?
Alright, let's say you are using the efficient tensor contraction method to apply each of these matrix products one at a time.
You have to apply at least three contractions, each being $O(N)$ (or really $O(nN)$ if doing each qubit in parallel; it doesn't change the analysis).
So let's say that's $O(m^3+N)$ to do the matrix multiplication then contraction, or $O(N+N)$ to do more than one contraction.
Doing mulitple contractions is preferable when $N < m^3$.
Well, $N=m^n$.
In other words, when working with qubit-wise matrices,
  matrix-matrix multiplication is faster than matrix-vector multiplication except when $n=1$.
In other words, there is no reason whatsoever to use the `Prediag` method over `Rotate`.

One last note: the code permits a keyword argument `suzukiorder` which can be 0, 1, or 2.
`suzukiorder=2` uses one of those higher-order product formulae I've alluded to a couple times;
  we can definitely use it here because we're factoring operators corresponding to the same time; the time-ordering operator won't interfere.
Personally I like `suzukiorder=2` a lot because it's simple, symmetric, and drops Trotter error by another degree ( $O(Δt^2)→O(Δt^3)$ ),
  but the fact is since the Trotter time evolution has $O(Δt^2)$ error, it doesn't really make much difference.
The `suzukiorder=0` option is, however, a bit more interesting.
To be clear, there isn't any such thing as a 0-order product formula; this is just short-hand for a related technique:

The product formula $e^{\hat A + \hat B}≈e^{\hat A}·e^{\hat B}$ is actually exact when the commutator $\[\hat A, \hat B\]$ is zero.
More interestingly, the formula $e^{\hat A + \hat B}≈e^{\hat A}·e^{\hat B}·e^{-\[\hat A, \hat B\]/2}$ is exact when the commutator is _scalar_.
This is in fact true of our canonical operators $\hat Q$ and $\hat P$: $\[\hat Q, \hat P\] = -2i$.
So, simply including that commutator as an extra phase factor in the first-order product formula should be enough to make it exact!
This is what `suzukiorder=0` does.
To be clear, there's still no reason to bother since `Rotate` is just as exact and generally faster,
  and anyways you still have Trotter error from the time evolution formula, but it's _neat_.

Unfortunately, it turns out that, practically speaking, this _doesn't_ work, because, truncated to $m$ modes, the commutator $\[\mathbf{Q}, \mathbf{P}\]$ is _not quite_ scalar: it is a diagonal operator with $-2i$ in every entry except the _very last_ one.
Presumably this method will approximately work for significantly larger $m$, but I reiterate one last time that there is no point. ;)
