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
Associate with each transmon qubit $q$ a set of vibrational modes $|0⟩$, $|1⟩$, $|2⟩$, etc.
If we approximate these vibrational modes as residing in a harmonic potential, we can say they are spaced apart by $ω_q$ energy units.
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

By truncating each transmon to just $m$ modes, we can represent $|ψ⟩$ and $\hat H$ at any given time with a finite vector and a finite matrix.
Thus, Schrödinger's equation is a linear first-order system of (many) differential equations,
  and in principle it can be integrated numerically.
  
Alternatively, we can solve this system symbolically, by momentarily pretending that $|ψ⟩$ and $\hat H$ are scalar...

$$ |ψ(t)⟩ = e^{ -i \int_0^t \hat H dt } |ψ(0)⟩ $$

When $H$ is time-independent, solving the time-evolution of $|ψ⟩$ just reduces to characterizing the matrix-exponential operator

$$ \hat U(t) ≡ e^{-it \hat H} $$

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
Since $H_0$ is time-independent, it makes sense to define the matrix-exponential operator $\hat U_0(t)≡e^{-it \hat H_0}$.
Now,

$$\begin{array}{rl}
  |ψ_I(t)⟩ &≡ \hat U(-t) |ψ(t)⟩ \\
    V_I(t) &≡ \hat U(-t) \hat V(t) \hat U(t)
\end{array}$$

Roughly speaking, $|ψ_I⟩$ represents the state $|ψ⟩$ as if it were not subject to the time-independent $\hat H_0$.
Perhaps more to the point, you can solve for $|ψ⟩$ simply by evolving $|ψ_I⟩$ in time subject _only_ to the time-independent $\hat H_0$.

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

$$ \int \hat H(t) dt → \sum_{i=1}^r \hat H_i Δt $$

We need to decide what $\hat H_i$ means!

Just think of the good-old Riemann sums you used when you were first learning what an "integral" _is_.
We had _left_-handed sums, _right_-handed sums, and then a more symmetric _trapezoidal rule_.

The left-handed sum would evaluate a function at the "left-handed" time point, ie. the earlier time for each step.
In that case, $\hat H_i ≡ \hat H((i-1)Δt)$.
The right-handed sum uses the "right-handed" time point, producing the somewhat-more-elegant definition $\hat H_i ≡ \hat H(iΔt)$.
Finally, the trapezoidal rule _averages_ the function at both time points:

$$ \hat H_i ≡ \frac{1}{2} \left\[ \hat H((i-1)Δt) + \hat H(iΔt) \right] $$

This one tends to be significantly more accurate, so it's the one we're going to use.
(Like in Trotterization, there are higher-order formulae that peform significantly better with fewer $r$,
  but I don't know how they interact with the time-ordering operator in the next step,
  so I've elected to not pursue them.)

Now let's see what effect this has on our Trotterization method:

$$\begin{array}{rl}
  \vec ψ(t) &= e^{ -i Δt \mathbf{V_r} } … e^{ -i Δt \mathbf{V_2} } · e^{ -i Δt \mathbf{V_1} } \vec ψ(0) \\
            &= e^{ -i \frac{Δt}{2} \left\[\mathbf{V}(T)+\mathbf{V}(T-Δt)\right\] } … e^{ -i \frac{Δt}{2} \left\[\mathbf{V}(2Δt)+\mathbf{V}(Δt)\right\] } · e^{ -i \frac{Δt}{2} \left\[\mathbf{V}(Δt)+\mathbf{V}(0)\right\] } \vec ψ(0) \\
            &= \left( e^{ -i \frac{Δt}{2} \mathbf{V}(T) } · e^{ -i \frac{Δt}{2} \mathbf{V}(T-Δt) } \right) … \left( e^{ -i \frac{Δt}{2} \mathbf{V}(2Δt) } · e^{ -i \frac{Δt}{2} \mathbf{V}(Δt) } \right) · \left( e^{ -i \frac{Δt}{2} \mathbf{V}(Δt) } · e^{ -i \frac{Δt}{2} \mathbf{V}(0) } \right) \vec ψ(0) \\
            &= e^{ -i \frac{Δt}{2} \mathbf{V}(T) } · e^{ -i Δt \mathbf{V}(T-Δt) } … e^{ -i Δt \mathbf{V}(2Δt) } · e^{ -i Δt \mathbf{V}(Δt) } · e^{ -i \frac{Δt}{2} \mathbf{V}(0) } \vec ψ(0) \\
\end{array}$$

Unpacking the new terms in the exponentials via Trotterization, and being careful to keep all terms time-ordered,
  we discover that the trapezoidal rule asks us to perform a time evolution at all $r+1$ time points in our discrete grid,
  but to treat the very first and last points as though the time-spacing $Δt$ were halved.