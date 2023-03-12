# Energy Gradients
##### 3/12/2023

In another note (`timeevolution.md`), we have discussed algorithms to simulate time evolution of a quantum-computational state under some controllable drive Hamiltonian.
The next step in ctrl-VQE is to optimize the parameters - let's call them $\bar x$ - of the drive Hamiltonian, in order to minimize the expectation value $E$ of an observable, typically a molecular Hamiltonian.
To that end, we want to characterize the gradient vector $\frac{\partial E}{\partial \bar x}$.

### Energy Functions

The energy function we are going to minimize in this note may be wriiten explicitly as:

$$ E = ⟨\psi|\hat O|\psi⟩ $$

The observable $\hat O$ is typically a fermionic molecular Hamiltonian, originally written as a qubit operator in the "computational subspace" but then projected up into a bosonic Hilbert space.

The wavefunction $|\psi⟩$ is the time-evolved state we studied in `timeevolution.md`.
Please see that note for the full treatment; here I will present the bare minimum needed to calculate the gradient:

#### Time Evolution, Abridged

We start with a reference state $|\psi_{\rm ref}⟩$, and evolve it under a Hamiltonian consisting of both static ($H_0$) and drive ($\hat V(t)$) terms, for a total time $T$.
We control the drive terms, so all $\bar x$ dependence resides in $\hat V$.
In the interaction picture, the Hamiltonian can be written as:

$$ \tilde H(t) = e^{iH_0t} \hat V(t) e^{-iH_0t} $$

NOTE: Tunable couplings such as the ones considered by our collaborators at NIST would simply be categorized as drive terms, with a different drive operator from the ones we have normally been considering.
The formalism presented in this note will accommodate these systems without issue.

Evolution under the above Hamiltonian is given by:

$$ \tilde U = \exp\left(-i \int_0^T \tilde H(t) dt\right) $$

The integral is converted into a sum by:

$$ \int_0^T \tilde H(t) dt \approx \sum_{j=0}^{r} \tilde H(t_j) \tau_j $$

Using the trapezoidal rule for time integration, we define a fixed time-step $\tau\equiv T/r$, so that $t_j=j\tau$ and:

$$ \tau_j \equiv \left\{\begin{array}{ll}
    \tau    & 0 < j < r \\
    \tau/2  & j \in \{ 0, r \}
\end{array}\right.$$

Trotterizing our approximation to the integral, we may now expand $\tilde U$ as:

$$ \tilde U = \prod_{j=0}^{r} e^{-i \tilde H(t_j) \tau_j} $$

NOTE: The product places higher times on the left, to respect time-ordering.
I don't know the propert notation to write that explicitly.

Now, because the rotation into the interaction frame is unitary, we can bring all the $e^{iH_0t}$ factors out of the exponential.
Substituting in our definition of $t_j$, adjacent factors $e^{-iH_0t_j}e^{-H_0t_{j-1}}$ can be combined into the time-independent rotation $e^{-iH_0\tau}$, and the final factor $e^{-iH_0t_0}$ is identity:

$$ \tilde U = e^{i H_0 T} \prod_{j=1}^{r} \left( e^{-i \hat V(t_j) \tau_j} e^{-i H_0 \tau} \right) e^{-i \hat V(t_0) \tau_0} $$

The evolved wavefunction $|\psi⟩$ is (switching back to the lab frame):

$$ \begin{aligned} |\psi⟩ &\equiv e^{-i H_0 T} \tilde U |\psi_0⟩ \\
    &= \prod_{j=1}^{r} \left( e^{-i \hat V(t_j) \tau_j} e^{-i H_0 \tau} \right) e^{-i \hat V(t_0) \tau_0} |\psi_{\rm ref}⟩
\end{aligned} $$

NOTE: In practice, in order for a zero-pulse drive to result in a static energy, the static molecular Hamiltonian of interest would be written as an interaction operator $\tilde O$.
In that case, the lab operator $\hat O$ would be time-dependent, including factors of $e^{\pm iH_0 T}$.
As long as evolution time $T$ is considered constant, this has no bearing on the gradient calculations.

All dependence on control parameters $\bar x$ is confined to the $\hat V(t_j)$ operators.

#### Alternative Energy Functions

In the near future (but not in this note), we will also consider the following more experimentally relevant functionals:
- Projected: $ E = ⟨\psi| \Pi \hat O \Pi |\psi⟩ $
- Normalized: $ E = ⟨\psi| \Pi \hat O \Pi |\psi⟩ / ⟨\psi|\Pi|\psi⟩ $

The operator $\Pi$ projects out any components of $|\psi⟩$ outside the computational subspace, ie. "leaked states", representing measurements on $|\psi⟩$ which _must_ result in either 0 or 1.
In fact, the projected energy function introduces no new dependence on $\bar x$, so the results from this note apply to it perfectly well.
The normalized energy function will involve more calculus.

### Energy Gradients

Now we wish to consider partial derivatives of the energy with respect to individual control parameters $x$.

$$ \partial_x E = ⟨\psi|\hat O \left(\partial_x |\psi⟩\right) + {\rm h.t.} $$

To break down the braket, let us define the following kets:

$$ \begin{aligned}
|\lambda_i⟩ &\equiv \prod_{j=i+1}^{r} \left( e^{i H_0 \tau} e^{i \hat V(t_j) \tau_j} \right) \hat O |\psi⟩ \\
|\psi_i⟩ &\equiv \prod_{j=1}^{i} \left( e^{-i \hat V(t_j) \tau_j} e^{-i H_0 \tau} \right) e^{-i \hat V(t_0) \tau_0} |\psi_{\rm ref}⟩
\end{aligned} $$

NOTE: In $|\lambda_i⟩$ only, the product places larger times on the right.

Note that $|\psi_r⟩=|\psi⟩$, $|\psi_0⟩ \ne |\psi_{\rm ref}⟩$, and $E=⟨\lambda_i|\psi_i⟩$ for all $i$.
More importantly:

$$ ⟨\psi|\hat O \left(\partial_x |\psi⟩\right) = \sum_{i=0}^r ⟨\lambda_i| \left( -i \tau_i \partial_x \hat V(t_i) \right) |\psi_i⟩ $$

This is the most important equation in this note.
Pretend there is a box around it.

Next, we're going to break this down for two special forms of the drive operator.