# Energy Gradients
##### 3/12/2023

In another note (`timeevolution.md`), we have discussed algorithms to simulate time evolution of a quantum-computational state under some controllable drive Hamiltonian.
The next step in ctrl-VQE is to optimize the parameters - let's call them $\bar x$ - of the drive Hamiltonian, in order to minimize the expectation value $E$ of an observable, typically a molecular Hamiltonian.
To that end, we want to characterize the gradient vector $\partial E / \partial \bar x$.

### Energy Functions

The energy function we are going to minimize in this note may be wriiten explicitly as:

$$ E = ⟨\psi|\hat O|\psi⟩ $$

The observable $\hat O$ is typically a fermionic molecular Hamiltonian, originally written as a qubit operator in the "computational subspace" but then projected up into a bosonic Hilbert space.

The wavefunction $|\psi⟩$ is the time-evolved state we studied in `timeevolution.md`.
Please see that note for the full treatment; here I will present the bare minimum needed to calculate the gradient:

#### Time Evolution, Abridged

We start with a reference state $|\psi_{\rm ref}⟩$, and evolve it under a Hamiltonian consisting of both static ($H_0$) and drive ( $\hat V(t)$ ) terms, for a total time $T$.
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
I don't know the proper notation to write that explicitly.

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

In the near future (but not yet in this note, as of 3/12/2023), we will also consider the following more experimentally relevant functions:
- Projected: $E = ⟨\psi| \Pi \hat O \Pi |\psi⟩$
- Normalized: $E = ⟨\psi| \Pi \hat O \Pi |\psi⟩ / ⟨\psi|\Pi|\psi⟩$

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

NOTE: Using the projected energy function, the definition of $|\lambda_i⟩$ is very slightly modified, but nothing else that follows is different.

Note that $|\psi_r⟩=|\psi⟩$, $|\psi_0⟩ \ne |\psi_{\rm ref}⟩$, and $E=⟨\lambda_i|\psi_i⟩$ for all $i$.
More importantly:

$$ ⟨\psi|\hat O \left(\partial_x |\psi⟩\right) = \sum_{i=0}^r ⟨\lambda_i| \left( -i \tau_i \partial_x \hat V(t_i) \right) |\psi_i⟩ $$

Thus, the energy derivative becomes:

$$ \partial_x E = - \sum_{i=0}^r \tau_i \left[ ⟨\lambda_i| \left( i \partial_x \hat V(t_i) \right) |\psi_i⟩ + {\rm h.t.} \right] $$

This is the most important equation in this note.
Pretend there is a box around it.

Next, we're going to break this down for two special forms of the drive operator.

#### Control Signals

Let us decompose our drive operators into the following form:

$$ \hat V(t) \equiv \sum_k f_k(\bar x_k, t) \cdot \hat V_k(t) $$

This notation conveys the following three ideas:
1. The drive Hamiltonian is a homogenous sum of individual drive operators $\hat V_k$, each modulated by a _real_ scalar signal $f_k$. This constraint ensures each term is Hermitian. Complex scalar signals can be decomposed into two separate terms, with two distinct drive operators. (We'll do this explicitly in the next section.)
2. All dependence on control parameters is relegated to the scalar signals $f_k$. This simplifies the calculus in a beautiful way, and the results turn out to still be useful in the slightly more general form of $\hat V(t)$ we will consider in the next section.
3. Each scalar signal depends on a _disjoint set_ of control parameters $\{x_k\}$, more conveniently written as a vector $\bar x_k$.

The last constraint isn't hard to relax, but it makes the calculus easier to manage: if $x_{kl}$ is an element of $\bar x_k$, the operator $\partial_{x_{kl}} \hat V(t)$ simplifies to:

$$ \partial_{x_{kl}} \hat V(t) = (\partial_{x_{kl}} f_k)|_t \cdot \hat V_k(t) $$

Substituting this expression into our formula for $\partial_x E$:

$$ \partial_{x_{kl}} E = - \sum_{i=0}^r \tau_i \cdot (\partial_{x_{kl}} f_k)|_{t_i} \cdot \left[ ⟨\lambda_i| \left( i \hat V_k(t_i) \right) |\psi_i⟩ + {\rm h.t.} \right] $$

Now we see the convenience of specifying $f_k$ is real; it factors out of the brackets.
The quantity $(\partial_{x_{kl}} f_k)|_{t_i}$ depends on the precise form of $f_k$ but will generally be easy to compute analytically with elementary calculus.

Thus, the quantity in brackets _does not depend on l_, and only needs to be calculated once for the entire set of parameters $\bar x_k$.
Let us go ahead and explicitly define this quantity as the _gradient signal_ $\vec \phi_k$ which is a time series with individual elements $\phi_{ik}$:

$$ \phi_{ik} \equiv \left[ ⟨\lambda_i| \left( i \hat V_k(t_i) \right) |\psi_i⟩ + {\rm h.t.} \right] $$

For Hilbert spaces tractable in a classical computer, each of the time series $\vec \phi_k$ can be computed in a single sweep, giving access to _all_ the gradients $\partial_{x_{kl}} E$:

$$ \partial_{x_{kl}} E = - \sum_{i=0}^r \tau_i \cdot (\partial_{x_{kl}} f_k)|_{t_i} \cdot \phi_{ik} $$

#### RWA Approximation with Fixed Channels

Now let us consider the drive operator of greatest interest:

$$ \hat V(t) \equiv \sum_k \left[ \Omega_k(\bar x_k, t) \cdot e^{i\nu_k t} \hat a_{q_k} + {\rm h.t.} \right] $$

NOTE: We have traditionally been considering one channel per qubit, but this formalism allows multiple channels per qubit, hence the qubit index $q_k$.

This form differs from that in the previous section in the following ways:
1. The (constant) drive frequency $\nu_k$ is one of our control parameters, but it is not a part of the scalar signal $\Omega_k$.
2. We will allow the scalar signal $\Omega_k$ to be complex.

We will start by insisting that $\Omega_k(\bar x_k, t) \equiv \alpha_k(\bar x_{\alpha k}, t) + i \beta_k(\bar x_{\beta k}, t)$, where the real scalar signals $\alpha_k$ and $\beta_k$ depend on disjoint sets of parameters.
As before, it isn't hard to relax the disjoint constraint but it's probably true in practice and it's definitely easier to write down.
(_Technically_, we could fold $\cos(\nu_k t)$ into $\alpha_k$ and $\sin(\nu_k t)$ into $\beta_k$, but then the parameter sets _aren't_ disjoint, so I'm just going to treat $\nu_k$ separately.)

The important thing to note is that we now have _two_ distinct drive operators for each $k$:

$$ \begin{aligned}
\hat V(t) &= \sum_k \left[ \alpha_k(\bar x_{\alpha k}, t) \cdot \hat V_{\alpha k}(t) + \beta_k(\bar x_{\beta k}, t) \cdot \hat V_{\beta k}(t) \right] \\
\hat V_{\alpha k}(t) &\equiv e^{i\nu_k t} \hat a_{q_k} + {\rm h.t.} \\
\hat V_{\beta k}(t) &\equiv i e^{i\nu_k t} \hat a_{q_k} + {\rm h.t.}
\end{aligned} $$

Thus, we also have two distinct gradient signals for each $k$:

$$ \begin{aligned}
\phi_{i\alpha k} &\equiv \left[ ⟨\lambda_i| \left( i \hat V_{\alpha k}(t_i) \right) |\psi_i⟩ + {\rm h.t.} \right] \\
\phi_{i\beta k} &\equiv \left[ ⟨\lambda_i| \left( i \hat V_{\beta k}(t_i) \right) |\psi_i⟩ + {\rm h.t.} \right]
\end{aligned} $$

For all the parameters appearing in the scalar signals $\Omega_k$, the energy derivatives are hardly distinct from the previous section.

$$ \begin{aligned}
\partial_{x_{\alpha kl}} E = - \sum_{i=0}^r \tau_i \cdot (\partial_{x_{\alpha kl}} \alpha_k)|_{t_i} \cdot \phi_{i\alpha k} \\
\partial_{x_{\beta kl}} E = - \sum_{i=0}^r \tau_i \cdot (\partial_{x_{\beta kl}} \beta_k)|_{t_i} \cdot \phi_{i\beta k}
\end{aligned} $$

Now let us consider the gradient with respect to the channel frequency $\nu_k$:

$$ \begin{aligned}
\partial_{\nu_k} \hat V(t) &= \alpha_k(\bar x_{\alpha k}, t) \cdot \partial_{\nu_k} \hat V_{\alpha k}(t) + \beta_k(\bar x_{\beta k}, t) \cdot \partial_{\nu_k} \hat V_{\beta k}(t) \\
\partial_{\nu_k} \hat V_{\alpha k}(t) &= t \cdot \left[ i e^{i\nu_k t} \hat a_{q_k} + {\rm h.t.} \right] = t \cdot \hat V_{\beta_k}(t) \\
\partial_{\nu_k} \hat V_{\beta k}(t) &= -t \cdot \left[ e^{i\nu_k t} \hat a_{q_k} + {\rm h.t.} \right] = -t \cdot \hat V_{\alpha_k}(t)
\end{aligned} $$

Conveniently, the exact same operators turn out to be relevant - they are now just modulated by an additional factor of t, which can be factored out, and they sort of swap places.
Thus, we may calculate the analytical gradient of the energy with respect to pulse frequencies for essentially zero overhead, using the same gradient signals we have already obtained for complex scalar signals:

$$ \partial_{\nu_k} E = - \sum_{i=0}^r \tau_i \cdot t_i \cdot \left( \alpha_k|_{t_i} \cdot \phi_{i\beta k} - \beta_k|_{t_i} \cdot \phi_{i\alpha k} \right) $$

We now have analytical formulae for every control parameter, _if_ we can only measure the gradient signals $\vec \phi_{\alpha k}$, $\vec \phi_{\beta k}$!
