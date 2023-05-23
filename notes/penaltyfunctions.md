# Penalty Functions
##### 5/22/2023

Control-VQE is much closer to hardware than its gate-model counterpart.
The capabilities of the experimental apparatus itself
    become relevant in the implementation of the algorithm.
For example, the strength of the microwave pulses used to drive qubits
    is limited, either by the capabilities of the waveform generator
    or by the fact that stronger pulses would fry the machine.
Thus, in optimizing these pulses,
    ctrl-VQE must accommodate certain constraints.

For example, let's say we are optimizing a square pulse,
    where our variational parameters $\Omega_i$ are the amplitudes in each time window.
We are trying to minimize an energy $E(\vec\Omega)$.
The note `energygradient.md` provides analytical expressions for the energy,
    meaning we can perform a gradient-based optimization.
However, let's say we require each $|\Omega_i| < \Omega_0$,
    where $\Omega_0$ is the first amplitude which is "too much" for our device.
How do we incorporate this constraint into the optimization?

## Penalty Functions

The solution is to use _penalty functions_.
Let's say we have a function $f(\vec\Omega)$ which is zero if all our constraints are satisfied,
    but some positive value if any are not (ie. if any $\Omega_i > \Omega_0$).
Now, rather than minimizing $E(\vec\Omega)$ directly,
    let us define a new loss function $\mathcal{L}(\vec\Omega)$:

$$ \mathcal{L}(\vec\Omega) = E(\vec\Omega) + \lambda f(\vec\Omega) $$

The _penalty function_ $f(\vec\Omega)$ serves to increase the value of $\mathcal{L}$
    if any $\Omega_i$ violates its constraint.
The idea is that $\mathcal{L}$ should be minimized only when $f(\vec\Omega)=0$.
This isn't necessarily true, though,
    since an out-of-bounds $\Omega_i$ may lower the energy $E(\vec\Omega)$
    by more than the value of $\lambda f(\vec\Omega)$.
Indeed, in analytical optimization, the "Lagrange multiplier" $\lambda$
    acts as another variational parameter whose gradient $\partial_\lambda \mathcal{L}$
    must vanish alonside each $\partial_{\Omega_i} \mathcal{L}$,
    _ensuring_ $f(\vec\Omega)=0$.
This may be a more robust approach (I'm honestly not sure),
    but in our case it isn't essential.
We want to operate above the "minimal evolution time",
    which by definition implies there exists a legal choice of $\vec\Omega$
    for which $E$ is minimized.
Therefore, as long as $\lambda>0$,
    a perfect optimizer must indeed find a $\vec\Omega$ for which $f(\vec\Omega)=0$.
Of course, we don't have perfect optimizers,
    so we'll usually set $\lambda\sim1$, and try to design $f(\vec\Omega)$
    in a clever way that helps the optimizer out.

## Design Principles

**Disclaimer:** I am no expert in optimization.
I assume that experts have thought much harder about what makes a good penalty function than I have.
What follows is simply how I understand the problem, and how I have thought to solve it.


An intuitive (yet unfortunate) choice of penalty function may be
    $f_{\rm BAD}(\vec\Omega) = \sum_i \Theta(|\Omega_i|-\Omega_0)$

Here $\Theta(x)$ represents the Heaviside step function:

$$ \Theta(x) = \left\{ \begin{array}{lcl}
    0 & ~ & x < 0 \\
    1/2 & ~ & x = 0 \\
    1 & ~ & x > 0
\end{array}\right. $$

Its derivative $\partial_x \Theta$ is usually said to be the Dirac-delata function $\delta(x)$:

$$ \delta(x) \equiv \partial_x \Theta \left\{ \begin{array}{lcl}
    0 & ~ & x \ne 0 \\
    \infty & ~ & x = 0
\end{array}\right. $$

The above penalty function $f_{\rm BAD}$ has the required property that $f_{\rm BAD}(x)=0$ when $\Omega_i<\Omega_0$,
    and $f_{\rm BAD}(x) > 0$ when $\Omega_i \ge \Omega_0$.
It is not, however, very clever.

In order to perform a gradient-based optimization,
    we need to find $\partial_{\Omega_i} \mathcal{L} = \partial_{\Omega_i} E + \partial_{\Omega_i} f_{\rm BAD} $.
We know the first term from `energygradient.md`.
We can obtain the second term directly:

$$ \partial_{\Omega_i} f_{\rm BAD} = \delta(|\Omega_i|-\Omega_0) \partial_{\Omega_i} (|\Omega_i|-\Omega_0) $$

This gradient is not very helpful to the optimizer for at least two reasons.
First is the obvious singularity at $|\Omega_i|-\Omega_0$: infinities are generally not good for algorithms.
(Actually, one can't even say if the infinity should be positive or negative, due to the cusp in $|\Omega_i|$:
    absolute values are generally not good for optimizations.)
Second is the fact that the gradient when $|\Omega_i|>\Omega_0$ is _zero_.
In which direction should $\Omega_i$ be changed, to bring $f_{\rm BAD}$ closer to zero?
The gradient doesn't tell us.

The lessons we need to learn from $f_{\rm BAD}$ are these:
1. The penalty function $f(\vec\Omega)$ should be continuous
    (ie. no singularities in $\partial_{\Omega_i} f$).
2. The gradient $\partial_{\Omega_i} f$ should be positive if $\Omega_i > \Omega_0$,
    and negative if $\Omega_i < -\Omega_0$.

In addition to this, we should probably avoid absolute values $|\Omega_i|$ wherever possible,
    and we should design the magnitude of the gradient $\partial_{\Omega_i} f$
    to be _larger_ the further away $\Omega_i$ is from a legal value.

The simplest function satisfying these rules would be something like
    $f_{\rm MEH}(\vec\Omega) = \sum_i {\rm ramp}(\Omega_i^2-\Omega_0^2)$, where the ramp function is:

$$ {\rm ramp}(x) = \left\{ \begin{array}{lcl}
    0 & ~ & x < 0 \\
    x & ~ & x \ge 0
\end{array}\right. $$

The gradient of $f_{\rm MEH}$ is:

$$ \partial_{\Omega_i} f_{\rm MEH} = \left\{ \begin{array}{lcl}
    0 & ~ & |\Omega_i| < \Omega_0 \\
    2\Omega_i & ~ & |\Omega_i| \ge \Omega_0
\end{array}\right. $$

I call this function "meh" because, while we have removed the singularity from the gradient,
    there is still a glaring discontinuity,
    which will probably muck the optimizer up near saturation, where $|\Omega_i|\sim\Omega_0$.
You can imagine the energy gradient kicking a legal $\Omega_i$ over the top,
    while in the next iteration the penalty gradient kicks it back.
They'll just go back and forth forever.
We should avoid kicks; we need the penalty gradient to _slide_ it back.
That means the function $f$ ought to be _differentiable_, ie. its derivative is continuous.
For quasi-Newton methods like BFGS,
    the optimization routine calculates something like a second-derivative,
    so it is probably best to ensure that the gradient is itself differentiable.
Indeed, I conjecture one might as well pick $f$ to be _smooth_,
    meaning _every_ derivative is continuous.

This sounds daunting.
We absolutely need a piecewise-function, so that $f(\vec\Omega)=0$ for any legal $\vec\Omega$.
How can we make a function which is smooth at this interface?

Well, here's one way:
    let $f_{\rm OK}(\vec\Omega) = \sum_i {\rm smoothramp}(\Omega_i^2-\Omega_0^2)$, where the smoothed-out ramp function is:

$$ {\rm smoothramp}(x) = \left\{ \begin{array}{lcl}
    0 & ~ & x < 0 \\
    x \cdot e^{-1/x} & ~ & x \ge 0
\end{array}\right. $$

I've just added in a factor $e^{-1/x}$, which serves to smooth out the transition at $x\sim0$.

Defining $x_i \equiv \Omega_i^2 - \Omega_0^2$, the gradient of $f_{\rm OK}$ is:

$$ \partial_{\Omega_i} f_{\rm OK} = \left\{ \begin{array}{lcl}
    0 & ~ & |\Omega_i| < \Omega_0 \\
    2\Omega_i \cdot e^{-1/x_i} \cdot(1 + 1/x_i) & ~ & |\Omega_i| \ge \Omega_0
\end{array}\right. $$

This function has everything we need, and it should work perfectly well.
(Note that the apparent singularity from $1/x_i$ is neutralized by the faster falling $e^{-1/x_i}$.)

I call this function merely "ok" simply because it looks ugly to me.
The next section introduces an equally suitable function
    which also happens to satisfy my personal aesthetic,
    and is of course the one I actually implemented in code.
It will be presented in a quite different way, but it's actually very similar.
The differences reduce to just two things:
1. I avoid the very ugly $\Omega_i^2-\Omega_0^2$ terms
    by explicitly placing a different "ramp" on each bound (ie. $\pm\Omega_0$),
2. I replace the linear ramp with an exponential function $e^x$,
    which looks better next to the smoothing factor.

## My Choice: Smoothed Exponential

My choice of penalty function is:

$$ f(\vec\Omega) \equiv \sum_i g(\Omega_i/\Omega_0) $$

The element-wise filter function $g(x)$ is:

$$ g(x) \equiv \left\{ \begin{array}{lcl}
    h(-x-1) & ~ & x < -1 \\
    0 & ~ & -1 \le x \le 1 \\
    h(x-1) & ~ & x > 1
\end{array} \right. $$

The smoothed exponential $h(y)$ is defined as:

$$ h(y) \equiv \exp(y-1/y) $$

Its derivative is:

$$ \partial_y h = h(y) \cdot (1+1/y^2) $$

Think of $x$ as a normalized amplitude, and $y$ as a normalized displacement from the bounds.
In that spirit, define $x_i \equiv \Omega_i/\Omega_0$, and $y_i \equiv |x_i| - 1$.

Then the derivative of $g(x_i)$ is:

$$ \partial_{x_i} g = \left\{ \begin{array}{lcl}
    -h(y_i) \cdot (1+1/y_i^2) & ~ & x_i < -1 \\
    0 & ~ & -1 \le x_i \le 1 \\
    +h(y_i) \cdot (1+1/y_i^2) & ~ & x_i > 1
\end{array} \right. $$

The gradient of $f(\vec\Omega)$ is:

$$ \partial_{\Omega_i} f = \partial_{x_i} g / \Omega_0 $$