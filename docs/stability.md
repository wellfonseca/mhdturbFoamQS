# Stability of the quasi-static MHD scheme

The time step of `mhdturbFoamQS` is not limited by the magnetic field, and this
note shows why. The result is a von Neumann analysis of the discretised
Lorentz force, together with the numerical confirmation described at the end.

## Governing equations

With a prescribed, constant, divergence-free and curl-free field
$\mathbf{B}_0$, the low-$Rm$ system is

$$
\frac{\partial \mathbf{U}}{\partial t} + (\mathbf{U}\cdot\nabla)\mathbf{U}
  = -\nabla\frac{p}{\rho} + \nu\nabla^2\mathbf{U}
  + \frac{1}{\rho}\mathbf{j}\times\mathbf{B}_0,
\qquad
\nabla\cdot\mathbf{U} = 0,
$$

$$
\mathbf{j} = \sigma\left(-\nabla\phi + \mathbf{U}\times\mathbf{B}_0\right),
\qquad
\nabla\cdot\mathbf{j} = 0
\;\Longrightarrow\;
\nabla\cdot(\sigma\nabla\phi) = \nabla\cdot(\sigma\,\mathbf{U}\times\mathbf{B}_0).
$$

Using $(\mathbf{U}\times\mathbf{B}_0)\times\mathbf{B}_0
= \mathbf{B}_0(\mathbf{U}\cdot\mathbf{B}_0) - |\mathbf{B}_0|^2\mathbf{U}$,
the Lorentz force splits into three pieces:

$$
\frac{1}{\rho}\mathbf{j}\times\mathbf{B}_0
  = \underbrace{-\frac{\sigma}{\rho}\nabla\phi\times\mathbf{B}_0
    + \frac{\sigma}{\rho}\mathbf{B}_0(\mathbf{U}\cdot\mathbf{B}_0)}_{\text{explicit}}
  \underbrace{-\alpha\mathbf{U}}_{\text{implicit}},
\qquad
\alpha \equiv \frac{\sigma|\mathbf{B}_0|^2}{\rho}.
$$

The last term is the stiff one: it is a pure sink with the relaxation time
$1/\alpha$, and it is the only piece that grows with $|\mathbf{B}_0|^2$ in a way
that can violate an explicit time-step limit. In the solver it is discretised
implicitly as `fvm::Sp(sigma*magSqr(B0)/rho, U)`, while the other two pieces are
formed explicitly from the previous step's fields.

## Von Neumann analysis

Write a perturbation as $\hat{\mathbf{U}}^n e^{\mathrm{i}\mathbf{k}\cdot\mathbf{x}}$.
The potential equation is elliptic, so its Fourier symbol is algebraic: from
$-k^2\hat{\phi} = \mathrm{i}\mathbf{k}\cdot(\hat{\mathbf{U}}\times\mathbf{B}_0)$,

$$
\hat{\phi} = -\frac{\mathrm{i}}{k^2}\,\mathbf{k}\cdot(\hat{\mathbf{U}}\times\mathbf{B}_0).
$$

Take $\mathbf{B}_0 = B_0\mathbf{e}_z$ and, without loss of generality,
$\mathbf{k} = k\mathbf{e}_x$. Then
$\mathbf{k}\cdot(\hat{\mathbf{U}}\times\mathbf{B}_0) = kB_0\hat{U}_y$,
$\hat{\phi} = -\mathrm{i}B_0\hat{U}_y/k$, and the explicit pair contributes

$$
-\mathrm{i}k\hat{\phi}\,\mathbf{B}_0\;\text{-part} + B_0(\hat{\mathbf{U}}\cdot\mathbf{B}_0)
  = B_0^2\left(\hat{U}_y\mathbf{e}_y + \hat{U}_z\mathbf{e}_z\right),
$$

so the complete explicit Lorentz symbol is
$\alpha(0,\hat{U}_y,\hat{U}_z)$. The discretised momentum equation is therefore

$$
\left(1+\alpha\Delta t\right)\hat{\mathbf{U}}^{n+1}
  = \hat{\mathbf{U}}^{n} + \alpha\Delta t\left(0,\hat{U}_y^n,\hat{U}_z^n\right),
$$

giving the amplification factors

| component | $G$ | behaviour |
|---|---|---|
| $\hat{U}_x$ (across $\mathbf{B}_0$ and $\mathbf{k}$) | $1/(1+\alpha\Delta t)$ | decays, $\lvert G\rvert < 1$ for every $\Delta t$ |
| $\hat{U}_y$ (across $\mathbf{B}_0$, along $\mathbf{k}$) | $1$ exactly | neutral |
| $\hat{U}_z$ (along $\mathbf{B}_0$) | $1$ exactly | neutral |

**The magnetic terms impose no time-step limit on the discretisation used here:
$\lvert G\rvert\le 1$ for every $\Delta t$.**

Two remarks make the result more than an algebraic accident:

1. *Why the split is safe.* The two explicit pieces and the implicit sink act on
   the same coefficient $\alpha$. Wherever both are present, as in $\hat{U}_y$
   and $\hat{U}_z$, they cancel exactly and leave $G=1$; where only the sink is
   present, as in $\hat{U}_x$, it is implicit and therefore contractive. Splitting
   an operator into an implicit and an explicit part with the *same* eigenvalue
   is what makes the scheme unconditionally stable rather than merely less
   restrictive.
2. *What the explicit treatment would give.* Had the sink also been explicit,
   the $\hat{U}_x$ factor would be $G=1-\alpha\Delta t$, stable only while
   $\alpha\Delta t\le 2$, i.e.
   $\Delta t \le 2\rho/(\sigma|\mathbf{B}_0|^2)$. This is the classical explicit
   damping limit, and it is the reason the original formulation failed. The
   first version of the verification case ran at $Ha=10$ with
   $\Delta t = 0.2$ s against a limit of $0.12$ s, and the duct campaign at
   $10$ T likewise needs $\Delta t<0.12$ s while running at $0.2$ s; both
   diverged. The failure is not a question of tuning: at fixed $\Delta t$ the
   limit is crossed as soon as the field is strong enough.

## The other terms

For completeness, the remaining operators carry their own standard conditions:

| term | discretisation | von Neumann condition |
|---|---|---|
| convection, first-order upwind | explicit | $\mathrm{Co}=|\mathbf{U}|\Delta t/\Delta x \le 1$ |
| convection, central (`Gauss linear`) | explicit | neutrally stable for the semi-discrete operator; in practice bounded by the pressure correction and by physical diffusion, and the cases run at $\mathrm{Co}\lesssim 1$ |
| molecular and turbulent diffusion | implicit | none |
| pressure (PISO) | implicit, 3 correctors | none |
| potential $\phi$ | implicit, once per step from $\mathbf{U}^{n+1}$ | none; the one-step lag enters the explicit Lorentz symbol above |

There is also **no Alfvén-wave condition**, and this is a property of the
formulation rather than of the discretisation: with $\mathbf{B}_0$ prescribed
there are no magnetic waves to resolve, so neither the Alfvén speed
$B_0/\sqrt{\mu_0\rho}$ nor the magnetic diffusion time appears in the stability
limit. Only the flow Courant number remains, which is what the case settings
show.

## Numerical confirmation

The `channelHartmann` case at $Ha=10$ has $\alpha = 4.18\times10^{-2}$ s$^{-1}$,
so the explicit limit is $\Delta t<47.9$ s. Two builds were run on the same case
with $\Delta t = 500$ s, ten times above that limit:

| step | implicit sink (this solver) | explicit sink |
|---|---|---|
| $t=500$ | $1$ | $1$ |
| $t=1000$ | $2.4\times10^{-2}$ | $0.41$ |
| $t=1500$ | $1.0\times10^{-3}$ | floating-point exception |
| $t=5000$ | $6.1\times10^{-12}$, clean exit | — |

(the figures are the initial residual of $U_x$). The implicit build decays
eleven orders of magnitude and finishes; the explicit build stalls and then
aborts with a floating-point exception, exactly as the analysis predicts.
