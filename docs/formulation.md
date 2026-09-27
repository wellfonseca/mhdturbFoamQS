# Formulation, discretisation and verification

## 1. Physical model

Incompressible, electrically conducting Newtonian fluid in the
low-magnetic-Reynolds-number regime. Lengths are in metres, $\mathbf{B}_0$ in
tesla, $\phi$ in volts.

$$
R_m = \mu_0\sigma U L \ll 1
$$

so the induced magnetic field is neglected and the applied field
$\mathbf{B}_0(\mathbf{x})$ enters as a prescribed, divergence-free, curl-free
field.

| Equation | Form |
|---|---|
| Mass | $\nabla\cdot\mathbf{u} = 0$ |
| Momentum | $\dfrac{\partial\mathbf{u}}{\partial t} + \nabla\cdot(\mathbf{u}\mathbf{u}) = -\nabla p + \nabla\cdot(\nu_\text{eff}\nabla\mathbf{u}) + \dfrac{1}{\rho}\mathbf{j}\times\mathbf{B}_0$ |
| Current | $\mathbf{j} = \sigma\left(-\nabla\phi + \mathbf{u}\times\mathbf{B}_0\right)$ |
| Charge | $\nabla\cdot\mathbf{j} = 0 \;\Rightarrow\; \nabla^2\phi = \nabla\cdot(\mathbf{u}\times\mathbf{B}_0)$ |

### Dimensionless groups

With $D$ a characteristic length, $U$ the bulk velocity and $B_0$ the applied
field magnitude:

$$
Re = \frac{UD}{\nu}, \qquad
Ha = B_0 D\sqrt{\frac{\sigma}{\rho\nu}}, \qquad
N = \frac{\sigma B_0^2 D}{\rho U}, \qquad
R_m = \mu_0\sigma U D .
$$

$N$ (the interaction parameter) measures the Lorentz force against inertia;
$Ha$ measures it against viscous diffusion.

### Exact solution used as reference

For a **circular pipe** of radius $a$ with insulating walls and a uniform
transverse field, the fully developed solution is that of Shercliff (1953),
Uhlenbusch & Fischer (1961) and Gold (1962). In the present notation it solves

$$
\nabla^2 w - Ha_a^2\left(w + \frac{\partial\varphi}{\partial x}\right) = -K,
\qquad
\nabla^2\varphi = -\frac{\partial w}{\partial x},
$$

with $w=0$ and $\partial\varphi/\partial r = 0$ on $r=1$, $Ha_a = aB_0\sqrt{\sigma/\rho\nu} = Ha/2$
and $K = (-\mathrm{d}p/\mathrm{d}x)a^2/(\rho\nu \bar U)$. Because the problem is
linear in $K$, its solution fixes the friction factor directly,

$$
f\,Re = \frac{8}{\langle w\rangle},
$$

with $\langle w\rangle$ the area average of the solution obtained with $K=1$.
The limit $Ha\to0$ gives $f\,Re = 64$ exactly.

For a **plane channel** of half-width $h$ the same analysis gives the closed form

$$
\frac{u(y)}{\bar U} = \frac{1 - \cosh(Ha\, y/h)/\cosh(Ha)}{1 - \tanh(Ha)/Ha},
\qquad
f\,Re = \frac{32\,Ha^{2}}{1-\tanh(Ha)/Ha},
$$

with $Re = \bar U D_h/\nu$ and $D_h = 4h$; the limit $Ha\to0$ gives $f\,Re = 96$.

## 2. Why the induction formulation fails with a non-uniform field

The full-induction solvers evaluate the Lorentz force from the total field
$\mathbf{B}=\mathbf{B}_0+\mathbf{b}$ as

$$
\mathbf{f} = -\nabla\cdot\left(\frac{\mathbf{B}\mathbf{B}}{\mu_0}\right)
+ \nabla\left(\frac{|\mathbf{B}|^2}{2\mu_0}\right).
$$

Each term is of order $B_0^2/(\mu_0 L)$. For a **constant** $\mathbf{B}_0$ the
discrete divergence and gradient annihilate it exactly, leaving only the
$O(R_m)$ induced contribution — which is why those solvers work for uniform
fields. For a **non-uniform** $\mathbf{B}_0$ the cancellation is only as good as
the truncation error of the two different discrete operators, and the residual

$$
f_\text{spurious} \sim \left(\frac{\Delta x}{L}\right)^2\frac{B_0^2}{\mu_0 L\rho}
$$

is a genuine, non-physical body force.

Measured example (the configuration studied in the accompanying work): pipe of
$D_h = 0.152$ m, $B_0=10$ T, $\Delta x \approx 2.4$ mm, $\mu_0 = 4\pi\times10^{-7}$,
$\rho = 10^3$ kg/m³, $\sigma = 70.9$ S/m, $U = 27.5$ m/s:

| Quantity | Value |
|---|---|
| $B_0^2/(2\mu_0\rho)$ | $3.98\times10^4$ m²/s² |
| physical Lorentz acceleration $\sigma B_0^2U/\rho$ | $1.9\times10^2$ m²/s² |
| spurious acceleration from the discrete cancellation | $\sim10^3$ m²/s² |
| outcome | floating-point exception within two time steps |

The spurious term is thus several times the physical force. `mhdturbFoamQS`
removes it by construction: the current is obtained from the potential, and no
order-$B_0^2$ quantity is ever subtracted from another.

## 3. Discretisation

Finite-volume, collocated, second order. The momentum equation is solved with
the standard PISO pressure–velocity coupling of the OpenFOAM framework
(Issa 1986; Weller et al. 1998):

1. momentum predictor with the Lorentz force evaluated from the previous
   potential;
2. pressure correctors, with `pFinal` used on the last non-orthogonal iteration;
3. explicit `fvOptions` correction.

The potential equation is then assembled and solved,

$$
\sum_f \left(\nabla\phi\cdot\mathbf{S}_f\right)
= \sum_f \left(\mathbf{u}\times\mathbf{B}_0\right)\cdot\mathbf{S}_f ,
$$

with a reference cell for the arbitrary constant.

### Insulating walls: $\mathbf{j}\cdot\mathbf{n}=0$

Integrating $\nabla\cdot\mathbf{j}=0$ over a boundary cell, the wall face
contributes $j_n = -\partial\phi/\partial n + (\mathbf{u}\times\mathbf{B}_0)\cdot\mathbf{n}$.
Requiring $j_n=0$ gives the Neumann condition
$\partial\phi/\partial n = (\mathbf{u}\times\mathbf{B}_0)\cdot\mathbf{n}$.
In the discrete system this is imposed **consistently on both sides**: $\phi$ is
given `zeroGradient` (so the wall contributes no Laplacian flux) *and* the
boundary flux of $\mathbf{u}\times\mathbf{B}_0$ is excluded from the right-hand
side. Omitting the second half silently imposes a conducting wall instead.

### Cell-centred current

$$
\mathbf{j}_P \simeq
\frac{1}{V_P}\sum_f \left(\mathbf{j}\cdot\mathbf{S}_f\right)\mathbf{C}_f
- \left[\frac{1}{V_P}\sum_f \left(\mathbf{j}\cdot\mathbf{S}_f\right)\right]\mathbf{C}_P .
$$

The Lorentz force is then $\sigma\,\mathbf{j}_P\times\mathbf{B}_{0,P}$.

## 4. Verification procedure

1. **Plane-channel Hartmann case** (`cases/channelHartmann`), closed-form
   solution, runs in seconds. Checks the potential solve, the insulating
   condition, the current reconstruction and the Lorentz force.
2. **Circular pipe, periodic, versus Gold/Shercliff**, for $Ha = 0, 1, 5, 10$,
   which for $Ha=0$ must also reproduce the exact $f\,Re = 64$.
3. **Grid and time-step refinement** of both.
4. **Cross-check against the full-induction solver** in its valid regime
   (uniform field), where both formulations must agree.

Items 2–4 require a running OpenFOAM installation and are performed in the
accompanying study.
