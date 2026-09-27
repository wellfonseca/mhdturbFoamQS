# mhdturbFoamQS

**A quasi-static MHD solver for incompressible turbulent flows with non-uniform
applied magnetic fields**, for OpenFOAM 6.

`mhdturbFoamQS` solves the incompressible Navier–Stokes equations coupled to a
low-magnetic-Reynolds-number MHD model through the **electric-potential
(quasi-static) formulation**, with the applied magnetic field $\mathbf{B}_0$
prescribed as a general, spatially varying field. It supports the standard
OpenFOAM turbulence closures unchanged — no magnetic term is added to the
closures.

It exists because the **full magnetic-induction formulation** (OpenFOAM's
`mhdFoam`, and the turbulent variant `mhdturbFoam`) **fails when the applied
field is not uniform**. That failure is not a matter of tuning: it is structural,
and it is documented below.

---

## 1. Governing equations

Incompressible, low-magnetic-Reynolds-number MHD:

$$
\nabla\cdot\mathbf{u} = 0
$$

$$
\frac{\partial \mathbf{u}}{\partial t}
+ \nabla\cdot(\mathbf{u}\mathbf{u})
= -\nabla p + \nabla\cdot\left(\nu_\text{eff}\,\nabla \mathbf{u}\right)
+ \frac{1}{\rho}\,\mathbf{j}\times\mathbf{B}_0
$$

$$
\mathbf{j} = \sigma\left(-\nabla\phi + \mathbf{u}\times\mathbf{B}_0\right),
\qquad
\nabla\cdot\mathbf{j} = 0
\;\;\Longrightarrow\;\;
\nabla^2\phi = \nabla\cdot\left(\mathbf{u}\times\mathbf{B}_0\right)
$$

where $\phi$ is the electric potential, $\mathbf{j}$ the current density,
$\mathbf{B}_0(\mathbf{x})$ the **prescribed** magnetic field (not evolved), and
$\nu_\text{eff}$ the effective viscosity of whichever turbulence closure is
selected. Electrically insulating walls are imposed as
$\mathbf{j}\cdot\mathbf{n} = 0$.

The numerical implementation computes the cell-centred current density from the
face fluxes by the standard surface-integral identity,

$$
\mathbf{j}_P \simeq
\frac{1}{V_P}\sum_f \left(\mathbf{j}\cdot\mathbf{S}_f\right)\mathbf{C}_f
- \left[\frac{1}{V_P}\sum_f \left(\mathbf{j}\cdot\mathbf{S}_f\right)\right]\mathbf{C}_P ,
$$

so that the Lorentz force follows directly from $\mathbf{j}$ and $\mathbf{B}_0$.

## 2. Why not the full induction formulation

`mhdFoam` and its turbulent derivatives evaluate the Lorentz force in
**Maxwell-stress form** from the *total* field $\mathbf{B} = \mathbf{B}_0 + \mathbf{b}$:

$$
\mathbf{f} = -\nabla\cdot\left(\frac{\mathbf{B}\mathbf{B}}{\mu_0}\right)
+ \nabla\left(\frac{|\mathbf{B}|^2}{2\mu_0}\right).
$$

For a **uniform** $\mathbf{B}_0$ the discrete operators annihilate the imposed
part exactly, and only the (small, $O(R_m)$) induced contribution survives — so
the formulation works. For a **non-uniform** $\mathbf{B}_0$ the two terms are
each of order $B_0^2/(\mu_0\rho)$ — about $4\times10^4\ \mathrm{m^2/s^2}$
for $B_0 = 10$ T — and their discrete cancellation is only as good as the
truncation error. The residual spurious force scales as

$$
f_\text{spurious} \sim \left(\frac{\Delta x}{L}\right)^2 \frac{B_0^2}{\mu_0 L \rho},
$$

which for the meshes and fields of a real magnet arrangement is of the order of
$10^3\ \mathrm{m^2/s^2}$ — **several times larger than the physical Lorentz
force** ($\approx 200\ \mathrm{m^2/s^2}$ in the same conditions). The result is a
non-physical body force that destabilises the solution.

This solver avoids the cancellation entirely: the current is obtained from a
Poisson problem for the potential, so no large term is ever subtracted from
another.

## 3. Relation to existing codes

`mhdturbFoamQS` is a combination of two existing, independently published pieces
plus a small number of additions. See `NOTICE` for the full attribution.

| Component | Origin |
|---|---|
| Turbulent MHD solver skeleton (`mhdturbFoam`: PISO + `divDevReff` + closures) | FOSSEE case-study project, R. Radhakrishnan (2019), based on OpenFOAM's `mhdFoam` |
| Electric-potential (quasi-static) coupling | `epotFoam`, A. Tassone (2016), Chalmers OpenFOAM course |

Additions made here:

1. $\mathbf{B}_0$ promoted from a constant `dimensionedVector` to a
   **`volVectorField`**, so that arbitrary non-uniform applied fields can be
   prescribed (e.g. the field of a set of permanent magnets).
2. Insulating-wall condition $\mathbf{j}\cdot\mathbf{n}=0$ implemented by
   **excluding the boundary flux of $\mathbf{u}\times\mathbf{B}_0$** from both
   sides of the potential equation, consistently with the `zeroGradient`
   condition on $\phi$.
3. Support for `fvOptions` in the momentum equation (periodic forcing,
   `meanVelocityForce`, constant pressure-gradient sources).
4. Reference cell for the singular magnetic-pressure/potential Poisson problem
   in fully periodic or closed domains.

## 4. Building

Requires OpenFOAM 6 (org) and a C++11 compiler.

```bash
mkdir -p $WM_PROJECT_USER_DIR/applications/solvers/magnetohydrodynamics
cp -r src $WM_PROJECT_USER_DIR/applications/solvers/magnetohydrodynamics/mhdturbFoamQS
cd $WM_PROJECT_USER_DIR/applications/solvers/magnetohydrodynamics/mhdturbFoamQS
wmake
```

The executable is installed in `$FOAM_USER_APPBIN`.

## 5. Running

The case must provide, in `0/`:

| Field | Description |
|---|---|
| `U`, `p` | velocity and kinematic pressure |
| `B0` | **prescribed** applied field (`volVectorField`, `fixedValue` on magnet patches) |
| `PotE` | electric potential (`zeroGradient` on insulating walls) |
| turbulence fields | `nuTilda`(SA), `k`/`epsilon`, `k`/`omega` — as usual |

and in `system/fvSolution` a `PotE` entry (solver) plus a `PotE` sub-dictionary
with `PotERefCell`/`PotERefValue`.

```bash
decomposePar
mpirun -np N mhdturbFoamQS -parallel
reconstructPar -latestTime
```

## 6. Verification

The reference for verification is the classical analytical solution of the
fully developed MHD flow in a circular pipe with insulating walls
(Shercliff 1953; Gold 1962; Uhlenbusch & Fischer 1961), which gives
$f\,Re$ as a function of the Hartmann number independently of $Re$.

`cases/channelHartmann` contains a small plane-channel Hartmann case with the
exact solution, used to verify the implementation in seconds. Results are
reported in the accompanying paper and in `docs/formulation.md`.

## 7. Status and limitations

- Valid for $R_m = \mu_0\sigma U L \ll 1$; the induced field is neglected by
  construction. In the applications reported so far $R_m \approx 3.7\times10^{-4}$.
- The applied field $\mathbf{B}_0$ must be divergence-free and curl-free
  ($\nabla\times\mathbf{B}_0 = 0$) for the initial state to be current-free;
  the solver does not enforce this.
- Turbulence closures are used unmodified. This is a deliberate modelling
  choice, not an oversight: it is precisely the question the accompanying study
  addresses.

## 8. Citation

See `CITATION.cff`. If you use this solver, please cite the accompanying paper
**and** the two upstream works listed in `NOTICE`.

The repository is at <https://github.com/wellfonseca/mhdturbFoamQS>.
A permanent archive with a DOI is deposited at Zenodo (see `CITATION.cff`).

## 9. License

GNU General Public License v3.0 — see `LICENSE`. This is mandatory: the code
derives from GPL-licensed OpenFOAM, FOSSEE and Chalmers sources.
