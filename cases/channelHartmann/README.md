# `channelHartmann` — plane-channel verification against the exact solution

A plane channel with the **same geometry as the OpenFOAM 6 tutorial**
`tutorials/electromagnetics/mhdFoam/hartmann`: length 20 m, height 2 m
(half-width $h = 1$ m) and thickness 0.1 m. The layout follows the same OpenFOAM
conventions as that tutorial — `constant/fvOptions`, `#includeFunc` function
objects driven by dictionaries in `system/`, and an `Allrun` that ends by
sampling the solution profile.

The one deliberate difference is that the streamwise direction is **cyclic**
instead of inlet/outlet. With an inlet the flow is still developing, and the
closed-form Hartmann solution — which describes *fully developed* flow — would
not apply; with a cyclic patch the flow is fully developed by construction and
the exact solution holds. The mesh follows the OpenFOAM convention for a fully
developed periodic channel (cf. `planarPoiseuille`): a single cell in the
streamwise direction, since nothing varies along $x$, and 100 cells across the
channel so that the Hartmann layer ($h/Ha = 0.1$ m, i.e. 5 cells at $Ha = 10$)
is resolved.

The flow is driven by a constant pressure gradient imposed with
`vectorSemiImplicitSource` in `constant/fvOptions`
($G = 1.309928\times10^{-5}$ m/s², chosen so that $Re = 100$), so the bulk
velocity is an *output*, not an input.

## Exact solution

With $Ha = hB_0\sqrt{\sigma/\rho\nu}$, $D_h = 4h$ and $Re = \bar U D_h/\nu$:

$$
\frac{u(y)}{\bar U} = \frac{1-\cosh(Ha\,y/h)/\cosh(Ha)}{1-\tanh(Ha)/Ha},
\qquad
f\,Re = \frac{32\,Ha^{2}}{1-\tanh(Ha)/Ha},
\qquad
f = \frac{2D_h G}{\bar U^2}.
$$

Limits: $f\,Re \to 96$ as $Ha\to0$ (plane Poiseuille) and $f\,Re \to 32Ha^{2}$
for large $Ha$.

## Expected values

| $Ha$ | $B_0$ (T) | $f\,Re$ exact | $f\,Re$ computed | deviation |
|---|---|---|---|---|
| 0.0 | 0.0000 | 96.0000 | 96.0066 | $+0.007$ % |
| 0.5 | 0.0384 | 105.5887 | 105.5804 | $-0.008$ % |
| 1.0 | 0.0768 | 134.2249 | 134.2008 | $-0.018$ % |
| 5.0 | 0.3839 | 999.9773 | 999.6654 | $-0.031$ % |
| 10.0 | 0.7678 | 3555.5556 | 3553.5961 | $-0.055$ % |

The computed values are the ones produced by the sampled profile described
below; the largest deviation is 0.055 %, at the strongest field. The five runs
take about 20 s in total.

## Running

Requires OpenFOAM 6 and the solver installed (`wmake` in `../../src`).

```bash
bash Allrun      # blockMesh + mhdturbFoamQS + sample the profile
```

`Allrun` fails loudly if a step fails, unlike the stock OpenFOAM script, which
carries on and leaves the user with a "successful" run and no results.

`0/B0` ships with $B_0 = 0$ (plain plane Poiseuille). To run a Hartmann case, set
both the `internalField` and the wall `value` of `0/B0` to $(0, B_0, 0)$ with the
$B_0$ from the table above, then rerun. `Allclean` removes the results.

`Allrun` writes three things, all of them through standard OpenFOAM function
objects declared with `#includeFunc` in `system/controlDict`:

| Output | Dictionary | What it is |
|---|---|---|
| `postProcessing/residuals/` | `system/residuals` | residual history of `p`, `U` and `PotE` |
| `postProcessing/volFieldValue/` | `system/volFieldValue` | volume averages of `U` and `p`, every 20 steps |
| `postProcessing/sample/` | `system/sample` | the velocity profile across the channel at each written time |

The profile is a line of 100 points placed **exactly on the 100 cell centres**,
so the sampled data reproduce the field without interpolation error. The file
`postProcessing/sample/<time>/centreProfile_U.xy` has four columns —
distance, $U_x$, $U_y$, $U_z$ — and the average of $U_x$ over the 100 points is
the bulk velocity $\bar U$ used above.

## Notes

- `0/B0` is a `volVectorField` with `internalField` and wall values equal to
  $(0, B_0, 0)$; both have to be edited for each $Ha$.
- `0/PotE` uses `zeroGradient` on the walls, together with the exclusion of the
  boundary flux of $\mathbf{u}\times\mathbf{B}_0$ on the non-coupled patches in
  the solver — the two halves of the insulating condition. Changing one without
  the other imposes a conducting wall and the verification will fail. The flux
  is kept on cyclic patches: there it is internal, and dropping it would make
  the right-hand side inconsistent with the Laplacian.
- The viscous transient of this geometry is long: the slowest channel mode
  decays with a time constant $4h^2/(\pi^2\nu) \approx 970$ s, which is why the
  case runs to $t = 8000$ s. The steady state itself does not depend on the time
  step; a larger $\Delta t$ only changes how the transient is resolved.
- The time step is not limited by the magnetic damping: the stiff part
  $-\sigma|\mathbf{B}_0|^2\mathbf{U}$ of the Lorentz force is discretised
  implicitly in the momentum equation.
- To inspect the field in ParaView, run `foamToVTK`. It names the files after
  the **time index** (the step count), not the physical time.
