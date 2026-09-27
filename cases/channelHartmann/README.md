# `channelHartmann` — verification against the exact Hartmann solution

A plane channel of half-width $h = 0.05$ m, periodic in the streamwise direction
(one cell), with a uniform magnetic field normal to the walls and **insulating**
walls. It is the cheapest case that exercises the whole coupling: the potential
Poisson solve, the insulating-wall condition $\mathbf{j}\cdot\mathbf{n}=0$, the
cell-centred current reconstruction and the Lorentz force.

The flow is driven by a constant pressure gradient imposed with
`vectorSemiImplicitSource` ($G = 0.0105$ m/s²), so the bulk velocity is an
*output*, not an input.

## Exact solution

With $Ha = hB_0\sqrt{\sigma/\rho\nu}$, $D_h = 4h$ and $Re = \bar U D_h/\nu$:

$$
\frac{u(y)}{\bar U} = \frac{1-\cosh(Ha\,y/h)/\cosh(Ha)}{1-\tanh(Ha)/Ha},
\qquad
f\,Re = \frac{32\,Ha\tanh(Ha)}{1-\tanh(Ha)/Ha},
\qquad
f = \frac{2D_h G}{\bar U^2}.
$$

Limits: $f\,Re \to 96$ as $Ha\to0$ (plane Poiseuille) and $f\,Re \to 32Ha$ for
large $Ha$.

## Expected values

| $Ha$ | $B_0$ (T) | $f\,Re$ exact |
|---|---|---|
| 0.0 | 0.0000 | 96.0000 |
| 0.5 | 0.7678 | 97.5887 |
| 1.0 | 1.5355 | 102.2249 |
| 5.0 | 7.6776 | 199.9773 |
| 10.0 | 15.3551 | 355.5556 |

## Running

Requires OpenFOAM 6 and the solver installed (`wmake` in `../../src`).

```bash
bash Allrun      # blockMesh + mhdturbFoamQS
```

`0/B0` ships with $B_0 = 0$ (plain plane Poiseuille). To run a Hartmann case, set
both the `internalField` and the wall `value` of `0/B0` to $(0, B_0, 0)$ with the
$B_0$ from the table above, then rerun. The bulk velocity is written as a volume
average of `U` in `postProcessing/volFieldValue1/`; take its last value and form
$f\,Re$ with the relations given above. `Allclean` removes the results.

The case is deliberately tiny (100 cells, one cell in the streamwise direction):
one $Ha$ runs in seconds, which is what makes it usable as a routine regression
check.

## Notes

- `0/B0` is a `volVectorField` with `internalField` and wall values equal to
  $(0, B_0, 0)$; both have to be edited for each $Ha$.
- `0/PotE` uses `zeroGradient` on all walls, together with the exclusion of the
  boundary flux of $\mathbf{u}\times\mathbf{B}_0$ in the solver — the two halves
  of the insulating condition. Changing one without the other imposes a
  conducting wall and the verification will fail.
- The channel is resolved with 100 uniform cells across the half-width, which
  resolves the Hartmann layer ($h/Ha = 5$ mm at $Ha = 10$) with several cells.
