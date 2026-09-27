#!/usr/bin/env python3
"""
make_case.py -- write the plane-channel Hartmann verification case.

Geometry: x in [0, L] (one cell, cyclic), y in [-h, h] (walls), z (one cell,
empty). Uniform magnetic field B0 = (0, B0, 0), normal to the walls, with
electrically insulating walls. The flow is driven by a fixed pressure gradient
imposed with vectorSemiImplicitSource.

Exact solution (see README.md):

    u/Ubar = [1 - cosh(Ha y/h)/cosh(Ha)] / [1 - tanh(Ha)/Ha]
    f*Re   = 32 Ha tanh(Ha) / [1 - tanh(Ha)/Ha],   Re = Ubar*Dh/nu, Dh = 4h
    Ha     = h B0 sqrt(sigma/(rho nu))

Limit: Ha -> 0 gives f*Re = 96 (plane Poiseuille).

Usage: python3 make_case.py [directory]
"""
import math
import os
import sys

L, H, T = 0.05, 0.05, 0.002        # length, half-width, thickness [m]
NY = 100                           # cells across the half-width
G = 0.0105                         # imposed pressure gradient [m/s^2]
RHO, NU, SIGMA = 1000.0, 4.1792e-4, 70.9
MU0 = 1.2566e-6
HB = math.sqrt(SIGMA / (RHO * NU))  # Ha = B0 * H * HB
DEST = sys.argv[1] if len(sys.argv) > 1 else 'channelHartmann'


def write(path, text):
    p = os.path.join(DEST, path)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w').write(text)


BLOCKMESH = """FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
convertToMeters 1;
vertices
(
    (0 -%(H)s 0) (%(L)s -%(H)s 0) (%(L)s %(H)s 0) (0 %(H)s 0)
    (0 -%(H)s %(T)s) (%(L)s -%(H)s %(T)s) (%(L)s %(H)s %(T)s) (0 %(H)s %(T)s)
);
blocks
(
    hex (0 1 2 3 4 5 6 7) (1 %(NY)d 1) simpleGrading (1 1 1)
);
edges ();
boundary
(
    left    { type cyclic; neighbourPatch right; faces ((0 4 7 3)); }
    right   { type cyclic; neighbourPatch left;  faces ((1 2 6 5)); }
    bottom  { type wall;   faces ((0 1 5 4)); }
    top     { type wall;   faces ((3 7 6 2)); }
    front   { type empty;  faces ((0 3 2 1)); }
    back    { type empty;  faces ((4 5 6 7)); }
);
mergePatchPairs ();
""" % dict(H=H, L=L, T=T, NY=NY)

FVSCHEMES = """FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; div(phi,U) Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
"""

FVSOLUTION = """FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers
{
    p      { solver GAMG; smoother GaussSeidel; tolerance 1e-11; relTol 0; }
    pFinal { solver GAMG; smoother GaussSeidel; tolerance 1e-11; relTol 0; }
    PotE { solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0; }
    U    { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-11; relTol 0; }
}
PISO { pRefCell 0; pRefValue 0; nCorrectors 3; nNonOrthogonalCorrectors 0; }
PotE { PotERefCell 0; PotERefValue 0; }
"""

CONTROL = """FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }
application     mhdturbFoamQS;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         80;
deltaT          0.2;
writeControl    timeStep;
writeInterval   1000;
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
runTimeModifiable true;
functions
{
    volFieldValue1
    {
        type            volFieldValue;
        libs            ( "libfieldFunctionObjects.so" );
        writeControl    timeStep;
        writeInterval   20;
        log             true;
        writeFields     false;
        operation       volAverage;
        weightField     none;
        fields          ( U p );
    }
}
"""

FVOPTIONS = """FoamFile { version 2.0; format ascii; class dictionary; object fvOptions; }
momentumSource
{
    type            vectorSemiImplicitSource;
    selectionMode   all;
    fields          (U);
    volumeMode      specific;
    injectionRateSuSp { U ((%(G)s 0 0) 0); }
}
""" % dict(G=G)

TRANSPORT = """FoamFile { version 2.0; format ascii; class dictionary; object transportProperties; }
transportModel  Newtonian;
rho     rho [1 -3 0 0 0 0 0] %(RHO)s;
nu      nu  [0 2 -1 0 0 0 0] %(NU)s;
mu      mu  [1 1 -2 0 0 -2 0] %(MU0)s;
sigma   sigma [-1 -3 3 0 0 2 0] %(SIGMA)s;
""" % dict(RHO=RHO, NU=NU, MU0=MU0, SIGMA=SIGMA)

TURBULENCE = """FoamFile { version 2.0; format ascii; class dictionary; object turbulenceProperties; }
simulationType  laminar;
"""


def field(name, cls, dims, internal, wall_bc, extra=''):
    """Build a field file with the patch layout of this case."""
    return """FoamFile { version 2.0; format ascii; class %(cls)s; object %(name)s; }
dimensions      %(dims)s;
internalField   uniform %(internal)s;
boundaryField
{
    left    { type cyclic; }
    right   { type cyclic; }
    bottom  { %(wall_bc)s }
    top     { %(wall_bc)s }
    front   { type empty; }
    back    { type empty; }
%(extra)s}
""" % dict(cls=cls, name=name, dims=dims, internal=internal,
           wall_bc=wall_bc, extra=extra)


U_FIELD = field('U', 'volVectorField', '[0 1 -1 0 0 0 0]', '(0 0 0)',
                'type noSlip;')
P_FIELD = field('p', 'volScalarField', '[0 2 -2 0 0 0 0]', '0',
                'type zeroGradient;')
POTE_FIELD = field('PotE', 'volScalarField', '[1 2 -3 0 0 -1 0]', '0',
                   'type zeroGradient;')


def b0_field(b):
    return field('B0', 'volVectorField', '[1 0 -2 0 0 -1 0]', '(0 %s 0)' % b,
                 'type fixedValue; value uniform (0 %s 0);' % b)


def main():
    os.makedirs(DEST, exist_ok=True)
    write('system/blockMeshDict', BLOCKMESH)
    write('system/fvSchemes', FVSCHEMES)
    write('system/fvSolution', FVSOLUTION)
    write('system/controlDict', CONTROL)
    write('system/fvOptions', FVOPTIONS)
    write('constant/transportProperties', TRANSPORT)
    write('constant/turbulenceProperties', TURBULENCE)
    write('0/U', U_FIELD)
    write('0/p', P_FIELD)
    write('0/PotE', POTE_FIELD)
    write('0/B0', b0_field(0.0))
    print("case written to", DEST)
    print("Ha = h B0 sqrt(sigma/(rho nu)) = B0 * %.6f" % (H * HB))
    for ha in (0.0, 0.5, 1.0, 5.0, 10.0):
        b = ha / (H * HB)
        t = math.tanh(ha)
        fRe = 96.0 if ha == 0 else 32.0 * ha * t / (1.0 - t / ha)
        print("   Ha=%5.2f  B0=%9.4f T   f*Re(exact)=%9.4f" % (ha, b, fRe))
    print("   G = %.6f m/s2 ; Dh = %.4f m" % (G, 4 * H))


if __name__ == '__main__':
    main()
