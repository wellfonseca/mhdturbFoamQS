#!/usr/bin/env python3
"""
gera_canal_hartmann.py -- cria o caso de canal plano (Hartmann) para validar o
solver quase-estatico mhdturbFoamQS contra a solucao EXACTA.

Geometria: x in [0, L] (1 celula, ciclico), y in [-h, h] (paredes), z (1 celula, empty)
Campo B0 = (0, B0, 0) uniforme, PERPENDICULAR as paredes, paredes isolantes.
Escoamento forcado por gradiente de pressao fixo G (vectorSemiImplicitSource).

Solucao exacta (u(y) e f*Re):
    u/Ubar = [1 - cosh(Ha y/h)/cosh(Ha)] / [1 - tanh(Ha)/Ha]
    f*Re   = 32 Ha tanh(Ha) / [1 - tanh(Ha)/Ha],   Re = Ubar*D_h/nu, D_h = 4h
    Ha = h B0 sqrt(sigma/(rho nu))
Limites: Ha->0 -> f*Re = 96 (Poiseuille plano).

Uso: python3 gera_canal_hartmann.py [dir]
"""
import math
import os
import sys

L, H, T = 0.05, 0.05, 0.002          # meio-canal H, espessura T
NY = 100
G = 0.0105                            # gradiente de pressao [m/s2]
RHO, NU, SIGMA = 1000.0, 4.1792e-4, 70.9
MU0 = 1.2566e-6
HB = math.sqrt(SIGMA / (RHO * NU))
DEST = sys.argv[1] if len(sys.argv) > 1 else 'canal_hartmann'


def w(path, txt):
    p = os.path.join(DEST, path)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w').write(txt)


BM = """FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
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

FVSOL = """FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
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

FVO = """FoamFile { version 2.0; format ascii; class dictionary; object fvOptions; }
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

TURB = """FoamFile { version 2.0; format ascii; class dictionary; object turbulenceProperties; }
simulationType  laminar;
"""


def campo(nome, cls, dims, interno, paredes, extra=''):
    return """FoamFile { version 2.0; format ascii; class %(cls)s; object %(nome)s; }
dimensions      %(dims)s;
internalField   uniform %(interno)s;
boundaryField
{
    left    { type cyclic; }
    right   { type cyclic; }
    bottom  { %(paredes)s }
    top     { %(paredes)s }
    front   { type empty; }
    back    { type empty; }
%(extra)s}
""" % dict(cls=cls, nome=nome, dims=dims, interno=interno, paredes=paredes, extra=extra)


U = campo('U', 'volVectorField', '[0 1 -1 0 0 0 0]', '(0 0 0)', 'type noSlip;')
P = campo('p', 'volScalarField', '[0 2 -2 0 0 0 0]', '0', 'type zeroGradient;')
POTE = campo('PotE', 'volScalarField', '[1 2 -3 0 0 -1 0]', '0', 'type zeroGradient;')


def b0_file(B):
    return campo('B0', 'volVectorField', '[1 0 -2 0 0 -1 0]', '(0 %s 0)' % B,
                 'type fixedValue; value uniform (0 %s 0);' % B)


def main():
    os.makedirs(DEST, exist_ok=True)
    w('system/blockMeshDict', BM)
    w('system/fvSchemes', FVSCHEMES)
    w('system/fvSolution', FVSOL)
    w('system/controlDict', CONTROL)
    w('system/fvOptions', FVO)
    w('constant/transportProperties', TRANSPORT)
    w('constant/turbulenceProperties', TURB)
    w('0/U', U)
    w('0/p', P)
    w('0/PotE', POTE)
    w('0/B0', b0_file(0.0))
    print("caso criado em", DEST)
    print("Ha = h*B0*sqrt(sigma/(rho nu)) = B0 * %.6f" % (H * HB))
    for Ha in (0.0, 0.5, 1.0, 5.0, 10.0):
        b = Ha / (H * HB)
        t = math.tanh(Ha)
        fRe = 96.0 if Ha == 0 else 32.0 * Ha * t / (1.0 - t / Ha)
        print("   Ha=%5.2f  B0=%9.4f T   f*Re(exacto)=%9.4f" % (Ha, b, fRe))
    print("   G = %.6f m/s2 ; D_h = %.4f m" % (G, 4 * H))


if __name__ == '__main__':
    main()
