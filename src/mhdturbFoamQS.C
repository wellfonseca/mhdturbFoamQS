/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  6
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    mhdturbFoamQS

Description
    Solver MHD turbulento incompressivel na APROXIMACAO QUASE-ESTATICA
    (baixo Rm), com modelos de turbulencia standard:

        ddt(U) + div(phi,U) + divDevReff(U) = -grad(p) + (1/rho) j x B0
        div(U) = 0
        lap(PotE) = div(u x B0)          (conservacao da carga)
        j = sigma (-grad(PotE) + u x B0)

    A forca de Lorentz e' obtida directamente da densidade de corrente, e NAO
    pela tensao de Maxwell do campo total. Isso e' essencial quando o campo
    imposto B0 e' NAO-UNIFORME (imas alternados): na formulacao de inducao
    completa os dois termos da tensao de Maxwell sao O(B0^2/mu0) e o seu
    cancelamento discreto deixa de ser exacto, gerando uma forca espuria
    ~ (dx/L)^2 B0^2/(mu0 L rho) que domina a forca fisica.

    Paredes ELECTRICAMENTE ISOLANTES: j.n = 0 em todas as fronteiras, o que
    se implementa excluindo o fluxo de u x B0 nas faces de fronteira.

    Valido para Rm = mu0 sigma U L << 1 (aqui Rm ~ 3.7e-4).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "MHD turbulento incompressivel, aproximacao quase-estatica (baixo Rm),"
        " com modelos de turbulencia standard."
    );

    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // forca de Lorentz com o potencial do passo anterior (predictor)
        volVectorField lorentz
        (
            "lorentz",
            sigma*(-fvc::grad(PotE) ^ B0) + sigma*((U ^ B0) ^ B0)
        );

        #include "CourantNo.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          + turbulence->divDevReff(U)
          - (1.0/rho)*lorentz
         ==
            fvOptions(U)
        );

        UEqn.relax();
        fvOptions.constrain(UEqn);

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        {
            volScalarField rAU(1.0/UEqn.A());
            surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + rAUf*fvc::ddtCorr(U, phi)
            );

            constrainPressure(p, U, phiHbyA, rAUf);

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve(mesh.solver(p.select(piso.finalNonOrthogonalIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"
            p.relax();

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }

        // --- potencial electrico (quase-estatico)
        {
            surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();

            // paredes isolantes: j.n = 0 -> o fluxo de u x B0 na fronteira nao
            // entra nem no laplaciano (PotE com zeroGradient) nem no segundo
            // membro, o que da j.n = -snGrad(PotE).n + (u x B0).n = 0.
            surfaceScalarField psiubInt("psiubInt", psiub);
            forAll(psiubInt.boundaryFieldRef(), patchi)
            {
                psiubInt.boundaryFieldRef()[patchi] = 0.0;
            }

            fvScalarMatrix PotEEqn
            (
                fvm::laplacian(PotE) == fvc::div(psiubInt)
            );
            PotEEqn.setReference(potERefCell, potERefValue);
            PotEEqn.solve();

            surfaceScalarField jn
            (
                "jn",
                -(fvc::snGrad(PotE)*mesh.magSf()) + psiubInt
            );
            surfaceVectorField jnv("jnv", jn*mesh.Cf());
            volVectorField jfinal
            (
                "jfinal",
                fvc::surfaceIntegrate(jnv)
              - (fvc::surfaceIntegrate(jn)*mesh.C())
            );
            jfinal.correctBoundaryConditions();

            lorentz = sigma*(jfinal ^ B0);
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
