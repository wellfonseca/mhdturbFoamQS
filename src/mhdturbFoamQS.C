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
    Incompressible turbulent MHD solver in the QUASI-STATIC (low-Rm)
    approximation, using the standard OpenFOAM turbulence models:

        ddt(U) + div(phi,U) + divDevReff(U) = -grad(p) + (1/rho) j x B0
        div(U) = 0
        lap(PotE) = div(u x B0)          (charge conservation)
        j = sigma (-grad(PotE) + u x B0)

    The Lorentz force is obtained directly from the current density, NOT from
    the Maxwell stress of the total field. This matters when the applied field
    B0 is NON-UNIFORM (e.g. alternating magnets): in the full-induction
    formulation the two Maxwell-stress terms are O(B0^2/mu0) and their discrete
    cancellation is no longer exact, leaving a spurious force of order
    (dx/L)^2 B0^2/(mu0 L rho) that dominates the physical force.

    Electrically insulating walls: j.n = 0 on all boundaries, implemented by
    excluding the boundary flux of u x B0.

    Valid for Rm = mu0 sigma U L << 1.

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
        "Incompressible turbulent MHD solver in the quasi-static (low-Rm)"
        " approximation, using the standard OpenFOAM turbulence models."
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

        // Lorentz force from the potential of the previous step (predictor)
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

        // --- electric potential (quasi-static)
        {
            surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();

            // Insulating walls: j.n = 0. The boundary flux of u x B0 is
            // excluded from BOTH sides of the potential equation: from the
            // Laplacian (PotE is zeroGradient) and from the right-hand side.
            // This gives j.n = -snGrad(PotE).n + (u x B0).n = 0.
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
