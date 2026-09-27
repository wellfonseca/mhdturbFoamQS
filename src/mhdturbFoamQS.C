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

        // Explicit part of the Lorentz force, evaluated from the fields of the
        // previous time step (predictor). The stiff sink -sigma*|B0|^2*U of the
        // complete force is deliberately left out: it is discretised
        // implicitly in UEqn (see below), which removes the explicit damping
        // limit dt < 2*rho/(sigma*|B0|^2) that otherwise forces a very small
        // time step at high Hartmann number. Using
        //     (U x B0) x B0 = B0*(U.B0) - |B0|^2*U
        // the explicit part is sigma*(-grad(PotE) x B0) + sigma*B0*(U.B0).
        volVectorField lorentz
        (
            "lorentz",
            sigma*(-fvc::grad(PotE) ^ B0) + sigma*B0*(U & B0)
        );

        #include "CourantNo.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          + turbulence->divDevReff(U)
          + fvm::Sp(sigma*magSqr(B0)/rho, U)
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
        //
        // The outer corrector loop (piso.correct(), nCorrectors in fvSolution)
        // is essential. rAU, HbyA and phiHbyA must be recomputed at EVERY
        // corrector; with a single corrector the pressure-velocity coupling is
        // only first order and the solution diverges even at Courant numbers
        // well below one. This was verified on the channelHartmann case, where
        // removing the loop makes the solution blow up at t = 3.4 s while the
        // standard structure below remains stable.
        while (piso.correct())
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
                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        fvOptions.correct(U);

        // --- electric potential (quasi-static)
        {
            surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();

            // Insulating walls: j.n = 0. The boundary flux of u x B0 is
            // excluded from BOTH sides of the potential equation: from the
            // Laplacian (PotE is zeroGradient on the walls) and from the
            // right-hand side. This gives j.n = -snGrad(PotE).n + (u x B0).n = 0.
            //
            // Only NON-coupled patches are zeroed. On cyclic and processor
            // patches the flux is an internal one and must be kept:
            // fvm::laplacian couples the two halves of a cyclic patch, so
            // dropping the same contribution from the right-hand side would
            // make the two sides of the equation inconsistent and put a
            // spurious source in the cells next to the periodic plane.
            surfaceScalarField psiubInt("psiubInt", psiub);
            forAll(psiubInt.boundaryFieldRef(), patchi)
            {
                if (!mesh.boundary()[patchi].coupled())
                {
                    psiubInt.boundaryFieldRef()[patchi] = 0.0;
                }
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

            // Explicit part of the force for the predictor of the next step:
            // the complete force from the reconstructed current, minus the
            // implicit sink already carried by UEqn.
            lorentz = sigma*(jfinal ^ B0) + sigma*magSqr(B0)*U;
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
