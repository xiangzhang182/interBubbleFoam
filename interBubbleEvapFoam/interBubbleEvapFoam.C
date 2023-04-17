/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interBubbleEvapFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, non-isothermal immiscible fluids with
    phase-change (evaporation-condensation) between a fluid and its vapour.
    Uses a VOF (volume of fluid) phase-fraction based interface capturing
    approach. Also with Lagrangian modeling of dispersed phase bubbles.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixtureEThermo.H"
#include "temperaturePhaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
//#include "CorrectPhi.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"

#include "basicThermalBubbleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, non-isothermal immiscible fluids with"
        " phase-change,"
        " using VOF phase-fraction based interface capturing.\n"
        "dispersed vapor bubbles tracked with Lagrangian method."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    volScalarField& T = thermo->T();

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence

        volScalarField divU("divU", fvc::div(fvc::absolute(phi, U)));

        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Update Lagrangian cloud
        mu = thermo->mu();
        bubbles.evolve(); //Update state of bubbles, motion, etc.

        //Continuous phases volume fraction
        alphac = max(1.0 - bubbles.theta(), alphacMin);
        alphac.correctBoundaryConditions();

        ddt_alphac = fvc::ddt(alphac);

        //Additional Lagrangian field calculations:
        Info<< "Continuous phase-1 volume fraction = "
            << alphac.weightedAverage(mesh.Vsc()).value()
            << "  Min(alphac) = " << min(alphac).value()
            << "  Max(alphac) = " << max(alphac).value()
            << endl;

        alphacf = fvc::interpolate(alphac);
        alphaRhoPhic = alphacf*rhoPhi;

        alphaPhic = alphacf*phi;
        alphacRho = alphac*rho;

        fvVectorMatrix cloudSU(bubbles.SU(U));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector(cloudSU.dimensions()/dimVolume, Zero),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();

        cloudSU.source() = vector::zero;



        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            mixture->correct();

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            interface.correct();

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = alpha1*rho1 + alpha2*rho2; //Why is this done again at the end of the loop?


        //Handle bubble popping
        //Rate of change (increase) of alpha_c due to bubbles popping and producing continuous gas phase
        bubbles.HandleBubblePopping();

        const volScalarField alphac_old = alphac;
        volScalarField VolPoppedSpecific = bubbles.VolPopped();
        VolPoppedSpecific.ref() /= mesh.V();
        alphac += VolPoppedSpecific;

        alpha1 *= alphac_old / alphac;



        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
