/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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
    interBubbleFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
//#include "CorrectPhi.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "fvcSmooth.H"

/*Addition of include files from reactingParcelFoam solver*/
//#include "turbulentFluidThermoModel.H"
//#include "basicKinematicCollidingCloud.H"
#include "basicBubbleCloud.H"
/*
#include "surfaceFilmModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pressureControl.H"
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
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

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

    	mu = turbulence->muEff();

        //CHW - Addition of command to add bubbles from reactingParcelFoam
        bubbles.evolve();

        //Rate of change (increase) of alpha_c due to bubbles popping and producing continuous gas phase
        volScalarField DDT_Popped("DDT_Popped", bubbles.VolPopped()/runTime.deltaT()  );
        DDT_Popped.ref() /= mesh.V();


        //CHW - Addition of update of continuous phase volume fraction field based on DPMFoam.C
        alphac = max(1.0 - bubbles.theta(), alphacMin);
        alphac.correctBoundaryConditions();

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

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

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
        
        
        runTime.write();

        runTime.printExecutionTime(Info);

    }
    
    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //


/*
	    scalar&       Dh,     // initial horizontal bubble diameter 0.568     mm   (from experiment)
	    scalar&       Dv,     // initial vertical bubble diameter   0.568     mm
	    const scalar  sigma,  // surface tension                    19.7      g/s2
	    scalar&       deltaDv,// vertical diameter change           0.05      mm   (Guess)
	    const scalar  n,      // film interface mobility            2 
	    const scalar  nu,     // viscosity                          0.000006  m2/s
	    const scalar  rho,    // liquid density                     913       kg/m3
	    const scalar  h0,     // initial film thickness             0.1       mm   (Guess)
        const scalar  Cm,     // mass coefficient                   0.85   
        vector&       Ub,     // bubble impact velocity             0.0024    m/s  (from experiment)

        scalar&       Rf,     // liquid film radius                 
        scalar&       Ek,     // kinetic energy
        scalar&       tcal,   // calculated liquid film drainage time
        scalar&       deltaP, // capillary presure
        scalar&       h,      // liquid film thickness
        vector&       Fdyn,   // temporary force 
        scalar&       Vb,     // bubble volume
        
        scalar Rf     =  sqr((Ek*Dh)/(2*3.14*sigma*deltaDv));
        scalar Ek     =  0.5*Cm*rho*Vb*pow2(Ub);
        scalar tcal   =  (3*nu*rho*pow2(Rf)*(1/pow2(h)-1/pow2(h0)))/(4*n*deltaP);
        scalar deltaP =  Fdyn/(3.14*pow2(Rf));
        scalar h      =  
        scalar Fdyn   =  Ek/deltaDv;
        scalar Vb     =  1.33*3.14*pow2(0.5*Dh)*0.5*Dv; 

   

*/
