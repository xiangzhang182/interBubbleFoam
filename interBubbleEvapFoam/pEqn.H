{


    //Instead, use rAUc and rAUcf created in UEqn
//    volScalarField rAU("rAU", 1.0/UEqn.A());
//    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));


    volVectorField HbyA(constrainHbyA(rAUc*UEqn.H(), U, p_rgh));
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + alphacf*fvc::interpolate(rho*rAUc)*fvc::ddtCorr(U, phi)
    );
    
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUcf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUcf);

    Pair<tmp<volScalarField>> vDot = mixture->vDot();
    const volScalarField& vDotc = vDot[0]();
    const volScalarField& vDotv = vDot[1]();

    while (pimple.correctNonOrthogonal())
    {
    
       surfaceScalarField Dp("Dp", alphacf*rAUcf);
    
       fvScalarMatrix p_rghEqn
       (
           fvm::laplacian(Dp, p_rgh)
          ==
            ddt_alphac
          + fvc::div(alphacf*phiHbyA)
          + vDotv + vDotc

        );


/*
//Rattner - block of sanity checks

Info<< endl << "================================" << endl;
//Total liquid mass flow in:
const fvBoundaryMesh& bMesh = mesh.boundary();
//Inlet patch
label TopPatchID  = bMesh.findPatchID("Top");
const fvPatch& TopPatch = bMesh[TopPatchID];

const scalarField phi_out = phi.boundaryField()[TopPatchID];
Info<<"V_dot total out: " << gSum(phi_out) << endl;

const scalar TotalEvapVolumeChange = gSum( vDotv * mesh.V().field() );
Info<<"Total Evap Volume change: " << TotalEvapVolumeChange << endl;
const scalar TotalCondVolumeChange = gSum( vDotc * mesh.V().field() );
Info<<"Total Cond Volume change: " << TotalCondVolumeChange << endl;


Pair<tmp<volScalarField>> mDotPC = mixture->mDot();
const volScalarField& mDote = mDotPC[1]();

const scalar TotalLiquidVolumeEvaporated = gSum( mDote.field() * mesh.V().field()  ) / rho1.value();
Info<<"Total Volume of Liquid Evaporated: " << TotalLiquidVolumeEvaporated << endl;


//Info<< "VDotv = " << gMax(vDotv.field() + vDotc.field()) << endl;
vDotNet = vDotv + vDotc;

Info<< "================================" << endl << endl;
*/



        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux()/alphacf;

            U = HbyA + rAUc*fvc::reconstruct((phig - p_rghEqn.flux()/alphacf)/rAUcf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    p = p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
