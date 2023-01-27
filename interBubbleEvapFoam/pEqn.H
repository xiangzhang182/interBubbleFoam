{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    Pair<tmp<volScalarField>> vDot = mixture->vDot();
    const volScalarField& vDotc = vDot[0]();
    const volScalarField& vDotv = vDot[1]();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          ==
            -vDotv - vDotc
        );

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


        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
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
