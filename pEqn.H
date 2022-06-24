{
    volVectorField HbyA(constrainHbyA(rAUc*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + alphacf*fvc::interpolate(rho*rAUc)*fvc::ddtCorr(U, phi)
    );

    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);


//    if (p_rgh.needReference())               
//    {
        fvc::makeRelative(phiHbyA, U);  //Note - this line, and the subsequent make absolute don't exist in the MPPICinterFoam, I'm not sure why, but it may have to do with dynamic mesh support?
//        adjustPhi(phiHbyA, U, p_rgh);        
//        fvc::makeAbsolute(phiHbyA, U);
//    }

//CHW-Revised surface tension force to not include bubbles
    surfaceScalarField phig
    (
        phicForces +
        (
            fvc::interpolate(mixture.sigmaK())*fvc::snGrad(alpha1) ///alphac)
          - ghf*fvc::snGrad(rho)
        )*rAUcf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUcf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        surfaceScalarField Dp("Dp", alphacf*rAUcf);
    
        // CHW - Added ddt(alphac) term from DPMFoam pEqn.H
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(Dp, p_rgh)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            //surfaceScalarField A = p_rghEqn.flux() / 0.5;
        
            phi = phiHbyA - p_rghEqn.flux()/alphacf;

            p_rgh.relax();

            U = HbyA + rAUc*fvc::reconstruct((phig - p_rghEqn.flux()/alphacf)/rAUcf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    //fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

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

//    if (!correctPhi)
//    {
//        rAUc.clear();
//    }
}
