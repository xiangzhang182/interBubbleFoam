if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

//ddt_ac_sum = gSum( mag( fvc::ddt(alphac)().internalField() )() );
//Info<< "2.1 ddt(alphac) " << ddt_ac_sum << endl;
//Info<< "   nold times " << alphac.nOldTimes() << endl;
//Info<< "   nowalphac " << gSum( alphac.internalField() ) << endl;
//Info<< "   oldalphac " << gSum( alphac.oldTime().internalField() ) << endl;
//Info<< "   alphac time index " << alphac.timeIndex() << endl;
//Info<< "   global time index " << runTime.timeIndex() << endl;
    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
//ddt_ac_sum = gSum( mag( fvc::ddt(alphac)().internalField() )() );
//Info<< "2.2 ddt(alphac) " << ddt_ac_sum << endl;
//Info<< "   nold times " << alphac.nOldTimes() << endl;
//Info<< "   nowalphac " << gSum( alphac.internalField() ) << endl;
//Info<< "   oldalphac " << gSum( alphac.oldTime().internalField() ) << endl;
//Info<< "   alphac time index " << alphac.timeIndex() << endl;
//Info<< "   global time index " << runTime.timeIndex() << endl;
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;

//Info<< "After alphac time index " << alphac.timeIndex() << endl;
//Info<< "After p_rgh time index " << p_rgh.timeIndex() << endl;
//Info<< "After U time index " << U.timeIndex() << endl;
//Info<< " global time index " << runTime.timeIndex() << endl;
