    fvVectorMatrix UEqn
    (
        fvm::ddt(alphacRho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      + turbulence->divDevRhoReff(rho, U)
      ==
      fvOptions(rho, U)
      + cloudSU
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);


    volScalarField rAUc(1.0/UEqn.A());
    surfaceScalarField rAUcf(fvc::interpolate(rAUc));


    surfaceScalarField phicForces
    (
        (fvc::interpolate(rAUc*cloudVolSUSu) & mesh.Sf())
    );


    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                phicForces/rAUcf
                +
                (
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
        
        fvOptions.correct(U);
    }
