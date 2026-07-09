{

  //https://www.cfd-online.com/Forums/openfoam-solving/251323-teqn-compressibleinterfoam.html





    int iter = 0;


    dimensionedScalar residual("residual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual("meanResidual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar residual2("residual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual2("meanResidual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));





    TRHS*=0.0;//=TRHSzero;



    do{
iter++;
epsilon1.storePrevIter();



TRHS=          LatentHeat*((
                fvc::ddt(rho,epsilon1)
              + fvc::div(mixture.rhoPhi(), epsilon1)
              )


);

TRHS.correctBoundaryConditions();





    // Compressible work terms (p*div(U), kinetic energy) removed:
    // gas phases use rhoConst EOS, so Ma << 1 and these terms are negligible.
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(mixture.rhoPhi(), T)
      - fvm::laplacian(mixture.alphaEff(turbulence->mut()), T)

        ==
          mixture.rCv()*laser.deposition()
        - mixture.rCv()*TRHS
        - mass_dot
    );

    TEqn.relax();
    TEqn.solve();


Tcorr=(TLiquidus-TSolidus)*epsilon1+TSolidus;
Tcorr.correctBoundaryConditions();


epsilon1=max(min(epsilon1+(epsilonRel*mixture.Cv()/LatentHeat)*(T-Tcorr),scalar(1)),scalar(0));
epsilon1.correctBoundaryConditions();


residual = gMax(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField()));
meanResidual = gSum(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField())*mesh.V())/gSum(mesh.V().field());


        Info << "Correcting epsilon1, mean residual = " << meanResidual.value()
             << ", max residual = " << residual.value()
             << endl;


}while((iter < minTCorr || (residual > epsilonTol)) && iter <= maxTCorr );




    // Clip the entire T field to Tvap of Cu (2835 K) — applied globally.
    //
    // Rationale: this solver has no vapour phases and no recoil pressure. When
    // T would exceed Tvap, the real physics is vaporisation + vapour plume
    // formation, which we deliberately omit (GTAW conduction mode, no keyhole).
    // Without a latent heat sink at Tvap, the energy equation drives T
    // arbitrarily high, producing:
    //   (a) huge Marangoni gradients at the metal surface (dσ/dT × ΔT → large),
    //       causing unphysical velocities (80–170 m/s observed),
    //   (b) overheated air cells (rhoConst + tiny ρCp → T runaway from stray
    //       ray energy), reinforcing (a).
    // Capping at Tvap_Cu = 2835 K mimics the latent heat plateau without adding
    // vapour phases. Any excess energy above 2835 K is discarded — this is
    // physically justified because in conduction mode T should never reach Tvap.
    // If the cap is active on metal cells, it means the power is too high and
    // should be reduced; the cap prevents numerical blowup in the meantime.
    {
        const dimensionedScalar TvapCu("TvapCu", dimTemperature, 2835.0);
        T = min(T, TvapCu);
        T.correctBoundaryConditions();
    }

    mixture.correct();
}
