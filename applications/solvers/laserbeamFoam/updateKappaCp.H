{
    const scalarField& alpha1I = alpha1;
    const scalarField& TI = T;
    const scalarField& epsilon1I = epsilon1;
    scalarField& kappaI = kappa;
    scalarField& cpI = cp;
    const scalar TvapValue = Tvap.value();
    forAll(kappaI, celli)
    {
        // poly_kappa_liquid is only fitted up to ~4000K; clamp the argument
        // at Tvap so cells that numerically overshoot boiling hold at
        // kappa(Tvap) instead of following the polynomial into negative
        // values, which would make the diffusion term unconditionally
        // unstable and runaway further.
        const scalar TliquidKappa = min(TI[celli], TvapValue);

        kappaI[celli] =
            (
                max(min(alpha1I[celli], 1.0), 0.0)
               *(
                    useLiquidKappa
                  ? (1.0 - epsilon1I[celli])*polykappa_m.value(TI[celli])
                  + epsilon1I[celli]*polykappa_m_liquidPtr().value(TliquidKappa)
                  : polykappa_m.value(TI[celli])
                )
            )
          + (
                (1.0 - max(min(alpha1I[celli], 1.0), 0.0))
               *polykappa_g.value(TI[celli])
            );

        // Defense-in-depth: never let kappa reach zero/negative regardless
        // of which branch produced it.
        kappaI[celli] = max(kappaI[celli], 1.0);

        cpI[celli] =
            (
                max(min(alpha1I[celli], 1.0), 0.0)
               *polycp_m.value(TI[celli])
            )
          + (
                (1.0 - max(min(alpha1I[celli], 1.0), 0.0))
               *polycp_g.value(TI[celli])
            );
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();

    if (usePolySigma)
    {
        scalarField& dSigmadTI = dSigmadT;
        forAll(TI, celli)
        {
            dSigmadTI[celli] = poly_sigmaPtr->derivative(TI[celli]);
        }
        dSigmadT.correctBoundaryConditions();
    }
}
