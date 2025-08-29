{
    Info << "DEBUG: Entering TEqn.H" << endl;
    int iter = 0;
    scalarField residual(epsilon1.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;
    scalar Tresidual = 1.0;

    TRHS *= 0.0;

    Info << "DEBUG: Before rhoCp/rhophicp assignment" << endl;
    rhoCp = rho*cp;
    rhophicp = fvc::interpolate(cp)*rhoPhi;

    // const volTensorField gradU(fvc::grad(U));
    // volTensorField tau
    // (
    //     ((mixture.mu()*gradU) + (mixture.mu()*gradU.T()))
    //   - ((2.0/3.0)*mixture.mu()*(fvc::div(phi))*I)
    // );

    do
    {
        Info << "DEBUG: Start of TEqn iteration, iter=" << iter << endl;
        epsilon1.storePrevIter();
        T.storePrevIter();

        Info << "DEBUG: Before TRHS assignment" << endl;
        TRHS =
            LatentHeat*((fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)));

        Info << "DEBUG: Before Qv assignment" << endl;
        Qv =
            0.82*LatentHeatVap*Mm*p0
           *Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)))
           /Foam::pow(2*M_PI*Mm*R*T,0.5);

        Info << "DEBUG: Before damperSwitch block" << endl;
        if (damperSwitch)
        {
            scalarField& thermalDamperI = thermalDamper;
            const scalarField& rhoCpI = rhoCp;
            const scalarField& TI = T;
            forAll(thermalDamperI, celli)
            {
                thermalDamperI[celli] =
                    2.0*rhoCpI[celli]
                   /(
                       (rho1.value()*polycp_m.value(TI[celli]))
                     + (rho2.value()*polycp_g.value(TI[celli]))
                   );
            }

            thermalDamper.correctBoundaryConditions();
            Info << "DEBUG: After damperSwitch block" << endl;
        }

        Info << "DEBUG: Before TEqn matrix assembly" << endl;
        Info << "T min: " << gMin(T) << " max: " << gMax(T) << endl;
        Info << "rhoCp min: " << gMin(rhoCp) << " max: " << gMax(rhoCp) << endl;
        Info << "rhophicp min: " << gMin(rhophicp) << " max: " << gMax(rhophicp) << endl;
        Info << "kappa min: " << gMin(kappa) << " max: " << gMax(kappa) << endl;
        Info << "Qv min: " << gMin(Qv) << " max: " << gMax(Qv) << endl;
        Info << "gradAlpha min: " << gMin(gradAlpha) << " max: " << gMax(gradAlpha) << endl;
        Info << "thermalDamper min: " << gMin(thermalDamper) << " max: " << gMax(thermalDamper) << endl;
        Info << "TRHS min: " << gMin(TRHS) << " max: " << gMax(TRHS) << endl;
        Info << "laser.deposition() min: " << gMin(laser.deposition()) << " max: " << gMax(laser.deposition()) << endl;

        // Optionally, check for zeros or negative values in fields that should be positive
        label zeroRhoCp = 0, negRhoCp = 0;
        const scalarField& rhoCpF = rhoCp.primitiveField();
        forAll(rhoCpF, i)
        {
            if (rhoCpF[i] == 0.0) zeroRhoCp++;
            if (rhoCpF[i] < 0.0) negRhoCp++;
        }
        Info << "rhoCp zero count: " << zeroRhoCp << " negative count: " << negRhoCp << endl;

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
          //- (tau && gradU) // Viscous dissipation
         ==
            fvModels.source(rhoCp,T)
          + laser.deposition() // Gaussian heat source
          - Qv*mag(gradAlpha)*thermalDamper
          - TRHS
        );
        Info << "DEBUG: After TEqn matrix assembly" << endl;
        // Debug: check diagonal and source before relax
        label zeroDiag = 0, nanDiag = 0, infDiag = 0, negDiag = 0;
        const scalarField& diagF = TEqn.diag();
        forAll(diagF, i)
        {
            if (diagF[i] == 0.0) zeroDiag++;
            if (diagF[i] < 0.0) negDiag++;
            if (std::isnan(diagF[i])) nanDiag++;
            if (std::isinf(diagF[i])) infDiag++;
            if (i < 10) Info << "TEqn.diag()[" << i << "] = " << diagF[i] << endl;
        }
        Info << "TEqn diagonal zero count: " << zeroDiag
             << " negative count: " << negDiag
             << " NaN count: " << nanDiag
             << " Inf count: " << infDiag << endl;

        // Also check the source term for similar issues
        label zeroSrc = 0, nanSrc = 0, infSrc = 0, negSrc = 0;
        const scalarField& srcF = TEqn.source();
        forAll(srcF, i)
        {
            if (srcF[i] == 0.0) zeroSrc++;
            if (srcF[i] < 0.0) negSrc++;
            if (std::isnan(srcF[i])) nanSrc++;
            if (std::isinf(srcF[i])) infSrc++;
            if (i < 10) Info << "TEqn.source()[" << i << "] = " << srcF[i] << endl;
        }
        Info << "TEqn source zero count: " << zeroSrc
             << " negative count: " << negSrc
             << " NaN count: " << nanSrc
             << " Inf count: " << infSrc << endl;

        // Check cell volumes for zero or negative values
        label zeroVol = 0, negVol = 0;
        const scalarField& vF = mesh.V();
        forAll(vF, i)
        {
            if (vF[i] == 0.0) zeroVol++;
            if (vF[i] < 0.0) negVol++;
            if (i < 10) Info << "mesh.V()[" << i << "] = " << vF[i] << endl;
        }
        Info << "Cell volume zero count: " << zeroVol
             << " negative count: " << negVol << endl;

        TEqn.relax();
        Info << "DEBUG: After TEqn.relax()" << endl;

        fvConstraints.constrain(TEqn);
        Info << "DEBUG: After fvConstraints.constrain(TEqn)" << endl;

        TEqn.solve();
        Info << "DEBUG: After TEqn.solve()" << endl;

        T.relax();
        Info << "DEBUG: After T.relax()" << endl;

        Tcorr = (TLiquidus - TSolidus)*epsilon1 + TSolidus;
        Info << "DEBUG: After Tcorr assignment" << endl;

        epsilon1 =
            max
            (
                min
                (
                    epsilon1 + (epsilonRel*cp/LatentHeat)*(T - Tcorr),
                    scalar(1)
                ),
                scalar(0)
            );
        Info << "DEBUG: After epsilon1 update" << endl;

        residual =
            mag
            (
                epsilon1.primitiveField()
              - epsilon1.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);
        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        //ddte1 = fvc::ddt(epsilon1);

        Info<< "Correcting epsilon1, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;
    }
    while
    (
        (
            ++iter < minTCorr
         // || (meanResidual > epsilonTol.value())
         || (maxResidual > epsilonTol.value())
        )
     && iter <= maxTCorr
    );

    Info << "DEBUG: After TEqn iteration loop" << endl;
    gradT = fvc::grad(T);
    Info << "DEBUG: After gradT assignment" << endl;
}





