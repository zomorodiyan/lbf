/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    laserbeamFoam

Description
    Ray-Tracing heat source implementation with two phase incompressible VoF
    description of the metallic substrate and shielding gas phase.

Authors
    Tom Flint, UoM.
    Philip Cardiff, UCD.
    Gowthaman Parivendhan, UCD.
    Joe Robson, UoM.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interfaceCompression.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "noPhaseChange.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "Polynomial.H"
#include "laserHeatSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        // Info << "DEBUG: Start of time loop" << endl;
        #include "readControls.H"
        // Info << "DEBUG: After readControls.H" << endl;
        #include "readDyMControls.H"
        // Info << "DEBUG: After readDyMControls.H" << endl;

        if (LTS)
        {
            #include "setRDeltaT.H"
            // Info << "DEBUG: After setRDeltaT.H" << endl;
        }
        else
        {
            #include "CourantNo.H"
            // Info << "DEBUG: After CourantNo.H" << endl;
            #include "alphaCourantNo.H"
            // Info << "DEBUG: After alphaCourantNo.H" << endl;
            #include "setDeltaT.H"
            // Info << "DEBUG: After setDeltaT.H" << endl;
        }

        fvModels.preUpdateMesh();
        // Info << "DEBUG: After fvModels.preUpdateMesh()" << endl;

        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        tmp<volScalarField> divU;

        if
        (
            correctPhi
         && !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
         && mesh.topoChanged()
        )
        {
            // Construct and register divU for correctPhi
            divU = new volScalarField
            (
                "divU0",
                fvc::div(fvc::absolute(phi, U))
            );
        }

        // Update the mesh for topology change, mesh to mesh mapping
        bool topoChanged = mesh.update();

        // Do not apply previous time-step mesh compression flux
        // if the mesh topology changed
        if (topoChanged)
        {
            talphaPhi1Corr0.clear();
        }

        runTime++;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // Info << "DEBUG: Start of PIMPLE loop" << endl;
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                if
                (
                    correctPhi
                 && !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
                 && !divU.valid()
                )
                {
                    // Construct and register divU for correctPhi
                    divU = new volScalarField
                    (
                        "divU0",
                        fvc::div(fvc::absolute(phi, U))
                    );
                }

                // Move the mesh
                mesh.move();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"

                        // Update rhoPhi
                        rhoPhi = fvc::interpolate(rho)*phi;
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }

                divU.clear();
            }

            fvModels.correct();
            // Info << "DEBUG: After fvModels.correct()" << endl;

            #include "alphaControls.H"
            // Info << "DEBUG: After alphaControls.H" << endl;
            #include "alphaEqnSubCycle.H"
            // Info << "DEBUG: After alphaEqnSubCycle.H" << endl;

            #include "updateProps.H"
            // Info << "DEBUG: After updateProps.H" << endl;

            // Update the laser deposition field
            // laser.updateDeposition
            // (
            //     alpha_filtered, n_filtered, electrical_resistivity
            // );

            // --- Use Gaussian heat source instead ---
            {
                const scalar time = runTime.value();
                forAll(laser.laserNames(), laserI)
                {
                    const word& laserName = laser.laserNames()[laserI];
                    vector currentLaserPosition = laser.timeVsLaserPosition()[laserI](time);
                    scalar currentLaserPower = laser.timeVsLaserPower()[laserI](time);

                    const dictionary& dict = laser.laserDicts()[laserI];
                    scalar laserRadius = 0.0;
                    if (dict.found("HS_a"))
                    {
                        laserRadius = readScalar(dict.lookup("HS_a"));
                    }
                    else if (dict.found("laserRadius"))
                    {
                        laserRadius = readScalar(dict.lookup("laserRadius"));
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "The laser radius should be specified via 'laserRadius' or 'HS_a'"
                            << exit(FatalError);
                    }

                    // --- Add parameter validation and debug output ---
                    if (laserRadius <= SMALL || currentLaserPower <= SMALL)
                    {
                        FatalErrorInFunction
                            << "Invalid laser parameters in solver: "
                            << "laserRadius=" << laserRadius << ", "
                            << "currentLaserPower=" << currentLaserPower << nl
                            << "Check your LaserProperties dictionary!" << exit(FatalError);
                    }
                    // Info << "Laser " << laserName
                    //      << ": radius=" << laserRadius
                    //      << ", power=" << currentLaserPower
                    //      << ", position=" << currentLaserPosition << endl;

                    laser.updateGaussianDeposition
                    (
                        alpha_filtered,
                        laserName,
                        currentLaserPosition,
                        currentLaserPower,
                        laserRadius
                    );
                }
            }
            // Info << "DEBUG: After updateGaussianDeposition" << endl;

            // Fix: Only print Q if it exists in the objectRegistry
            if (mesh.objectRegistry::foundObject<volScalarField>("Q"))
            {
                volScalarField& Q = mesh.lookupObjectRef<volScalarField>("Q");
                // Info << "Q min: " << gMin(Q) << " max: " << gMax(Q) << endl;
            }
            else
            {
                // Info << "Q field not found in objectRegistry." << endl;
            }
            // Info << "T min: " << gMin(T) << " max: " << gMax(T) << endl;
            // Info << "U min: " << gMin(U) << " max: " << gMax(U) << endl;
            // Info << "p min: " << gMin(p) << " max: " << gMax(p) << endl;

            turbulence.correctPhasePhi();
            // Info << "DEBUG: After turbulence.correctPhasePhi()" << endl;

            mixture.correct();
            // Info << "DEBUG: After mixture.correct()" << endl;

            #include "UEqn.H"
            // Info << "DEBUG: After UEqn.H" << endl;

            #include "TEqn.H"
            // Info << "DEBUG: After TEqn.H" << endl;

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                // Info << "DEBUG: Start of pressure corrector loop" << endl;
                #include "pEqn.H"
                // Info << "DEBUG: After pEqn.H" << endl;
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
                // Info << "DEBUG: After turbulence.correct() (turbCorr)" << endl;
            }
            // Info << "DEBUG: End of PIMPLE loop" << endl;
        }

        // Check the cells that have melted
        // Info << "DEBUG: Before meltHistory update" << endl;
        volScalarField alphaMetal = 
            mesh.lookupObject<volScalarField>("alpha.metal");
        condition = pos(alphaMetal - 0.5) * pos(epsilon1 - 0.5);
        meltHistory += condition;
        // Info << "DEBUG: After meltHistory update" << endl;

        runTime.write();
        // Info << "DEBUG: After runTime.write()" << endl;

        // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        //     << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        //     << nl << endl;
        // Info << "DEBUG: End of time loop" << endl;
    }

    // Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
