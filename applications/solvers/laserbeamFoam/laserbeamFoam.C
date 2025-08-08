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

    while (pimple.run(runTime))
    {
        #include "readControls.H"
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

        runTime++;

        // --- Energy-only simulation: Only update properties and heat source ---
        #include "updateProps.H"

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

                if (laserRadius <= SMALL || currentLaserPower <= SMALL)
                {
                    FatalErrorInFunction
                        << "Invalid laser parameters in solver: "
                        << "laserRadius=" << laserRadius << ", "
                        << "currentLaserPower=" << currentLaserPower << nl
                        << "Check your LaserProperties dictionary!" << exit(FatalError);
                }

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

        if (mesh.objectRegistry::foundObject<volScalarField>("Q"))
        {
            volScalarField& Q = mesh.lookupObjectRef<volScalarField>("Q");
        }

        #include "TEqn.H"

        volScalarField alphaMetal = 
            mesh.lookupObject<volScalarField>("alpha.metal");
        condition = pos(alphaMetal - 0.5) * pos(epsilon1 - 0.5);
        meltHistory += condition;

        runTime.write();
    }

    return 0;
}


// ************************************************************************* //
