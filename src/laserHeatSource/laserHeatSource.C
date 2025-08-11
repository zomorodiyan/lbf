/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "laserHeatSource.H"
#include "fvc.H"
#include "constants.H"
#include "findLocalCell.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laserHeatSource, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laserHeatSource::laserHeatSource
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "LaserProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    deposition_
    (
        IOobject
        (
            "Deposition",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("deposition", dimensionSet(1, -1, -3, -0, 0), -1.0)
    ),
    laserBoundary_
    (
        IOobject
        (
            "Laser_boundary", // rename to laserBoundary?
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    errorTrack_
    (
        IOobject
        (
            "errorTrack",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "errorTrack", dimensionSet(0, 0, 0, -0, 0), 0.0
        )
    ),
    rayNumber_
    (
        IOobject
        (
            "rayNumber",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rayNumber",dimensionSet(0, 0, 0, -0, 0),-1.0)
    ),
    rayQ_
    (
        IOobject
        (
            "rayQ",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rayQ", dimensionSet(1, 0, -3, 0, 0), scalar(0.0))
    ),
    yDim_
    (
        IOobject
        (
            "yDim",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("yDim", dimensionSet(0, 1, 0, 0, 0), 1.0)
    ),
    refineFlag_
    (
        IOobject
        (
            "refineflag", // rename refineFlag
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("refineflag", dimensionSet(0,0,0,0,0), 0.0)
    ),
    powderSim_(lookupOrDefault<Switch>("PowderSim", false)),
    effectiveRadius_(lookupOrDefault<scalar>("effectiveRadius", 1e-3)),
    laserNames_(0),
    laserDicts_(0),
    timeVsLaserPosition_(0),
    timeVsLaserPower_(0)
{
    // Initialise the laser power and position
    if (found("lasers"))
    {
        const PtrList<entry> laserEntries(lookup("lasers"));

        laserNames_.setSize(laserEntries.size());
        laserDicts_.setSize(laserEntries.size());
        timeVsLaserPosition_.setSize(laserEntries.size());
        timeVsLaserPower_.setSize(laserEntries.size());

        forAll(laserEntries, laserI)
        {
            laserNames_[laserI] = laserEntries[laserI].keyword();

            laserDicts_.set(laserI, new dictionary(laserEntries[laserI].dict()));

            timeVsLaserPosition_.set
            (
                laserI,
                new interpolationTable<vector>
                (
                    laserEntries[laserI].dict().subDict("timeVsLaserPosition")
                )
            );

            timeVsLaserPower_.set
            (
                laserI,
                new interpolationTable<scalar>
                (
                    laserEntries[laserI].dict().subDict("timeVsLaserPower")
                )
            );
        }

        if (found("timeVsLaserPosition"))
        {
            FatalErrorInFunction
                << "timeVsLaserPosition should not be defined in the main dict"
                << " if a list of lasers is provided" << exit(FatalError);
        }

        if (found("timeVsLaserPower"))
        {
            FatalErrorInFunction
                << "timeVsLaserPower should not be defined in the main dict"
                << " if a list of lasers is provided" << exit(FatalError);
        }
    }
    else
    {
        laserNames_.setSize(1);
        laserDicts_.setSize(1);
        timeVsLaserPosition_.setSize(1);
        timeVsLaserPower_.setSize(1);

        laserNames_[0] = "laser0";

        laserDicts_.set(0, new dictionary(*this));

        timeVsLaserPosition_.set
        (
            0,
            new interpolationTable<vector>(subDict("timeVsLaserPosition"))
        );

        timeVsLaserPower_.set
        (
            0,
            new interpolationTable<scalar>(subDict("timeVsLaserPower"))
        );
    }

    laserBoundary_ = fvc::average(laserBoundary_);

    // Give errors if the old input format is found

    if (found("HS_bg"))
    {
        FatalErrorInFunction
            << "'HS_bg' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }

    if (found("HS_lg"))
    {
        FatalErrorInFunction
            << "'HS_lg' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }

    if (found("HS_velocity"))
    {
        FatalErrorInFunction
            << "'HS_velocity' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }

    if (found("HS_Q"))
    {
        FatalErrorInFunction
            << "'HS_Q' is deprecated: please instead specify the laser "
            << "power in time via the laserPowereVsTime sub-dict"
            << exit(FatalError);
    }

    if (found("elec_resistivity"))
    {
        FatalErrorInFunction
            << "'elec_resistivity' is deprecated: resistivity is now "
            << "passed in from the solver as a field"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// --- Restore empty stub for updateDeposition to fix linker error ---
void Foam::laserHeatSource::updateDeposition
(
    const volScalarField& alphaFiltered,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in
)
{
    // Ray-tracing is disabled. This is a stub to satisfy the linker.
}


void Foam::laserHeatSource::updateGaussianDeposition
(
    const volScalarField& alphaFiltered,
    const word& laserName,
    const vector& currentLaserPosition,
    const scalar currentLaserPower,
    const scalar laserRadius
)
{
    const scalar laserHeight = 0.0001;
    const scalar tiltAngleDeg = 5.0;
    const scalar tiltAngleRad = tiltAngleDeg * constant::mathematical::pi / 180.0;

    vector axisDir(Foam::sin(tiltAngleRad), Foam::cos(tiltAngleRad), 0.0);
    axisDir /= mag(axisDir);

    const fvMesh& mesh = deposition_.mesh();
    const vectorField& cellCenters = mesh.C();

    const scalar pi = constant::mathematical::pi;
    // Normalize so that the peak (center) value integrates to the total power over the cross-section
    // Use: Power / (pi * laserRadius^2 * laserHeight)
    const scalar gaussianNorm = 2 * currentLaserPower / (pi * Foam::pow(laserRadius, 2.0) * laserHeight);

    deposition_ *= 0.0;

    forAll(cellCenters, celli)
    {
        vector d = cellCenters[celli] - currentLaserPosition;
        scalar axial = d & axisDir;

        if (axial >= 0 && axial <= laserHeight)
        {
            vector radialVec = d - axial * axisDir;
            scalar radial2 = magSqr(radialVec);

            if (radial2 <= Foam::pow(effectiveRadius_, 2.0))
            {
                scalar Q = gaussianNorm * Foam::exp(-radial2 / Foam::pow(laserRadius, 2.0));
                deposition_[celli] = Q;
            }
            else
            {
                deposition_[celli] = 0.0;
            }
        }
        else
        {
            deposition_[celli] = 0.0;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
// ************************************************************************* //
