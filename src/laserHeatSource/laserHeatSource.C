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
            Info<< "Reading laser " << laserNames_[laserI] << endl;

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

        // Check that a single laser is not also defined

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
        // There is no lists of lasers, just one

        laserNames_.setSize(1);
        laserDicts_.setSize(1);
        timeVsLaserPosition_.setSize(1);
        timeVsLaserPower_.setSize(1);

        laserNames_[0] = "laser0";

        // Copy the main dict
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

    // Update laserBoundary
    laserBoundary_ = fvc::average(laserBoundary_);

    if (debug)
    {
        errorTrack_.writeOpt() = IOobject::AUTO_WRITE;
        rayNumber_.writeOpt() = IOobject::AUTO_WRITE;
    }

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
    // Reset deposition field
    deposition_ *= 0.0;

    const fvMesh& mesh = deposition_.mesh();
    const vectorField& cellCenters = mesh.C();
    const scalarField& alphaFilteredI = alphaFiltered;

    // Find laser dictionary index
    label laserI = -1;
    forAll(laserNames_, i)
    {
        if (laserNames_[i] == laserName)
        {
            laserI = i;
            break;
        }
    }
    if (laserI == -1)
    {
        FatalErrorInFunction << "Laser name not found: " << laserName << exit(FatalError);
    }
    const dictionary& dict = laserDicts_[laserI];

    // Read V_incident and laserHeight from dictionary
    vector V_incident(dict.lookup("V_incident"));
    scalar laserHeight = dict.lookupOrDefault<scalar>("laserHeight", laserRadius); // fallback if not present

    // --- Add parameter validation to avoid floating point exceptions ---
    if (laserRadius <= SMALL || laserHeight <= SMALL || currentLaserPower <= SMALL)
    {
        FatalErrorInFunction
            << "Invalid laser parameters: "
            << "laserRadius=" << laserRadius << ", "
            << "laserHeight=" << laserHeight << ", "
            << "currentLaserPower=" << currentLaserPower << nl
            << "Check your LaserProperties dictionary!" << exit(FatalError);
    }

    // Normalize direction
    vector beamDir = V_incident/mag(V_incident);

    // The center of the Gaussian is at the middle of the segment
    vector gaussianCenter = currentLaserPosition + 0.5 * laserHeight * beamDir;

    // 3D Gaussian normalization factor
    const scalar pi = constant::mathematical::pi;
    const scalar norm = currentLaserPower / (Foam::pow(pi, 1.5) * laserRadius * laserRadius * laserHeight);

    forAll(cellCenters, celli)
    {
        // Only deposit in cells with sufficient material (alphaFiltered > 0.5)
        if (alphaFilteredI[celli] > 0.5)
        {
            // Vector from Gaussian center to cell
            vector d = cellCenters[celli] - gaussianCenter;

            // Axial distance (along beam direction, centered at middle of segment)
            scalar axial = (d & beamDir);

            // Radial distance (perpendicular to beam direction)
            scalar radial2 = magSqr(d - axial*beamDir);

            // Only deposit within the segment [0, laserHeight] along beamDir starting from currentLaserPosition
            scalar proj = (cellCenters[celli] - currentLaserPosition) & beamDir;
            if (proj >= 0 && proj <= laserHeight)
            {
                // 3D Gaussian profile
                scalar Q = norm
                    * Foam::exp(-radial2 / Foam::pow(laserRadius, 2.0))
                    * Foam::exp(-Foam::pow(axial, 2.0) / Foam::pow(laserHeight, 2.0));

                deposition_[celli] += Q;
            }
        }
    }

    deposition_.correctBoundaryConditions();

    Info << "Total 3D Gaussian Q deposited this timestep: "
         << fvc::domainIntegrate(deposition_).value() << endl;

    Info << "deposition_ min: " << gMin(deposition_) << " max: " << gMax(deposition_) << endl;
    Info << "Gaussian norm: " << norm
         << " | laserRadius: " << laserRadius
         << " | laserHeight: " << laserHeight
         << " | currentLaserPower: " << currentLaserPower << endl;

    if (norm > 1e10)
    {
        WarningInFunction
            << "Gaussian normalization factor is extremely large: " << norm << nl
            << "This is likely due to very small laserRadius or laserHeight." << nl
            << "Check your LaserProperties dictionary for physical values." << endl;
    }

    if (fvc::domainIntegrate(deposition_).value() < 1e-12)
    {
        WarningInFunction
            << "Total Gaussian Q deposited is extremely small (" << fvc::domainIntegrate(deposition_).value() << ")"
            << " compared to the simulation time scale (1e-8 s)." << nl
            << "Check if your laser power, spot size, or time step are physically reasonable." << endl;
    }

    // --- Additional debug: check for NaN or Inf in deposition_ ---
    label nanCount = 0, infCount = 0;
    forAll(deposition_, celli)
    {
        if (Foam::isnan(deposition_[celli]))
        {
            nanCount++;
        }
        if (Foam::isinf(deposition_[celli]))
        {
            infCount++;
        }
    }
    Info << "deposition_ NaN count: " << nanCount << " | Inf count: " << infCount << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
