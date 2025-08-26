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

#include "SortableList.H"
#include "RHFTable.H"

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
: IOdictionary
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
    effectiveRadius_(lookupOrDefault<scalar>("effectiveRadius", 0.001)),
    laserHeight_(lookupOrDefault<scalar>("laserHeight", 0.0001085)),
    absorptivity_(lookupOrDefault<scalar>("absorptivity", 0.48)),
    taperLength_(lookupOrDefault<scalar>("taperLength", 0.00001)),
    laserNames_(0),
    laserDicts_(0),
    timeVsLaserPosition_(0),
    timeVsLaserPower_(0),
    rhfTable_(),
    rhfTableFile_("")
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

    // Load RHF table if specified
    if (this->found("rhfTableFile"))
    {
        rhfTableFile_ = fileName(Foam::string(this->lookup("rhfTableFile")));
        rhfTable_.readCSV(rhfTableFile_);
    }
}
    

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar Foam::laserHeatSource::gaussianFunction(scalar r, scalar sigma) const
{
    return Foam::exp(-Foam::pow(r, 2.0) / (2.0 * Foam::pow(sigma, 2.0)));
}

scalar Foam::laserHeatSource::getRadiusAtPosition(scalar s, scalar effectiveRadius, scalar laserHeight, scalar taperLength) const
{
    if (s < 0 || s > laserHeight)
    {
        return 0.0;
    }

    scalar constantLength = laserHeight - taperLength;
    if (s <= constantLength)
    {
        return effectiveRadius;
    }
    else
    {
        // Linear taper to zero
        scalar taperProgress = (s - constantLength) / taperLength;
        return effectiveRadius * (1.0 - taperProgress);
    }
}

void Foam::laserHeatSource::toCartesian(scalar r, scalar theta, scalar s, scalar& x, scalar& y, scalar& z) const
{
    const scalar tiltAngleDeg = 5.0;
    const scalar tiltAngleRad = tiltAngleDeg * constant::mathematical::pi / 180.0;
    
    x = r * Foam::cos(theta);
    y = s * Foam::cos(tiltAngleRad) - r * Foam::sin(theta) * Foam::sin(tiltAngleRad);
    z = s * Foam::sin(tiltAngleRad) + r * Foam::sin(theta) * Foam::cos(tiltAngleRad);
}

scalar Foam::laserHeatSource::integrand(scalar r, scalar theta, scalar s, scalar effectiveRadius, scalar laserHeight, scalar taperLength, scalar laserRadius) const
{
    // Check if point is within object bounds
    scalar maxRadius = getRadiusAtPosition(s, effectiveRadius, laserHeight, taperLength);
    if (r > maxRadius)
    {
        return 0.0;
    }

    // Check if point is above y=0 surface (cut-off constraint)
    scalar x, y, z;
    toCartesian(r, theta, s, x, y, z);
    if (y < 0)
    {
        return 0.0;
    }

    // Calculate Gaussian value
    scalar sigma = laserRadius / 2.0;
    scalar gaussianVal = gaussianFunction(r, sigma);

    // Return with cylindrical Jacobian
    return gaussianVal * r;
}

scalar Foam::laserHeatSource::calculateGaussianIntegral(scalar effectiveRadius, scalar laserHeight, scalar taperLength, scalar laserRadius) const
{
    const label nS = 100;
    const label nR = 50;
    const label nTheta = 100;
    
    scalar integralSum = 0.0;
    const scalar pi = constant::mathematical::pi;
    
    // Create integration grids
    scalar ds = laserHeight / scalar(nS);
    scalar dtheta = 2.0 * pi / scalar(nTheta);
    
    for (label si = 0; si < nS; si++)
    {
        scalar s = -laserRadius + (scalar(si) + 0.5) * (laserHeight + laserRadius) / scalar(nS);
        scalar maxR = getRadiusAtPosition(s, effectiveRadius, laserHeight, taperLength);
        
        if (maxR > 0)
        {
            scalar dr = maxR / scalar(nR);
            
            for (label ri = 0; ri < nR; ri++)
            {
                scalar r = (scalar(ri) + 0.5) * dr;
                
                for (label thetai = 0; thetai < nTheta; thetai++)
                {
                    scalar theta = (scalar(thetai) + 0.5) * dtheta;
                    
                    // Check if point is above y=0 surface
                    scalar x, y, z;
                    toCartesian(r, theta, s, x, y, z);
                    if (y >= 0)
                    {
                        // Calculate integrand: Gaussian * r (cylindrical Jacobian)
                        scalar sigma = laserRadius / 2.0;
                        scalar gaussianVal = gaussianFunction(r, sigma);
                        scalar integrandVal = gaussianVal * r;
                        
                        // Add to integral with volume element
                        integralSum += integrandVal * dr * dtheta * ds;
                    }
                }
            }
        }
    }
    
    return integralSum;
}

void Foam::laserHeatSource::updateDeposition
(
    const volScalarField& /*alphaFiltered*/,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in
)
{
    // Ray-tracing is disabled. This is a stub to satisfy the linker.
}


void Foam::laserHeatSource::updateGaussianDeposition
(
    const volScalarField& /*alphaFiltered*/,
    const word& laserName,
    const vector& currentLaserPosition,
    const scalar currentLaserPower,
    const scalar laserRadius
)
{
    // Get current simulation time
    const scalar time = this->db().time().value();

    // Calculate and display the integral value for current parameters
    scalar integralValue = calculateGaussianIntegral(effectiveRadius_, laserHeight_, taperLength_, laserRadius);
    Info << " effectiveR: " << effectiveRadius_ << " Height: " << laserHeight_ << " taperLength: " << taperLength_ << " laserRadius:" << laserRadius << " Integral: " << integralValue << endl;

    // Get RHF value for this time (default to 1 if table not loaded)
    scalar rhf = 1.0;
    if (!rhfTableFile_.empty())
    {
        rhf = rhfTable_.value(time);
    }

    // Scale parameters by RHF squared (no enforced minima/maxima)
    const scalar rhf2 = Foam::pow(rhf, 2.0);

    const scalar scaledLaserRadius = laserRadius * rhf2;
    const scalar scaledLaserHeight = laserHeight_ * rhf2;
    const scalar scaledEffectiveRadius = effectiveRadius_ * rhf2;
    const scalar scaledAbsorptivity = absorptivity_ * rhf2;

    const scalar tiltAngleDeg = 5.0;
    const scalar tiltAngleRad = tiltAngleDeg * constant::mathematical::pi / 180.0;

    vector axisDir(Foam::sin(tiltAngleRad), Foam::cos(tiltAngleRad), 0.0);
    axisDir /= mag(axisDir);

    const fvMesh& mesh = deposition_.mesh();
    const vectorField& cellCenters = mesh.C();

    const scalar pi = constant::mathematical::pi;
    // Normalize so that the peak (center) value integrates to the total power over the cross-section
    // Use: Power / (pi * laserRadius^2 * laserHeight)
    const scalar gaussianNorm = 2 * currentLaserPower * scaledAbsorptivity / (pi * Foam::pow(scaledLaserRadius, 2.0) * scaledLaserHeight);
    Info << " gaussianNorm: " << gaussianNorm << endl;

    deposition_ *= 0.0;


    forAll(cellCenters, celli)
    {
        vector d = cellCenters[celli] - currentLaserPosition;
        scalar axial = d & axisDir;

        if (axial >= 0 && axial <= scaledLaserHeight)
        {
            // Calculate local radius for taper
            scalar localRadius = scaledLaserRadius;
            scalar localEffectiveRadius = scaledEffectiveRadius;
            if (taperLength_ > SMALL && axial > (scaledLaserHeight - taperLength_))
            {
                // Linear taper: radius decreases from full to zero over taperLength_
                scalar taperFrac = (scaledLaserHeight - axial) / taperLength_;
                taperFrac = Foam::max(0.0, Foam::min(1.0, taperFrac));
                localRadius = scaledLaserRadius * taperFrac;
                localEffectiveRadius = scaledEffectiveRadius * taperFrac;
            }

            vector radialVec = d - axial * axisDir;
            scalar radial2 = magSqr(radialVec);

            if (radial2 <= Foam::pow(localEffectiveRadius, 2.0))
            {
                // Use localRadius in Gaussian profile with standard form: exp(-r^2 / (2*sigma^2))
                scalar r = Foam::sqrt(radial2);
                scalar sigma = localRadius / 2.0;
                scalar Q = gaussianNorm * Foam::exp(-Foam::pow(r, 2.0) / (2.0 * Foam::pow(sigma, 2.0)));
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

} // End namespace Foam

// ************************************************************************* //
