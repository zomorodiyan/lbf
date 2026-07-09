/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      thermophysicalProperties.Cu;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Commercially pure copper Grade 110 (Soysal & Kou 2016)
// Sources: Mills 2002, Brillo & Egry 2005
//

phaseState      condensed;

thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleInternalEnergy;
}

TSolidus            1356.0;   // K  (pure Cu melting point)
TLiquidus           1358.0;   // K  (nearly isothermal — pure metal)
LatentHeat          2.05e5;   // J/kg
elec_conductivity   5.0e6;    // S/m liquid
magnetic_permeability 1.256e-6;

mixture
{
    specie
    {
        molWeight   63.5;
    }
    equationOfState
    {
        rho     8960;         // kg/m3 constant density
    }
    thermodynamics
    {
        Cp      495;          // J/(kg.K) liquid
        Hf      0;
    }
    transport
    {
        mu      4.0e-3;       // Pa.s liquid viscosity
        Pr      0.03;
    }
}

// ************************************************************************* //
