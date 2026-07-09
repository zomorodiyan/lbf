FoamFile { version 2.0; format ascii; class volScalarField; object alpha.Cu; }
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    topWall     { type zeroGradient; }
    bottomWall  { type zeroGradient; }
    leftWall    { type zeroGradient; }
    rightWall   { type zeroGradient; }
    front       { type zeroGradient; }
    back        { type zeroGradient; }
}
