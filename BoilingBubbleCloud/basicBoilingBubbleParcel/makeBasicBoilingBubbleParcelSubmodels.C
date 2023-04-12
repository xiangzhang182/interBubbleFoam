/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "basicBoilingBubbleCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic sub-models
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelCollisionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeParcelSurfaceFilmModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"


#include "InterfacePinningForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	makeParcelCloudFunctionObjects(basicBoilingBubbleCloud);

	// Kinematic sub-models
	makeParcelForces(basicBoilingBubbleCloud);

    //Additional special force for bubble pinning at the interface
    makeParticleForceModelType(InterfacePinningForce, basicBoilingBubbleCloud);

	makeParcelDispersionModels(basicBoilingBubbleCloud);
	makeParcelInjectionModels(basicBoilingBubbleCloud);

    makeParcelCollisionModels(collidingBubbleType);
	makeParcelPatchInteractionModels(basicBoilingBubbleCloud);
	makeParcelStochasticCollisionModels(basicBoilingBubbleCloud);
	makeParcelSurfaceFilmModels(basicBoilingBubbleCloud);

    // MPPIC sub-models
    makeMPPICParcelDampingModels(basicBoilingBubbleCloud);
    makeMPPICParcelIsotropyModels(basicBoilingBubbleCloud);
    makeMPPICParcelPackingModels(basicBoilingBubbleCloud);

// ************************************************************************* //
