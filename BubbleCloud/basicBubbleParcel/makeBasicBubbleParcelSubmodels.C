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

#include "basicBubbleCloud.H"

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

	makeParcelCloudFunctionObjects(basicBubbleCloud);

	// Kinematic sub-models
	makeParcelForces(basicBubbleCloud);

    //Additional special force for bubble pinning at the interface
    makeParticleForceModelType(InterfacePinningForce, basicBubbleCloud);

	makeParcelDispersionModels(basicBubbleCloud);
	makeParcelInjectionModels(basicBubbleCloud);

    makeParcelCollisionModels(collidingBubbleType);
	makeParcelPatchInteractionModels(basicBubbleCloud);
	makeParcelStochasticCollisionModels(basicBubbleCloud);
	makeParcelSurfaceFilmModels(basicBubbleCloud);

    // MPPIC sub-models
    makeMPPICParcelDampingModels(basicBubbleCloud);
    makeMPPICParcelIsotropyModels(basicBubbleCloud);
    makeMPPICParcelPackingModels(basicBubbleCloud);

// ************************************************************************* //
