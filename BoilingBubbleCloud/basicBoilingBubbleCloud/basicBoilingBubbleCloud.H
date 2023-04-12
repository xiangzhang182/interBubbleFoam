/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::basicBoilingBubbleCloud

Description
    Cloud class to introduce popping bubble parcels

\*---------------------------------------------------------------------------*/

#ifndef basicBoilingBubbleCloud_H
#define basicBoilingBubbleCloud_H

#include "Cloud.H"
#include "KinematicCloud.H"
#include "BoilingBubbleCloud.H"
#include "CollidingCloud.H"
#include "basicBoilingBubbleParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef BoilingBubbleCloud
    <
        CollidingCloud
        <        
            KinematicCloud
	        <
	            Cloud
	            <
	               basicBoilingBubbleParcel
	            >
	        >
	    >
	> basicBoilingBubbleCloud;
	
	typedef basicBoilingBubbleCloud::collidingCloudType collidingBubbleType;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
