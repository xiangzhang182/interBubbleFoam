/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::basicThermalBubbleParcel

Description
    Definition of basic Bubble parcel

SourceFiles
    basicThermalBubbleParcel.H

\*---------------------------------------------------------------------------*/

#ifndef basicThermalBubbleParcel_H
#define basicThermalBubbleParcel_H

#include "contiguous.H"
#include "particle.H"
#include "KinematicParcel.H"
#include "CollidingParcel.H"
#include "ThermalBubbleParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ThermalBubbleParcel< BubbleParcel< CollidingParcel< KinematicParcel< particle > > > > basicThermalBubbleParcel;

    //- Contiguous data for basicThermalBubbleParcel
    template<> struct is_contiguous<basicThermalBubbleParcel> : std::true_type {};
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
