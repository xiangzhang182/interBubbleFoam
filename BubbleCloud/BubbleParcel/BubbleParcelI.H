/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::BubbleParcel<ParcelType>::BubbleParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(owner, coordinates, celli, tetFacei, tetPti),
    intTime_(0.0)
{}


template<class ParcelType>
inline Foam::BubbleParcel<ParcelType>::BubbleParcel
(
    const polyMesh& owner,
    const vector& position,
    const label celli
)
:
    ParcelType(owner, position, celli),
    intTime_(0.0)
{}


template<class ParcelType>
inline Foam::BubbleParcel<ParcelType>::BubbleParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalar& intTime0,
    const typename ParcelType::constantProperties& constProps
)
:
    ParcelType
    (
        owner,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    intTime_(intTime0)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalar& Foam::BubbleParcel<ParcelType>::intTime() const
{
    return intTime_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleParcel<ParcelType>::intTime()
{
    return intTime_;
}


// ************************************************************************* //
