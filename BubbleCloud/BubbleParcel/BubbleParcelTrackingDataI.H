/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

template<class ParcelType>
template<class TrackCloudType>
inline Foam::BubbleParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),
    alpha_LInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.alpha_L()
        )
    ),
    n_hatInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.n_hat()
        )
    ),
    alpha_L_(Zero),
    n_hat_(0,0,0)
{
   
   
}

template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::BubbleParcel<ParcelType>::trackingData::alpha_LInterp() const
{
    return *alpha_LInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::vector>&
Foam::BubbleParcel<ParcelType>::trackingData::n_hatInterp() const
{
    return *n_hatInterp_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleParcel<ParcelType>::trackingData::alpha_L() const
{
    return alpha_L_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleParcel<ParcelType>::trackingData::alpha_L()
{
    return alpha_L_;
}

template<class ParcelType>
inline Foam::vector Foam::BubbleParcel<ParcelType>::trackingData::n_hat() const
{
    return n_hat_;
}


template<class ParcelType>
inline Foam::vector& Foam::BubbleParcel<ParcelType>::trackingData::n_hat()
{
    return n_hat_;
}



// ************************************************************************* //
