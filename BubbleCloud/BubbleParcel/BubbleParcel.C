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

#include "BubbleParcel.H"


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::BubbleParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    tetIndices tetIs = this->currentTetIndices();

    td.alpha_L() = td.alpha_LInterp().interpolate(this->coordinates(), tetIs);
    
    td.n_hat() = td.n_hatInterp().interpolate(this->coordinates(), tetIs);
}





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::BubbleParcel<ParcelType>::BubbleParcel
(
    const BubbleParcel<ParcelType>& p
)
:
    ParcelType(p),
    intTime_(p.intTime_)
{}


template<class ParcelType>
Foam::BubbleParcel<ParcelType>::BubbleParcel
(
    const BubbleParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    intTime_(p.intTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::BubbleParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{

    ParcelType::calc(cloud, td, dt);

    if ( td.alpha_L() < 0.95 )          //   0.9
    {    this->intTime_ += dt; }
    else
    {    this->intTime_ = 0; }

}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "BubbleParcelIO.C"

// ************************************************************************* //
