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
inline Foam::ThermalBubbleParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),
    
    Cp_(cloud.thermo().Cp()),   //ASR - note: original version copied from ThermoParcel has Cloud.thermo().thermo.Cp()
    kappa_(cloud.thermo().kappa()),
      
    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.T()
        )
    ),
    CpInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            Cp_
        )
    ),
    kappaInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            kappa_
        )
    ),
    Tc_(Zero),
    Cpc_(Zero)
{
   
   
}

template<class ParcelType>
inline const Foam::volScalarField&
Foam::ThermalBubbleParcel<ParcelType>::trackingData::Cp() const     
{
    return Cp_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::ThermalBubbleParcel<ParcelType>::trackingData::kappa() const   
{
    return kappa_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermalBubbleParcel<ParcelType>::trackingData::TInterp() const
{
    return *TInterp_;
}

template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermalBubbleParcel<ParcelType>::trackingData::CpInterp() const
{
    return *CpInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermalBubbleParcel<ParcelType>::trackingData::kappaInterp() const
{
    return *kappaInterp_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermalBubbleParcel<ParcelType>::trackingData::Tc() const
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermalBubbleParcel<ParcelType>::trackingData::Tc()
{
    return Tc_;
}

template<class ParcelType>
inline Foam::scalar Foam::ThermalBubbleParcel<ParcelType>::trackingData::Cpc() const
{
    return Cpc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermalBubbleParcel<ParcelType>::trackingData::Cpc()
{
    return Cpc_;
}

// ************************************************************************* //
