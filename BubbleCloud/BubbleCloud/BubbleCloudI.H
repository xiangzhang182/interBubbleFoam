/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::BubbleCloud<CloudType>&
Foam::BubbleCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}

template<class CloudType>
inline const Foam::volScalarField& Foam::BubbleCloud<CloudType>::alpha_L() const
{
    return alpha_L_;
}

template<class CloudType>
inline const Foam::volScalarField& Foam::BubbleCloud<CloudType>::VolPopped() const
{
    return VolPopped_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::BubbleCloud<CloudType>::n_hat() const
{
    return n_hat_;
}

template<class CloudType>
inline const Foam::immiscibleIncompressibleTwoPhaseMixture& Foam::BubbleCloud<CloudType>::mixture() const
{
    return mixture_;
}

// ************************************************************************* //
