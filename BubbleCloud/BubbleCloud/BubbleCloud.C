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

#include "BubbleCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::BubbleCloud<CloudType>::setModels()
{


}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BubbleCloud<CloudType>::BubbleCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false)
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }
}


template<class CloudType>
Foam::BubbleCloud<CloudType>::BubbleCloud
(
    BubbleCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name)
{}


template<class CloudType>
Foam::BubbleCloud<CloudType>::BubbleCloud
(
    const fvMesh& mesh,
    const word& name,
    const BubbleCloud<CloudType>& c
)
:
    CloudType(mesh, name, c)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BubbleCloud<CloudType>::~BubbleCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BubbleCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<BubbleCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::BubbleCloud<CloudType>::restoreState()
{
    this->cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::BubbleCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::BubbleCloud<CloudType>::info()
{
    CloudType::info();

    //Todo - print out some relevant info
}


// ************************************************************************* //
