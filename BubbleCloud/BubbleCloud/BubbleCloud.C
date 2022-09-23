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
    const immiscibleIncompressibleTwoPhaseMixture& mixture,
    bool  readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
    mixture_(mixture),
    alpha_L_(mixture.alpha1()),
    VolPopped_
    (
        IOobject
        (
            "VolPopped",
            rho.mesh().time().timeName(),
            rho.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho.mesh(),
        dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero)
    ),
    
    //To calculate pinning forces at the interface that keep the bubble from escaping, need access to interface normal field
    n_hat_
    (
        IOobject
        (
            "n_hat",
            rho.mesh().time().timeName(),
            rho.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho.mesh(),
        dimensionedVector("tmp",dimless,vector(0,0,0))
    )

    
    
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
    CloudType(c, name),
    mixture_(c.mixture()),
    alpha_L_(c.alpha_L()),
    VolPopped_
    (
        IOobject
        (
            "VolPopped",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero)
    ),
    n_hat_
    (
        IOobject
        (
            "n_hat",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedVector("tmp",dimless,vector(0,0,0))
    )
    
{}


template<class CloudType>
Foam::BubbleCloud<CloudType>::BubbleCloud
(
    const fvMesh& mesh,
    const word& name,
    const BubbleCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    mixture_(c.mixture()),
    alpha_L_(c.alpha_L()),
    VolPopped_
    (
        IOobject
        (
            "VolPopped",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero)
    ),
    n_hat_
    (
        IOobject
        (
            "n_hat",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedVector("tmp",dimless,vector(0,0,0))
    )
       
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
        //Calculate current interface normal field
        n_hat_ = mixture_.nearInterface()*fvc::grad(alpha_L_) / ( mag(fvc::grad(alpha_L_)) + mixture_.deltaN() );


        typename parcelType::trackingData td(*this);

        this->solve(*this, td);

        //Reset popped field to zero
        VolPopped_ = dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero);

        
        //Remove particles on the interface for too long
        for (parcelType& p : *this)
        {

            //Empirical fit for bubble popping time for test case.
            //In future, could be made more flexible/run time modifiable
            
            const scalar PopTime = ( 1668.8*( p.d() ) - 0.35795);

            if (p.intTime() > PopTime)
            {

               //Account for volume from popping
               VolPopped_[ p.cell() ] += p.volume();
               
          
               CloudType::deleteParticle(p);
               
            }
        
        }
        
        
    }
}


template<class CloudType>
void Foam::BubbleCloud<CloudType>::info()
{
    CloudType::info();

    //Todo - print out some relevant info
}


// ************************************************************************* //
