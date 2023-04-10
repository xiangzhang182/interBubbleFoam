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

#include "ThermalBubbleCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::setModels()
{


}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermalBubbleCloud<CloudType>::ThermalBubbleCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    const twoPhaseMixtureEThermo& thermo,
    const interfaceProperties& interface,
    bool  readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
    thermo_(thermo),
    cloudCopyPtr_(nullptr),
    interface_(interface),
    alpha_L_(thermo.alpha1()),    
    T_(thermo.T()),
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
    ),
    QTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":QTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime, Zero)
        )
    ),
    TCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":TCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime/dimTemperature, Zero)
        )
    ),
    VLTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":VLTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimVolume, Zero)
        )
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
    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ThermalBubbleCloud<CloudType>::ThermalBubbleCloud
(
    ThermalBubbleCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermo_(c.thermo()),
    cloudCopyPtr_(nullptr),
    interface_(c.interface()),
    alpha_L_(c.alpha_L()),
    T_(c.T()),
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
    ),
    QTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":QTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime/dimVolume, Zero)
        )
    ),
    TCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":TCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime/dimVolume/dimTemperature, Zero)
        )
    ),
    VLTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":VLTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimVolume, Zero)
        )
    )
    
{}


template<class CloudType>
Foam::ThermalBubbleCloud<CloudType>::ThermalBubbleCloud
(
    const fvMesh& mesh,
    const word& name,
    const ThermalBubbleCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    thermo_(c.thermo()),
    cloudCopyPtr_(nullptr),
    interface_(c.interface()),
    alpha_L_(c.alpha_L()),
    T_(c.T()),
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
    ),
    QTrans_(nullptr),
    TCoeff_(nullptr),
    VLTrans_(nullptr)
       
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermalBubbleCloud<CloudType>::~ThermalBubbleCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ThermalBubbleCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::restoreState()
{
    this->cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    QTrans_->field() = 0.0;
    TCoeff_->field() = 0.0;
    VLTrans_->field() = 0.0;
}

template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::relaxSources
(
    const ThermalBubbleCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    this->relax(QTrans_(), cloudOldTime.QTrans(), "T");
    this->relax(TCoeff_(), cloudOldTime.TCoeff(), "T");
    this->relax(VLTrans_(), cloudOldTime.VLTrans(), "T");
}

template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(QTrans_(), "T");
    this->scale(TCoeff_(), "T");
    this->scale(VLTrans_(), "T");
}

template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::preEvolve(td);
}


template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        //Calculate current interface normal field
        n_hat_ = interface_.nearInterface()*fvc::grad(alpha_L_) / ( mag(fvc::grad(alpha_L_)) + interface_.deltaN() );

        typename parcelType::trackingData td(*this);

        this->solve(*this, td);         
    }
}


//- Handle popping of bubbles
template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::HandleBubblePopping()
{
    if (this->solution().canEvolve())
    {
        //Reset popped field to zero
        VolPopped_ = dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero);

        //Remove particles on the interface for too long
        for (parcelType& p : *this)
        {

            //Empirical fit for bubble popping time for test case.
            //In future, could be made more flexible/run time modifiable
            
            const scalar PopTime = 0.1;    // ( 1668.8*( p.d() ) - 0.35795);   5.916 * pow((p.d()/0.00148),0.5); 

            if (p.intTime() > PopTime)
            {

               //Account for volume from popping
               VolPopped_[ p.cell() ] += p.volume();
               
          
               CloudType::deleteParticle(p);
               
            }
        
        }
        
        
        //If bubble fraction in a cell exceeds a limit, then pop all bubbles in that cell:
         const scalar alphac_thresh = 0.6;  
        
        //List of bubbles in each cell
        List<DynamicList<parcelType*>>& BubblesInCells = this->cellOccupancy();
        
        //Continuous-phase volume fraction
        const scalarField alphac = 1.0 - this->theta()().internalField();
        
        //scan over all cells j
        for (label j = 0; j < BubblesInCells.size(); j++ )
        {
            if (alphac[j] < alphac_thresh)
            {
                DynamicList<parcelType*>& BubblesToPop = BubblesInCells[j];
                for (parcelType* ptrBubble : BubblesToPop)
                {
                   //Account for volume from popping
                   VolPopped_[j] += ptrBubble->volume();
                                 
                   CloudType::deleteParticle(*ptrBubble);
                } 
                
                Info<<"Popped" << endl;
                
                //CHECK - Does delete particle already handle this?
                BubblesToPop.clear();
            }
        
        
        }
        
        
    }
}


template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ThermalBubbleCloud<CloudType>::info()
{
    CloudType::info();

    //Todo - print out some relevant info
}


// ************************************************************************* //
