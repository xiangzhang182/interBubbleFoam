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

Class
    Foam::BoilingBubbleCloud

Group
    grpLagrangianIntermediateClouds

Description
    Adds Bubble popping modelling to kinematic clouds

SourceFiles
    BoilingBubbleCloudI.H
    BoilingBubbleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef BoilingBubbleCloud_H
#define BoilingBubbleCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "volFields.H"

#include "interfaceProperties.H"
#include "twoPhaseMixtureEThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

//Rattner - none to perform here

/*---------------------------------------------------------------------------*\
                       Class BoilingBubbleCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BoilingBubbleCloud
:
    public CloudType
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::parcelType parcelType;

        //- Convenience typedef for this cloud type
        typedef BoilingBubbleCloud<CloudType> BoilingBubbleCloudType;

private:

    // Private data

        //- Cloud copy pointer
        autoPtr<BoilingBubbleCloud<CloudType>> cloudCopyPtr_;
        
    // Private Member Functions

        //- No copy construct
        BoilingBubbleCloud(const BoilingBubbleCloud&) = delete;

        //- No copy assignment
        void operator=(const BoilingBubbleCloud&) = delete;


protected:

    // Protected data

    // References to the cloud sub-models
	// Rattner - none needed

        //Reference to mixture and interface classes for VOF model
        const twoPhaseMixtureEThermo& thermo_;
        const interfaceProperties& interface_;

	    //References to the carrier/continuous phase field
	    const volScalarField& alpha_L_;
	    const volScalarField& T_;
	    
        //- Volume of bubbles popped this time step
        volScalarField VolPopped_;

	    //-Interface normal field
	    volVectorField n_hat_;

        // Initialisation

            //- Set cloud sub-models
            void setModels();

        // Sources

            //- Heat transfer [W/m^3]
            autoPtr<volScalarField::Internal> QTrans_;

            //- Coefficient for carrier phase T equation [W/(m^3*K)]
            autoPtr<volScalarField::Internal> TCoeff_;

            //- Volume source in continuous liquid from bubbles evaporating or condensing [m^3]
            autoPtr<volScalarField::Internal> VLTrans_;

        //Nucleation site and boiling properties - should move to a separate class with its own I/O properties later...
		const dictionary NucleateBoilingProperties;
		label NucleateBoilingPatchID = -1;
		label NucleationSiteCount_Global = 0;
		label NucleationSiteCount_Proc = 0;
		vectorList NucleationSite_Positions;
		scalarList NucleationSite_Radii;
		vectorList NucleationSite_Normals;
		labelList NucleationSite_BubbleIDs;
		
		//Protected nucleation methods
		void InitNucleateBoiling();

public:

    // Constructors

        //- Construct given carrier fluid fields
        BoilingBubbleCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            const twoPhaseMixtureEThermo& thermo,
            const interfaceProperties& interface,
            bool  readFields = true
        );

        //- Copy constructor with new name
        BoilingBubbleCloud
        (
            BoilingBubbleCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        BoilingBubbleCloud
        (
            const fvMesh& mesh,
            const word& name,
            const BoilingBubbleCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new BoilingBubbleCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new BoilingBubbleCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~BoilingBubbleCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const BoilingBubbleCloud& cloudCopy() const;

            //- Return const access to the mixture field
            inline const twoPhaseMixtureEThermo& thermo() const;
            inline const interfaceProperties& interface() const;
                        

            //- Return const access to the carrier volume fraction field
            inline const volScalarField& alpha_L() const;
            inline const volScalarField& T() const;

            //- Return const access to the interface normal field fraction field
            inline const volVectorField& n_hat() const;            
            
            //- Return const access to the carrier volume fraction field
            inline const volScalarField& VolPopped() const;
  
      
       // Sources
           
            // Enthalpy
           
            //- heat transfer [W/m^3]
            inline volScalarField::Internal& QTrans();

            //- heat transfer [W/m^3]
            inline const volScalarField::Internal& QTrans() const;

            //- Return coefficient for carrier phase T equation
            inline volScalarField::Internal& TCoeff();

            //- Return const coefficient for carrier phase T equation
            inline const volScalarField::Internal& TCoeff() const;

            //- Return Bubbles -> Liquid source term
            inline volScalarField::Internal& VLTrans();

            //- Return Bubbles -> Liquid source term
            inline const volScalarField::Internal& VLTrans() const;

            //- Return heat transfer source term [W/m^3]
            //- Note - modified to not require passing in enthalpy/temperature field
            inline tmp<fvScalarMatrix> TSource() const;

            //- Returns volume source term for pressure equation from bubbles condensing or evaporating
            inline tmp<volScalarField> VLSource() const;



        // Cloud evolution functions

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const BoilingBubbleCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve(const typename parcelType::trackingData& td);

            //- Evolve the cloud
            void evolve();
            
            //- Handle popping of bubbles
            void HandleBubblePopping();

        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        //- I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BoilingBubbleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BoilingBubbleCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
