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
    Foam::BubbleCloud

Group
    grpLagrangianIntermediateClouds

Description
    Adds Bubble popping modelling to kinematic clouds

SourceFiles
    BubbleCloudI.H
    BubbleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef BubbleCloud_H
#define BubbleCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

//Rattner - none to perform here

/*---------------------------------------------------------------------------*\
                       Class BubbleCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BubbleCloud
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
        typedef BubbleCloud<CloudType> BubbleCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<BubbleCloud<CloudType>> cloudCopyPtr_;


    // Private Member Functions

        //- No copy construct
        BubbleCloud(const BubbleCloud&) = delete;

        //- No copy assignment
        void operator=(const BubbleCloud&) = delete;


protected:

    // Protected data

        // References to the cloud sub-models
	//Rattner - none needed

	//In future, this is where some reference fields to carrier fluid need to be specified

        // Initialisation

            //- Set cloud sub-models
            void setModels();


public:

    // Constructors

        //- Construct given carrier gas fields
        BubbleCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        BubbleCloud
        (
            BubbleCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        BubbleCloud
        (
            const fvMesh& mesh,
            const word& name,
            const BubbleCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new BubbleCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new BubbleCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~BubbleCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const BubbleCloud& cloudCopy() const;


        // Cloud evolution functions

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the cloud
            void evolve();



        //- I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BubbleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BubbleCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //