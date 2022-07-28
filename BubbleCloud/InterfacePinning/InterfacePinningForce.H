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

Class
    Foam::InterfacePinningForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Calculates particle InterfacePinning force

SourceFiles
    InterfacePinningForceI.H
    InterfacePinningForce.C

\*---------------------------------------------------------------------------*/

#ifndef InterfacePinningForce_H
#define InterfacePinningForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class InterfacePinningForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InterfacePinningForce
:
    public ParticleForce<CloudType>
{
    // Private data

        //- Reference to the acceleration due to gravity
        const vector& g_;


public:

    //- Runtime type information
    TypeName("InterfacePinning");


    // Constructors

        //- Construct from mesh
        InterfacePinningForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        InterfacePinningForce(const InterfacePinningForce& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new InterfacePinningForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~InterfacePinningForce();


    // Member Functions

        // Access

            //- Return the acceleration due to gravity
            inline const vector& g() const;


        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InterfacePinningForceI.H"

#ifdef NoRepository
    #include "InterfacePinningForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
