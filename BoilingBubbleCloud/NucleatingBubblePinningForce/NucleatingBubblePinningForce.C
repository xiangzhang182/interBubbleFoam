/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "NucleatingBubblePinningForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NucleatingBubblePinningForce<CloudType>::NucleatingBubblePinningForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false) //I think I may need to include more information here to handle the nucleation site, etc.
{}


template<class CloudType>
Foam::NucleatingBubblePinningForce<CloudType>::NucleatingBubblePinningForce
(
    const NucleatingBubblePinningForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::NucleatingBubblePinningForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp F(Zero); //I think this initializes both Su and Sp to 0

    //Change F if the bubble is still pinned
    if ( p.Pinned()  )
    {
        //Maximum retaining force:
        const scalar sigma = owner().interface().
        const scalar F_max = sigma*p.D()*constant::mathematical::pi; //Replace with more accurate calculation depending on things like contact angles, etc.
        const vector x_eq = p.x_NucleationSite() + p.nhat_NucleationSite()*( p.d()/2.0 );

        const scalar k_spring = F_max / ( 0.2*p.d() );
        vector F_pinning = k_spring*( x_eq - p.position() );

        if ( mag(F_pinning) < F_max) //Bubble still pinned
        {
            F.Su() = F_pinning;
        }
        else
        {
            p.UnPin();
        }


    }
    
    
    
    // (AOB:Eq. 34)
    return F;
}


// ************************************************************************* //
