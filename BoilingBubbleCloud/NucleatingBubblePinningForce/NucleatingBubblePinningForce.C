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

	//Note - overall algorithm needs to be adjusted somewhat so that the bubble isn't bouncing all over the place

    //Change F if the bubble is still pinned
    if ( p.IsPinned  )
    {
        //Maximum retaining force:
		const tetIndices tetIs = p.currentTetIndices();
		const scalar sigma = td.sigmaInterp().interpolate(p.coordinates(), tetIs);
		
		const scalar mag_F_max = sigma*p.d()*constant::mathematical::pi; //Replace with more accurate calculation depending on things like contact angles, etc.
        const vector x_eq = p.NucleationSite_Position + p.NucleationSite_Normal*( p.d()/2.0 );

        const scalar k_spring = mag_F_max / ( 0.2*p.d() );
        vector F_pinning = k_spring*( x_eq - p.position() );
		//Limit the pinning force
		const scalar mag_F_pinning_0 = mag(F_pinning);
		if (mag_F_pinning_0 > mag_F_max)
		{
			F_pinning *= (mag_F_max/mag_F_pinning_0);
		}
		
		F.Su() = F_pinning;
        
		//Note - breaking of pinning is handled at the cloud level, not as part of this force
		//to reduce the need for coupled calculations with the nucleation sites	
    }
    
    return F;
}


// ************************************************************************* //
