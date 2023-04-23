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

#include "BoilingBubbleParcel.H"


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::BoilingBubbleParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    tetIndices tetIs = this->currentTetIndices();

    td.Cpc() = td.CpInterp().interpolate(this->coordinates(), tetIs);

    td.Tc() = td.TInterp().interpolate(this->coordinates(), tetIs);
}


//ASR - may need to add cellValueSourceCorrection around here




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::BoilingBubbleParcel<ParcelType>::BoilingBubbleParcel
(
    const BoilingBubbleParcel<ParcelType>& p
)
:
    ParcelType(p)
    // T_(p.T_),
    // Cp_(p.Cp_)
    
{}


template<class ParcelType>
Foam::BoilingBubbleParcel<ParcelType>::BoilingBubbleParcel
(
    const BoilingBubbleParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ASR - content copied from ThermoParcel
template<class ParcelType>
template<class TrackCloudType>
void Foam::BoilingBubbleParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{

    //Call the parent calculation function
    ParcelType::calc(cloud, td, dt);
    

    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar mass0 = this->mass();

    // Store T for consistent radiation source
    // const scalar T0 = this->T_; 

//Info<< "Cp = " << td.Cpc() << endl;
//Info<< "mu = " << td.muc() << endl;

    
    // Reynolds number
    const scalar Re = this->Re(td);
    
    // Prandtl number
    
    const tetIndices tetIs = this->currentTetIndices();
    const scalar kappa = td.kappaInterp().interpolate(this->coordinates(), tetIs);
    
    const scalar Pr = td.Cpc()*td.muc()/kappa;
    

    // Sources
    // ~~~~~~~

    // Explicit enthalpy source for particle
    // scalar Sh = 0.0;

    // Linearised heat transfer source coefficient
    scalar SpT = 0.0;

    // Heat transfer from the particle to the carrier phase
    scalar dQTrans = 0.0;

    // Continuous-phase volume source due to bubble evaporation or condensation
    scalar dVLTrans = 0.0;


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle size
    this->d_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Re,
            Pr,
            kappa,
            dQTrans,
            SpT,
            dVLTrans
        );

    //ASR TODO - for now, just transfer heat to carrier phase, later add mass transfer as well
    
    // Remove the particle when mass falls below minimum threshold
    if ( this->d() == 0.0 )
    {
        td.keepParticle = false;

        if (cloud.solution().coupled())
        {

            /*            
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(Y_, i)
            {
                scalar dmi = dm*Y_[i];
                label gid = composition.localToCarrierId(0, i);
                scalar hs = composition.carrier().Hs(gid, td.pc(), T0);

                cloud.rhoTrans(gid)[this->cell()] += dmi;
                cloud.hsTrans()[this->cell()] += dmi*hs;
            }
            cloud.UTrans()[this->cell()] += dm*U0;

            cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
            */
        }
    }

     //Info<< "D Particle = " << this->d() << endl;
      

    //  Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (cloud.solution().coupled())
    {
        // Update heat transfer
        cloud.QTrans()[this->cell()] += np0*dQTrans;

        // Update heat transfer
        cloud.TCoeff()[this->cell()] += np0*SpT;

        // Update the volume transfer between bubbles and continuous phase
        cloud.VLTrans()[this->cell()] += np0*dVLTrans;
    }

}


//ASR - need to update calcHeatTransfer below

//Returns new diameter of parcel

template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::BoilingBubbleParcel<ParcelType>::calcHeatTransfer
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    scalar& dQTrans,
    scalar& SpT,
    scalar& dVLTrans
)
{
    //Rattner - in the future, can add a switch to the cloud dictionary to activate/deactivate heat transfer
    //if (!cloud.heatTransfer().active())
    //{
    //    return this->d_;
    //}

    const scalar d = this->d();
    const scalar rho_V = this->rho();                      //  cloud.thermo().rho2().value() ? 
    const scalar rho_L = cloud.thermo().rho1().value();
    const scalar As = this->areaS(d);
    const scalar V = this->volume(d);
    const scalar m = rho_V*V;
    const scalar h_LV = cloud.thermo().Hf2().value() - cloud.thermo().Hf1().value();    // Latent heat of vaparization 

	const tetIndices tetIs = this->currentTetIndices();
    const scalar sigma = td.sigmaInterp().interpolate(this->coordinates(), tetIs);
Info<<"Bubble sigma = " << sigma << endl;

    // Calc heat transfer coefficient
    // Hard code in Ranz Marshall formula
//Info<<"Re = " << Re << endl;
//Info<<"Pr = " << Pr << endl;
    const scalar htc = (kappa/d)*(2.0 + 0.6 * pow(Re,0.5) * pow(Pr,0.33) );
//Info<<"HTC = " << htc << endl;
    // Assume explicit particle heat transfer for calculating new size
    // Q is positive for evaporation
    scalar Q = dt*As*htc*(td.Tc() - cloud.thermo().TSat().value() );
//Info<<"Q = " << Q << endl;
    //Limit so that bubble isn't over-condensed to negative mass
    Q = max(Q, -m*h_LV);
    
    const scalar dM = Q/h_LV;
    
    //Negative for condensation 
    dVLTrans = -dM/rho_L;       

    const scalar dnew = pow( max( (6*dM)/ (rho_V * constant::mathematical::pi) + pow(d,3.0), 0), 1.0/3.0 );

    
    //Positive for evaporation  
    dQTrans = Q/dt;

    //Define SpT to be implicit
    SpT = htc * As;

    return dnew;
}




// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "BoilingBubbleParcelIO.C"

// ************************************************************************* //
