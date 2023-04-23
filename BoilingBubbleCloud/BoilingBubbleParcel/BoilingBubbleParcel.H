/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::BubbleParcel

Group
    grpLagrangianIntermediateParcels

Description
    Wrapper around bubble parcels with phase change HT added in

SourceFiles
    BoilingBubbleParcelI.H
    BoilingBubbleParcelTrackingDataI.H
    BoilingBubbleParcel.C
    BoilingBubbleParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BoilingBubbleParcel_H
#define BoilingBubbleParcel_H

#include "particle.H"
#include "BubbleParcel.H"        
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"

#include "twoPhaseMixtureEThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class ParcelType>
class BoilingBubbleParcel;


// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const BoilingBubbleParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class BoilingBubbleParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class BoilingBubbleParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


    //- Need to add Boiling properties to the tracking data for bubbles
    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            //- Local copy of carrier specific heat field
            //  Cp not stored on carrier thermo, but returned as tmp<...>
            const volScalarField Cp_;

            //- Local copy of carrier Thermal conductivity field
            //  kappa not stored on carrier thermo, but returned as tmp<...>
            const volScalarField kappa_;

			//- Same story for sigma (surface tension
			const volScalarField sigma_;

            // Interpolators for continuous phase fields

                //- Temperature field interpolator
                autoPtr<interpolation<scalar>> TInterp_;
                
                //- Specific heat capacity field interpolator
                autoPtr<interpolation<scalar>> CpInterp_;

                //- Thermal conductivity field interpolator
                autoPtr<interpolation<scalar>> kappaInterp_;
				
				autoPtr<interpolation<scalar>> sigmaInterp_;
                

            // Cached continuous phase properties

                //- Liquid phase volume fraction
                scalar Tc_;
                
                //- Specific heat capacity [J/(kg.K)]
                scalar Cpc_;
                

    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member functions
        
        
            //- Return access to the locally stored carrier Cp field
            inline const volScalarField& Cp() const;

            //- Return access to the locally stored carrier kappa field
            inline const volScalarField& kappa() const;
			inline const volScalarField& sigma() const;

                
            //- Return const access to the interpolators
            inline const interpolation<scalar>& TInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase specific heat capacity field
            inline const interpolation<scalar>& CpInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase thermal conductivity field
            inline const interpolation<scalar>& kappaInterp() const;
            inline const interpolation<scalar>& sigmaInterp() const;


            //- Return the continuous phase props
            inline scalar Tc() const;

            //- Access the continuous phase props
            inline scalar& Tc();
            
            //- Return the continuous phase specific heat capacity
            inline scalar Cpc() const;

            //- Access the continuous phase specific heat capacity
            inline scalar& Cpc();
            
    };


protected:

    // Protected data
    //ASR - no additional state data for now     
    
    
    //- Calculate new particle diameter (vs. new particle Temperature from regular ThermoParcel)
    template<class TrackCloudType>
    scalar calcHeatTransfer
    (
        TrackCloudType& cloud,
        trackingData& td,
        const scalar dt,           // timestep
        const scalar Re,           // Reynolds number
        const scalar Pr,           // Prandtl number - surface
        const scalar kappa,        // Thermal conductivity - surface
        scalar& dQTrans,           // Heat transfer to carrier
        scalar& SpT,               // linearised heat transfer coefficient
        scalar& dVLTrans           // Volume source in continuous liquid from bubbles condensing/evaporating
    );
    
    

public:

    // Static data members

        //- Runtime type information
        TypeName("BoilingBubbleParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            ""               //  Empty here
        );

    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline BoilingBubbleParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline BoilingBubbleParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline BoilingBubbleParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const scalar& intTime0,        
            const typename ParcelType::constantProperties& constProps
        );

        //- Construct from Istream
        BoilingBubbleParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        BoilingBubbleParcel(const BoilingBubbleParcel& p);

        //- Construct as a copy
        BoilingBubbleParcel(const BoilingBubbleParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new BoilingBubbleParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new BoilingBubbleParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<BoilingBubbleParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<BoilingBubbleParcel<ParcelType>>
                (
                    new BoilingBubbleParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access



        // Main calculation loop 
        
            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);



          //- Correct cell values using latest transfer information
          //ASR - may need to bring item below back!
          /*
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );
           */


            //- Update parcel properties over the time interval 
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


       // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly  // = false
            ) const;

            //- Read particle fields as objects from the obr registry
            template<class CloudType>
            static void readObjects(CloudType& c, const objectRegistry& obr);

            //- Write particle fields as objects into the obr registry
            template<class CloudType>
            static void writeObjects(const CloudType& c, objectRegistry& obr);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const BoilingBubbleParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BoilingBubbleParcelI.H"
#include "BoilingBubbleParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "BoilingBubbleParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
