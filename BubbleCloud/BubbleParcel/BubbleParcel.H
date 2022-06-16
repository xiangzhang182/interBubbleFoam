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
    Wrapper around colliding parcel types to add bubble parcels that can pop at the interface

SourceFiles
    BubbleParcelI.H
    BubbleParcelTrackingDataI.H
    BubbleParcel.C
    BubbleParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BubbleParcel_H
#define BubbleParcel_H

#include "particle.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class ParcelType>
class BubbleParcel;


// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const BubbleParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class BubbleParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class BubbleParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


    //- Need to add the liquid phase volume fraction to the tracking data
    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Liquid phase volume fraction field interpolator
                autoPtr<interpolation<scalar>> alpha_LInterp_;


            // Cached continuous phase properties

                //- Liquid phase volume fraction
                scalar alpha_L_;


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
            //- Return const access to the interpolator for 
            inline const interpolation<scalar>& alpha_LInterp() const;

            //- Return the continuous phase liquid volume fraction
            inline scalar alpha_L() const;

            //- Access the continuous phase liquid volume fraction
            inline scalar& alpha_L();

    };


protected:

    // Protected data

        //- Time bubble has been at an interface
	scalar intTime_;
	
public:

    // Static data members

        //- Runtime type information
        TypeName("BubbleParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " intTime_"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline BubbleParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline BubbleParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline BubbleParcel
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
        BubbleParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        BubbleParcel(const BubbleParcel& p);

        //- Construct as a copy
        BubbleParcel(const BubbleParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new BubbleParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new BubbleParcel(*this, mesh));
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

            autoPtr<BubbleParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<BubbleParcel<ParcelType>>
                (
                    new BubbleParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to force
            inline const scalar& intTime() const;


            //- Return access to temperature
            inline scalar& intTime();

        // Tracking
            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
          /*
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );
           */

        

            //- update the parcel
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
                const bool namesOnly
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
            const BubbleParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BubbleParcelI.H"
#include "BubbleParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BubbleParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
