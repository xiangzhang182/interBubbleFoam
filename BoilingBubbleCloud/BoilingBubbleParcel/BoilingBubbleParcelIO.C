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

\*---------------------------------------------------------------------------*/

#include "BoilingBubbleParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::BoilingBubbleParcel<ParcelType>::propertyList_ =
    Foam::BoilingBubbleParcel<ParcelType>::propertyList();


template<class ParcelType>
const std::size_t Foam::BoilingBubbleParcel<ParcelType>::sizeofFields
(
    sizeof(BoilingBubbleParcel<ParcelType>) - sizeof(ParcelType)     // Note - Different from ThermoParcelIO.C
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::BoilingBubbleParcel<ParcelType>::BoilingBubbleParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat)
    //  T_(0.0),
    //  Cp_(0.0) 
    
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            //is >> intTime_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            //readRawScalar(is, &intTime_);

            is.endRawRead();
        }
        else
        {
            //is.read(reinterpret_cast<char*>(&intTime_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::BoilingBubbleParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    //IOField<scalar> intTime(c.fieldIOobject("intTime", IOobject::MUST_READ), valid);
    //c.checkFieldIOobject(c, intTime);

    label i = 0;
    for (BoilingBubbleParcel<ParcelType>& p : c)
    {
        //p.intTime_ = intTime[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::BoilingBubbleParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool valid = np;

    IOField<vector> NucleationSite_Position_IO(c.fieldIOobject("NucleationSite_Position", IOobject::NO_READ), np);
    IOField<vector> NucleationSite_Normal_IO(c.fieldIOobject("NucleationSite_Normal", IOobject::NO_READ), np);
    IOField<scalar> NucleationSite_Radius_IO(c.fieldIOobject("NucleationSite_Radius", IOobject::NO_READ), np);
    IOField<bool> IsPinned_IO(c.fieldIOobject("IsPinned", IOobject::NO_READ), np);

    label i = 0;

    for (const BoilingBubbleParcel<ParcelType>& p : c)
    {
        NucleationSite_Position_IO[i] = p.NucleationSite_Position;
        NucleationSite_Normal_IO[i] = p.NucleationSite_Normal;
        NucleationSite_Radius_IO[i] = p.NucleationSite_Radius;
        IsPinned_IO[i] = p.IsPinned;    
//        intTime[i] = p.intTime_;

        ++i;
    }

    NucleationSite_Position_IO.write(valid);
    NucleationSite_Normal_IO.write(valid);
    NucleationSite_Radius_IO.write(valid);
    IsPinned_IO.write(valid);

//    intTime.write(valid);
}


template<class ParcelType>
void Foam::BoilingBubbleParcel<ParcelType>::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    ParcelType::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        ParcelType::writeProperty(os, Name, Value, namesOnly, delim, filters)

        writeProp("NucleationSite_Position", NucleationSite_Position);
        writeProp("NucleationSite_Normal", NucleationSite_Normal);
        writeProp("NucleationSite_Radius", NucleationSite_Radius);
        writeProp("IsPinned", IsPinned);
        
//    writeProp("intTime", intTime_);

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::BoilingBubbleParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readFields(c);

    if (!c.size()) return;
    
//    auto& intTime = cloud::lookupIOField<scalar>("intTime", obr);

    label i = 0;
    for (BoilingBubbleParcel<ParcelType>& p : c)
    {
//        p.intTime_ = intTime[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::BoilingBubbleParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    auto& NucleationSite_Position_IO = cloud::createIOField<vector>("NucleationSite_Position", np, obr);
    auto& NucleationSite_Normal_IO = cloud::createIOField<vector>("NucleationSite_Normal", np, obr);
    auto& NucleationSite_Radius_IO = cloud::createIOField<scalar>("NucleationSite_Radius", np, obr);
    auto& IsPinned_IO = cloud::createIOField<bool>("IsPinned", np, obr);
 //    auto& intTime = cloud::createIOField<scalar>("intTime", np, obr);

    label i = 0;
    for (const BoilingBubbleParcel<ParcelType>& p : c)
    {
        NucleationSite_Position_IO[i] = p.NucleationSite_Position;
        NucleationSite_Normal_IO[i] = p.NucleationSite_Normal;
        NucleationSite_Radius_IO[i] = p.NucleationSite_Radius;
        IsPinned_IO[i] = p.IsPinned;    

//        intTime[i] = p.intTime_;

        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const BoilingBubbleParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.NucleationSite_Position
            << token::SPACE << p.NucleationSite_Normal
            << token::SPACE << p.NucleationSite_Radius
            << token::SPACE << p.IsPinned;

 //           << token::SPACE << p.intTime();
    }
    else
    {
        //Rattner - note this binary write mode is not updated for the boiling bubble parcel
       os  << static_cast<const ParcelType&>(p);
       os.write
       (
           reinterpret_cast<const char*>(&p.intTime_),
           BubbleParcel<ParcelType>::sizeofFields
       );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
