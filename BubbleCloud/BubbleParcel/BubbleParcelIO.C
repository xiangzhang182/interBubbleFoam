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

#include "BubbleParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::BubbleParcel<ParcelType>::propertyList_ =
    Foam::BubbleParcel<ParcelType>::propertyList();


template<class ParcelType>
const std::size_t Foam::BubbleParcel<ParcelType>::sizeofFields
(
    sizeof(BubbleParcel<ParcelType>) - sizeof(ParcelType)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::BubbleParcel<ParcelType>::BubbleParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    intTime_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> intTime_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &intTime_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&intTime_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::BubbleParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> intTime(c.fieldIOobject("intTime", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, intTime);

    label i = 0;
    for (BubbleParcel<ParcelType>& p : c)
    {
        p.intTime_ = intTime[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::BubbleParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool valid = np;

    IOField<scalar>
        intTime(c.fieldIOobject("intTime", IOobject::NO_READ), np);

    label i = 0;

    for (const BubbleParcel<ParcelType>& p : c)
    {
        intTime[i] = p.intTime_;

        ++i;
    }

    intTime.write(valid);
}


template<class ParcelType>
void Foam::BubbleParcel<ParcelType>::writeProperties
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

    writeProp("intTime", intTime_);

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::BubbleParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readFields(c);

    if (!c.size()) return;

    auto& intTime = cloud::lookupIOField<scalar>("intTime", obr);

    label i = 0;
    for (BubbleParcel<ParcelType>& p : c)
    {
        p.intTime_ = intTime[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::BubbleParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    auto& intTime = cloud::createIOField<scalar>("intTime", np, obr);

    label i = 0;
    for (const BubbleParcel<ParcelType>& p : c)
    {
        intTime[i] = p.intTime_;

        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const BubbleParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.intTime();
    }
    else
    {
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
