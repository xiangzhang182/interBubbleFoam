/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ThermalBubbleCloud<CloudType>&
Foam::ThermalBubbleCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}

template<class CloudType>
inline const Foam::volScalarField& Foam::ThermalBubbleCloud<CloudType>::alpha_L() const
{
    return alpha_L_;
}

template<class CloudType>
inline const Foam::volScalarField& Foam::ThermalBubbleCloud<CloudType>::T() const
{
    return T_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ThermalBubbleCloud<CloudType>::VolPopped() const
{
    return VolPopped_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::ThermalBubbleCloud<CloudType>::n_hat() const
{
    return n_hat_;
}

template<class CloudType>
inline const Foam::twoPhaseMixtureEThermo& Foam::ThermalBubbleCloud<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
inline const Foam::interfaceProperties& Foam::ThermalBubbleCloud<CloudType>::interface() const
{
    return interface_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermalBubbleCloud<CloudType>::QTrans()
{
    return *QTrans_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermalBubbleCloud<CloudType>::QTrans() const
{
    return *QTrans_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermalBubbleCloud<CloudType>::TCoeff()
{
    return *TCoeff_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermalBubbleCloud<CloudType>::TCoeff() const
{
    return *TCoeff_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermalBubbleCloud<CloudType>::VLTrans()
{
    return *VLTrans_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermalBubbleCloud<CloudType>::VLTrans() const
{
    return *VLTrans_;
}



//Rattner - convert to a temperature basis
template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermalBubbleCloud<CloudType>::TSource( ) const
{
/*
    DebugInfo
        << "hsTrans min/max = " << min(hsTrans()).value() << ", "
        << max(hsTrans()).value() << nl
        << "hsCoeff min/max = " << min(hsCoeff()).value() << ", "
        << max(hsCoeff()).value() << endl;
*/
    if (this->solution().coupled())
    {
        const volScalarField::Internal dV( this->mesh().V() );

        if (this->solution().semiImplicit("T"))
        {
            return
                QTrans()/dV
              - fvm::SuSp(TCoeff()/dV, T_)
              + TCoeff()*T_/dV;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(T_, dimEnergy/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            //fvm.source() = -QTrans()/dV;
            fvm.source() = -QTrans();

            return tfvm;
        }
    }

    //If the solution is not coupled, return an empty fvScalar field
    return tmp<fvScalarMatrix>::New(T_, dimEnergy/dimTime/dimVolume);
}


//Rattner - bubbles to liquid volume source
template<class CloudType>
inline Foam::tmp<Foam::volScalarField>  // Need to doublecheck here
Foam::ThermalBubbleCloud<CloudType>::VLSource( ) const
{

    //If solution is coupled, return a tmp vol scalar field for volume source, scaled by dt and cell size

    //If not coupled, return 0 tmp Vol Scalar field

    tmp<volScalarField> VLSource_
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":VLSource",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimTime, Zero)
        )
    );

    if (this->solution().coupled())
    {
        scalarField& VLSource_prim = VLSource_.ref().primitiveFieldRef();

        const scalarField& dV = this->mesh().V() ;
        const scalar dt = this->db().time().deltaTValue();
        
        VLSource_prim = VLTrans().field() /(dV*dt) ;
    }

    //If the solution is not coupled, return an empty fvScalar field
    return VLSource_;
}


// ************************************************************************* //
