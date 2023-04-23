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

#include "BoilingBubbleCloud.H"

#include "distributionModel.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::setModels()
{


}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoilingBubbleCloud<CloudType>::BoilingBubbleCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    const twoPhaseMixtureEThermo& thermo,
    const interfaceProperties& interface,
    bool  readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
    thermo_(thermo),
    cloudCopyPtr_(nullptr),
    interface_(interface),
	transportProperties_( rho.db().objectRegistry::lookupObject<IOdictionary>("transportProperties") ),
	sigmaPtr_(surfaceTensionModel::New( transportProperties_, rho.mesh())),
    alpha_L_(thermo.alpha1()),    
    T_(thermo.T()),
    VolPopped_
    (
        IOobject
        (
            "VolPopped",
            rho.mesh().time().timeName(),
            rho.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho.mesh(),
        dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero)
    ),
    
    //To calculate pinning forces at the interface that keep the bubble from escaping, need access to interface normal field
    n_hat_
    (
        IOobject
        (
            "n_hat",
            rho.mesh().time().timeName(),
            rho.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho.mesh(),
        dimensionedVector("tmp",dimless,vector(0,0,0))
    ),
    QTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":QTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime, Zero)
        )
    ),
    TCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":TCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime/dimTemperature, Zero)
        )
    ),
    VLTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":VLTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimVolume, Zero)
        )
    ),
	NucleateBoilingProperties
	(
	    this->particleProperties().subOrEmptyDict
        (
            "NucleateBoiling",
            keyType::REGEX,
            this->solution().active()
        )
	)

{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
		
		InitNucleateBoiling();
    }
    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::BoilingBubbleCloud<CloudType>::BoilingBubbleCloud
(
    BoilingBubbleCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermo_(c.thermo()),
    cloudCopyPtr_(nullptr),
    interface_(c.interface()),
	transportProperties_( c.alpha_L().db().objectRegistry::lookupObject<IOdictionary>("transportProperties") ),
	sigmaPtr_(surfaceTensionModel::New( transportProperties_, c.alpha_L().mesh())),
    alpha_L_(c.alpha_L()),
    T_(c.T()),
    VolPopped_
    (
        IOobject
        (
            "VolPopped",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero)
    ),
    n_hat_
    (
        IOobject
        (
            "n_hat",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedVector("tmp",dimless,vector(0,0,0))
    ),
    QTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":QTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime/dimVolume, Zero)
        )
    ),
    TCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":TCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime/dimVolume/dimTemperature, Zero)
        )
    ),
    VLTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":VLTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimVolume, Zero)
        )
    ),
	NucleateBoilingProperties
	(
	    this->particleProperties().subOrEmptyDict
        (
            "NucleateBoiling",
            keyType::REGEX,
            this->solution().active()
        )
	)
    
{}


template<class CloudType>
Foam::BoilingBubbleCloud<CloudType>::BoilingBubbleCloud
(
    const fvMesh& mesh,
    const word& name,
    const BoilingBubbleCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    thermo_(c.thermo()),
    cloudCopyPtr_(nullptr),
    interface_(c.interface()),
	transportProperties_( mesh.thisDb().objectRegistry::lookupObject<IOdictionary>("transportProperties") ),
	sigmaPtr_(surfaceTensionModel::New( transportProperties_, mesh) ),
    alpha_L_(c.alpha_L()),
    T_(c.T()),
    VolPopped_
    (
        IOobject
        (
            "VolPopped",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero)
    ),
    n_hat_
    (
        IOobject
        (
            "n_hat",
            alpha_L_.mesh().time().timeName(),
            alpha_L_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_L_.mesh(),
        dimensionedVector("tmp",dimless,vector(0,0,0))
    ),
    QTrans_(nullptr),
    TCoeff_(nullptr),
    VLTrans_(nullptr),
	NucleateBoilingProperties
	(
	    this->particleProperties().subOrEmptyDict
        (
            "NucleateBoiling",
            keyType::REGEX,
            this->solution().active()
        )
	)
       
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoilingBubbleCloud<CloudType>::~BoilingBubbleCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<BoilingBubbleCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::restoreState()
{
    this->cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    QTrans_->field() = 0.0;
    TCoeff_->field() = 0.0;
    VLTrans_->field() = 0.0;
}

template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::relaxSources
(
    const BoilingBubbleCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    this->relax(QTrans_(), cloudOldTime.QTrans(), "T");
    this->relax(TCoeff_(), cloudOldTime.TCoeff(), "T");
    this->relax(VLTrans_(), cloudOldTime.VLTrans(), "T");
}

template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(QTrans_(), "T");
    this->scale(TCoeff_(), "T");
    this->scale(VLTrans_(), "T");
}

template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::preEvolve(td);
	
	//Adds nucleating bubbles to any unoccupied sites
	if ( NucleateBoilingProperties.lookupOrDefault("active", false) )
	{	PopulateNucleationSites(); }
}


template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        //Calculate current interface normal field
        n_hat_ = interface_.nearInterface()*fvc::grad(alpha_L_) / ( mag(fvc::grad(alpha_L_)) + interface_.deltaN() );

        typename parcelType::trackingData td(*this);

        this->solve(*this, td);         
    }
}


//- Handle popping of bubbles
template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::HandleBubblePopping()
{
    if (this->solution().canEvolve())
    {
        //Reset popped field to zero
        VolPopped_ = dimensionedScalar(dimensionSet(0, 3, 0, 0, 0, 0, 0) , Zero);

        //Remove particles on the interface for too long
        for (parcelType& p : *this)
        {

            //Empirical fit for bubble popping time for test case.
            //In future, could be made more flexible/run time modifiable
            
            const scalar PopTime = 0.1;    // ( 1668.8*( p.d() ) - 0.35795);   5.916 * pow((p.d()/0.00148),0.5); 

            if (p.intTime() > PopTime)
            {

               //Account for volume from popping
               VolPopped_[ p.cell() ] += p.volume();
               
          
               CloudType::deleteParticle(p);
               
            }
        
        }
        
        
        //If bubble fraction in a cell exceeds a limit, then pop all bubbles in that cell:
         const scalar alphac_thresh = 0.6;  
        
        //List of bubbles in each cell
        List<DynamicList<parcelType*>>& BubblesInCells = this->cellOccupancy();
        
        //Continuous-phase volume fraction
        const scalarField alphac = 1.0 - this->theta()().internalField();
        
        //scan over all cells j
        for (label j = 0; j < BubblesInCells.size(); j++ )
        {
            if (alphac[j] < alphac_thresh)
            {
                DynamicList<parcelType*>& BubblesToPop = BubblesInCells[j];
                for (parcelType* ptrBubble : BubblesToPop)
                {
                   //Account for volume from popping
                   VolPopped_[j] += ptrBubble->volume();
                                 
                   CloudType::deleteParticle(*ptrBubble);
                } 
                
                Info<<"Popped" << endl;
                
                //CHECK - Does delete particle already handle this?
                BubblesToPop.clear();
            }
        
        
        }
        
        
    }
}


template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::info()
{
    CloudType::info();

    //Todo - print out some relevant info
}



//Nucleate boiling model functionality below
template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::InitNucleateBoiling()
{
	if ( NucleateBoilingProperties.lookupOrDefault("active", false) )
	{
		Info<< "Initializing nucleate boiling model" << endl;
		
		const word& PatchName = NucleateBoilingProperties.getWord("patch");
		NucleateBoilingPatchID = this->mesh().boundaryMesh().findPatchID(PatchName);
		if (NucleateBoilingPatchID < 0)
		{
			FatalErrorInFunction
				<< "Requested patch " << PatchName << " not found" << nl
				<< "Available patches are: " << this->mesh().boundaryMesh().names() << nl
				<< exit(FatalError);
		}
		
		
		//Now - based on dictionary info, randomly distribute nucleation sites:
		//Copying approach from patchinjection
		
		const polyPatch& patch = this->mesh().boundaryMesh()[NucleateBoilingPatchID];
		const pointField& points = patch.points();

		const labelList& cellOwners = patch.faceCells();
		
		// Triangulate the patch faces and create addressing
		DynamicList<label> triToFace(2*patch.size());
		DynamicList<scalar> triMagSf(2*patch.size());
		DynamicList<face> triFace(2*patch.size());
		DynamicList<face> tris(5);

		// Set zero value at the start of the tri area list
		triMagSf.append(0.0);

		forAll(patch, facei)
		{
			const face& f = patch[facei];

			tris.clear();
			f.triangles(points, tris);

			forAll(tris, i)
			{
				triToFace.append(facei);
				triFace.append(tris[i]);
				triMagSf.append(tris[i].mag(points));
			}
		}

		//Build up cumulative area fraction over triangles (from faces)
		for (label i = 1; i < triMagSf.size(); i++)
		{
			triMagSf[i] += triMagSf[i-1];
		}
				
		const scalar patchArea_Proc( triMagSf[triMagSf.size()-1] );
		
		
		const scalarField magSf(mag(patch.faceAreas()));
		const vectorField patchNormals = patch.faceAreas()/magSf;
			
		scalar patchArea_Global = patchArea_Proc;
		reduce(patchArea_Global, sumOp<scalar>());
		Info<< "Total boiling patch area: " << patchArea_Global << endl;
		Info<< "Nucleation sites (global): " << NucleationSiteCount_Global << endl;


		//Now create nucleation sites on area
		const scalar SiteDensity = NucleateBoilingProperties.lookupOrDefault("SiteDensity", 0.0);
		if (SiteDensity == 0)
		{
			FatalErrorInFunction
				<< "Undefined nucleate boiling site density" << nl
				<< exit(FatalError);
		}
		NucleationSiteCount_Proc = ceil( SiteDensity * patchArea_Proc );
		NucleationSiteCount_Global = NucleationSiteCount_Proc;
		reduce(NucleationSiteCount_Global, sumOp<label>());
		
		
		
		//Initialize nucleation site lists per processor with correct number of nucleation sites
		NucleationSite_Positions.resize(NucleationSiteCount_Proc, Zero);
		NucleationSite_Radii.resize(NucleationSiteCount_Proc, Zero);
		NucleationSite_Normals.resize(NucleationSiteCount_Proc, Zero);
		NucleationSite_BubbleIDs.resize(NucleationSiteCount_Proc, -1);
		NucleationSite_CellOwners.resize(NucleationSiteCount_Proc,-1);
		
		//Now, create nucleation sites, 1 by 1
		Random rnd;
		
		//Object for generating site size:
		const autoPtr<distributionModel> SizeDistribution =
			distributionModel::New
			(
				NucleateBoilingProperties.subDict("RadiusDistribution"),
				rnd
			);
		
		for (label i = 0; i < NucleationSiteCount_Proc; i++)
		{
			//Random number to generate face
			const scalar f = rnd.sample01<scalar>() * patchArea_Proc;
			
			//Find containing triangle
			label trii = 0;
			for (label j = 1; j < triMagSf.size(); j++)
			{
				if( triMagSf[j] >= f )
				{
					trii = j-1;
					break;
				}	
			}
			
			//Info<< "Site face number: " << triToFace[trii] << endl;
			
			//Now get a random point on that tri/face
			const face& tf = triFace[trii];
            const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);
			point p = tri.randomPoint(rnd);
			// Apply corrections to position for 2-D cases
            meshTools::constrainToMeshCentre(this->mesh(), p);
            NucleationSite_Positions[i] = p;
			//Info<< "Nucleation site: " << NucleationSite_Positions[i] << endl;
			
			//Indicate normal direction (into domain)
			NucleationSite_Normals[i] = -patchNormals[triToFace[trii]];
			//Info<< "Nucleation site normal: " << NucleationSite_Normals[i] << endl;
			
			//Set cells that own the site
			NucleationSite_CellOwners[i] = cellOwners[triToFace[trii]];
			
			//Set nucleation site size
			NucleationSite_Radii[i] = SizeDistribution->sample();
			//Info<< "Nucleation site size: " << NucleationSite_Radii[i] << endl;
		
			
		
		
		}		
		
		
		
		
		
		
	}
		
	
}


//Create nucleating bubbles for any unoccupied sites
template<class CloudType>
void Foam::BoilingBubbleCloud<CloudType>::PopulateNucleationSites()
{
	
	for (label i = 0; i < NucleationSiteCount_Proc; i++) //Iterate over all sites on processor
	{
		if ( NucleationSite_BubbleIDs[i] != -1)
		{ continue; } //Site already has a nucleating bubble in it

		//Properties for the new bubble
		label celli = NucleationSite_CellOwners[i];
		vector pos = NucleationSite_Positions[i] + NucleationSite_Radii[i]*NucleationSite_Normals[i];
		
		// Apply corrections to position for 2-D cases
        meshTools::constrainToMeshCentre(this->mesh(), pos);
		
		// Create a new parcel
        parcelType* pPtr = new parcelType(this->mesh(), pos, celli);

        // Check/set new parcel properties
        this->setParcelThermoProperties(*pPtr, 0.0);
		pPtr->U() = vector(Zero);
        pPtr->d() = 2.0*NucleationSite_Radii[i];
		pPtr->nParticle() = 1;

        // Check/set new parcel injection properties
        this->checkParcelProperties(*pPtr, 0, false);

		//Not sure which of the following are needed:
		//parcelsAdded++;
        //massAdded += pPtr->nParticle()*pPtr->mass();
		//pPtr->move(cloud, td, dt)
        //pPtr->typeId() = injectorID_;
        
		//Info about the nucleation site for this bubble
		pPtr->NucleationSite_Position =  NucleationSite_Positions[i];
		pPtr->NucleationSite_Normal = NucleationSite_Normals[i];
		pPtr->NucleationSite_Radius = NucleationSite_Radii[i];
		pPtr->IsPinned = true;
		
		this->addParticle(pPtr);
		NucleationSite_BubbleIDs[i] = pPtr->origId();
		Info<<"Added bubble ID: " << pPtr->origId() << endl;

	}
}

// ************************************************************************* //
