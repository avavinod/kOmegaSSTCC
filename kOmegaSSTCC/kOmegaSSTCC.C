/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "kOmegaSSTCC.H"
#include "fvOptions.H"
#include "bound.H"
#include <typeinfo>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volVectorField> kOmegaSSTCC<BasicTurbulenceModel>::rotRateMesh() const
{
    volVectorField rotRate
    (
        IOobject
        (
            "rotRate",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedVector( "dummy", dimensionSet(0,0,-1,0,0,0,0), vector::zero)
    );
    // const auto* MRFZones =
    // this->mesh_.cfindObject<IOMRFZoneList>("MRFProperties");

    // if (!MRFZones)
    // {
    //     return rotRate;
    // }
    // const IOMRFZoneList& MRF_(this->mesh_);

    // if (!MRF_.active())
    // {
    //     return rotRate;
    // }

    // for (label j = 0; j < MRF_.UPtrList<MRFZone>::size(); ++j)
    // {
    //     const MRFZone& mrf = MRF_.UPtrList<MRFZone>::operator[](j);
    //     const label& cellZoneID_ = this->mesh_.cellZones().findZoneID(mrf.name());
    //     if (cellZoneID_ == -1)
    //     {
    //         continue;
    //     }
    //     const labelList& cells = this->mesh_.cellZones()[cellZoneID_];
    //     const vector Omega = mrf.Omega();
    //     forAll(cells, i)
    //     {
    //         label celli = cells[i];
    //         rotRate[celli] += Omega ;
    //     }
    // }
    return rotRate;
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTCC<BasicTurbulenceModel>::onebyOmegaD3
(
    const volScalarField& S2,
    const volScalarField::Internal& sqrtOmega2
) const
{
    const volScalarField::Internal& omega_ = this->omega_();
    tmp<volScalarField::Internal>  D2(max(S2, 0.09*2.0*omega_*omega_));
    return scalar(1.0)/(sqrtOmega2 * D2 * sqrt(D2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTCC<BasicTurbulenceModel>::rTilda
(
    const volSymmTensorField& symmGradU,
    const volTensorField& Omega,
    const volTensorField& hodgeDualrotRateMesh,
    const volScalarField::Internal& onebyOmegaD3
) const
{
    tmp<volTensorField> twoOmegaS = 2.0 * (Omega & symmGradU);
    tmp<volSymmTensorField> DDtS
    (
        fvc::DDt(this->phi(), symmGradU)
    );
    tmp<volTensorField> leviCivitaSRotRate
    (
        hodgeDualrotRateMesh & symmGradU
    );
    tmp<volScalarField::Internal> resultingScalarInternal
    (
        twoOmegaS && (DDtS + (leviCivitaSRotRate & T(leviCivitaSRotRate)))
    );

    return  resultingScalarInternal * onebyOmegaD3;
}

// - Return square of strain rate
template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTCC<BasicTurbulenceModel>::fRotation
(
    const volScalarField& rStarByOnePlusrStar,
    const volScalarField::Internal& rTilda
) const
{
    tmp<volScalarField::Internal> resultingScalarInternal
    (
        (scalar(1.0) + cr1_) * scalar(2.0) * rStarByOnePlusrStar
    );
    return resultingScalarInternal * (scalar(1.0) - cr3_ * tanh(cr2_ * rTilda)) - cr1_;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class BasicTurbulenceModel>
kOmegaSSTCC<BasicTurbulenceModel>::kOmegaSSTCC
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        typeName
    ),

    cr1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "cr1",
            this->coeffDict_,
            1
        )
    ),
    cr2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "cr2",
            this->coeffDict_,
            2
        )
    ),
    cr3_
    (
        dimensionedScalar::getOrAddToDict
        (
            "cr3",
            this->coeffDict_,
            1
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}
    // deltaU_("deltaU", dimVelocity, SMALL),

    // ReThetat_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("ReThetat", alphaRhoPhi.group()),
    //         this->runTime_.timeName(),
    //         this->mesh_,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE,
    //         IOobject::REGISTER
    //     ),
    //     this->mesh_
    // ),

    // gammaInt_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("gammaInt", alphaRhoPhi.group()),
    //         this->runTime_.timeName(),
    //         this->mesh_,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE,
    //         IOobject::REGISTER
    //     ),
    //     this->mesh_
    // ),

    // gammaIntEff_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("gammaIntEff", alphaRhoPhi.group()),
    //         this->runTime_.timeName(),
    //         this->mesh_
    //     ),
    //     this->mesh_,
    //     dimensionedScalar(dimless, Zero)
    // )



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTCC<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        cr1_.readIfPresent(this->coeffDict());
        cr2_.readIfPresent(this->coeffDict());
        cr3_.readIfPresent(this->coeffDict());
        // ce2_.readIfPresent(this->coeffDict());
        // sigmaThetat_.readIfPresent(this->coeffDict());
        // cThetat_.readIfPresent(this->coeffDict());
        // this->coeffDict().readIfPresent("lambdaErr", lambdaErr_);
        // this->coeffDict().readIfPresent("maxLambdaIter", maxLambdaIter_);

        return true;
    }

    return false;
}

template<class BasicEddyViscosityModel>
void kOmegaSSTCC<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicEddyViscosityModel::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    // tmp<volTensorField> tgradU = fvc::grad(U);
    // const volScalarField  (this->S2(tgradU()));

    tmp<volTensorField> tgradU = fvc::grad(U);
    tmp<volSymmTensorField> symmGradU = symm(tgradU());
    tmp<volVectorField> rotRateMesh(this->rotRateMesh());
    tmp<volTensorField> hodgeDualrotRateMesh(*rotRateMesh);
    tmp<volTensorField> Omega(skew(tgradU()) + hodgeDualrotRateMesh);
    
    const volScalarField& Omega2(2*magSqr(Omega));
    const volScalarField& S2(2*magSqr(symmGradU));
    const volScalarField& sqrtS2(sqrt(S2));
    const volScalarField& sqrtOmega2(sqrt(Omega2));
    const volScalarField::Internal& onebyOmegaD3(this->onebyOmegaD3(S2, sqrtOmega2));
    const volScalarField& rStarByOnePlusrStar(sqrtS2/(sqrtS2+sqrtOmega2));
    const volScalarField::Internal& rTilda(this->rTilda(symmGradU,Omega,hodgeDualrotRateMesh,onebyOmegaD3));
    const volScalarField::Internal& fRotation(this->fRotation(rStarByOnePlusrStar,rTilda));
    const volScalarField::Internal& fr1(max(min(fRotation, 1.25), 0.0));


    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // - boundary condition changes a cell value
    // - normally this would be triggered through correctBoundaryConditions
    // - which would do
    //      - fvPatchField::evaluate() which calls
    //      - fvPatchField::updateCoeffs()
    // - however any processor boundary conditions already start sending
    //   at initEvaluate so would send over the old value.
    // - avoid this by explicitly calling updateCoeffs early and then
    //   only doing the boundary conditions that rely on initEvaluate
    //   (currently only coupled ones)

    //- 1. Explicitly swap values on coupled boundary conditions
    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();
    // omegaWallFunctions change the cell value! Make sure to push these to
    // coupled neighbours. Note that we want to avoid the re-updateCoeffs
    // of the wallFunctions so make sure to bypass the evaluate on
    // those patches and only do the coupled ones.
    this->omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    ////- 2. Make sure the boundary condition calls updateCoeffs from
    ////     initEvaluate
    ////     (so before any swap is done - requires all coupled bcs to be
    ////      after wall bcs. Unfortunately this conflicts with cyclicACMI)
    //omega_.correctBoundaryConditions();


    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    const volScalarField F1(this->F1(CDkOmega));
    const volScalarField F23(this->F23());

    {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = this->GbyNu(GbyNu0, F23(), S2());

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            alpha()*rho()*gamma*fr1*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*fr1*divU, this->omega_)
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
          + alpha()*rho()*beta*sqr(this->omegaInf_)
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvOptions(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    {
        // Turbulent kinetic energy equation
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, this->k_)
          + fvm::div(alphaRhoPhi, this->k_)
          - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
         ==
            alpha()*rho()*fr1*this->Pk(G)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*fr1*divU, this->k_)
          - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, tgradU()), this->k_)
          + alpha()*rho()*this->betaStar_*this->omegaInf_*this->kInf_
          + this->kSource()
          + fvOptions(alpha, rho, this->k_)
        );

        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(this->k_);
        bound(this->k_, this->kMin_);
    }


    this->correctNut(S2);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
