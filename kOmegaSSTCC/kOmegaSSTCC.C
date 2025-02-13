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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTCC<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSST<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut
    (
        S2
    );

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTCC<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
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
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


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


template<class BasicTurbulenceModel>
void kOmegaSSTCC<BasicTurbulenceModel>::correctProductionTerm()
{
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const tmp<volScalarField> tnu = this->nu();
    const volScalarField::Internal& nu = tnu()();
    const volScalarField::Internal& y = this->y_();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Fields derived from the velocity gradient
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal Omega(sqrt(2*magSqr(skew(tgradU()()))));
    const volScalarField::Internal S(sqrt(2*magSqr(symm(tgradU()()))));
    const volScalarField::Internal Us(max(mag(U()), deltaU_));
    const volScalarField::Internal dUsds((U() & (U() & tgradU()()))/sqr(Us));
    tgradU.clear();

    const volScalarField::Internal Fthetat(this->Fthetat(Us, Omega, nu));

    {
        const volScalarField::Internal t(500*nu/sqr(Us));
        const volScalarField::Internal Pthetat
        (
            alpha()*rho()*(cThetat_/t)*(1 - Fthetat)
        );

        // Transition onset momentum-thickness Reynolds number equation
        tmp<fvScalarMatrix> ReThetatEqn
        (
            fvm::ddt(alpha, rho, ReThetat_)
          + fvm::div(alphaRhoPhi, ReThetat_)
          - fvm::laplacian(alpha*rho*DReThetatEff(), ReThetat_)
         ==
            Pthetat*ReThetat0(Us, dUsds, nu) - fvm::Sp(Pthetat, ReThetat_)
          + fvOptions(alpha, rho, ReThetat_)
        );

        ReThetatEqn.ref().relax();
        fvOptions.constrain(ReThetatEqn.ref());
        solve(ReThetatEqn);
        fvOptions.correct(ReThetat_);
        bound(ReThetat_, 0);
    }

    const volScalarField::Internal ReThetac(this->ReThetac());
    const volScalarField::Internal Rev(sqr(y)*S/nu);
    const volScalarField::Internal RT(k()/(nu*omega()));

    {
        const volScalarField::Internal Pgamma
        (
            alpha()*rho()
           *ca1_*Flength(nu)*S*sqrt(gammaInt_()*Fonset(Rev, ReThetac, RT))
        );

        const volScalarField::Internal Fturb(exp(-pow4(0.25*RT)));

        const volScalarField::Internal Egamma
        (
            alpha()*rho()*ca2_*Omega*Fturb*gammaInt_()
        );

        // Intermittency equation
        tmp<fvScalarMatrix> gammaIntEqn
        (
            fvm::ddt(alpha, rho, gammaInt_)
          + fvm::div(alphaRhoPhi, gammaInt_)
          - fvm::laplacian(alpha*rho*DgammaIntEff(), gammaInt_)
        ==
            Pgamma - fvm::Sp(ce1_*Pgamma, gammaInt_)
          + Egamma - fvm::Sp(ce2_*Egamma, gammaInt_)
          + fvOptions(alpha, rho, gammaInt_)
        );

        gammaIntEqn.ref().relax();
        fvOptions.constrain(gammaIntEqn.ref());
        solve(gammaIntEqn);
        fvOptions.correct(gammaInt_);
        bound(gammaInt_, 0);
    }

    const volScalarField::Internal Freattach(exp(-pow4(RT/20.0)));
    const volScalarField::Internal gammaSep
    (
        min(2*max(Rev/(3.235*ReThetac) - 1, scalar(0))*Freattach, scalar(2))
       *Fthetat
    );

    gammaIntEff_ = max(gammaInt_(), gammaSep);
}


template<class BasicTurbulenceModel>
void kOmegaSSTCC<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Correct k and omega
    kOmegaSST<BasicTurbulenceModel>::correct();

    // Correct Production Term
    correctProductionTerm();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
