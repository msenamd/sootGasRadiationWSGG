/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::radiation::sootYield

Description
    This soot model is purely an state model. The ammount of soot produced is
    determined as a Triangle function of the mixture fraction (Burke–Schumann like)
    with maximum mass fraction etaSoot*Z_st

    etaSoot is prescribed by the user.

    The input dictionary reads like in the radiationProperties dictionary:

    sootModel sootYield<gasHThermoPhysics>;

    sootYieldCoeffs
    {
        fuel          	C3H8;       
    	fvMax         	0.5;
    }

SourceFiles
    sootYield.C

\*---------------------------------------------------------------------------*/

#ifndef sootYield_H
#define sootYield_H

#include "interpolationLookUpTable.H"
#include "sootModel.H"
#include "HashTable.H"

#include "fluidThermo.H"
#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class sootYield Declaration
\*---------------------------------------------------------------------------*/
template<class ThermoType>
class sootYield
:
    public sootModel
{

    // Static functions

        //- Check mixture in thermo
        static const reactingMixture<ThermoType>& checkThermo
        (
            const fluidThermo&
        );


    // Private data

        //- Soot volume fraction
        volScalarField fv_;

        //- Soot mass fraction
        volScalarField Ysoot_;

        //- Soot model dictionary
        dictionary coeffsDict_;

        //- Solution flag
        Switch solveSoot_;

        //- Soot volume fraction max.
        scalar sootYield_;

        //- Thermo package
        const fluidThermo& thermo_;

        const psiReactionThermo& thermo;


public:

    //- Runtime type information
    TypeName("sootYield");


    // Constructors

        //- Construct from components
        sootYield
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Destructor
    virtual ~sootYield();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();


        // Access

            //- Return Xsoot
            const volScalarField& soot() const
            {
                return fv_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sootYield.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
