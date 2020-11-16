/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "FixedList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ReactionRateSH::ReactionRateSH
(
    const speciesTable& st,
    const dictionary& dict
)
:
    A_(readScalar(dict.lookup("A"))),
    beta_(readScalar(dict.lookup("beta"))),
    Ename_(dict.lookupOrDefault<word>("Ename", "E")),
    TEname_(dict.lookupOrDefault<word>("TEname", "TE"))
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::ReactionRateSH::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    //const volScalarField& FieldE = mesh_.lookupObject<volScalarField>(Ename_);
    //const volScalarField& FieldTE = mesh_.lookupObject<volScalarField>(TEname_);

    //scalar TE = FieldTE[0];//FieldTE[i]//10000.0;
    //scalar E = FieldE[0];//FieldE[i]//0.1;

    scalar TE = 10000.0;
    scalar E = 0.1;

    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(TE) > vSmall)
    {
        ak *= exp(-TE/T);
    }

    ak *= E;

    return ak;
}


inline Foam::scalar Foam::ReactionRateSH::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    //const volScalarField& FieldE = mesh_.lookupObject<volScalarField>(Ename_);
    //const volScalarField& FieldTE = mesh_.lookupObject<volScalarField>(TEname_);

    //scalar TE = FieldTE[0];//FieldTE[i]//10000.0;
    //scalar E = FieldE[0];//FieldE[i]//0.1;

    scalar TE = 10000.0;
    scalar E = 0.1;

    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(TE) > vSmall)
    {
        ak *= exp(-TE/T);
    }

    ak *= E;

    return ak*(beta_+TE/T)/T;
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::ReactionRateSH::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline void Foam::ReactionRateSH::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::ReactionRateSH::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0;
}


inline void Foam::ReactionRateSH::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactionRateSH& lhrr
)
{
    lhrr.write(os);
    return os;
}


// ************************************************************************* //
