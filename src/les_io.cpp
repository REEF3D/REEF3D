/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"les_io.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

les_io::les_io(lexer *p, fdm *a) : strain(p,a), uprime(p), vprime(p), wprime(p)
{
}

les_io::~les_io()
{
}

void les_io::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{

}

double les_io::kinval(int ii, int jj, int kk)
{
    double val=0.0;

    return val;
}

double les_io::epsval(int ii, int jj, int kk)
{
    double val=0.0;

    return val;
}

double les_io::ccipol_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double les_io::ccipol_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double les_io::ccipol_a_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double les_io::ccipol_a_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

void les_io::gcupdate(lexer *p, fdm *a, ghostcell *pgc)
{

}

void les_io::kinget(int ii, int jj, int kk,double val)
{

}

void les_io::epsget(int ii, int jj, int kk,double val)
{

}

void les_io::ini(lexer* p, fdm*a, ghostcell* pgc)
{
}

void les_io::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void les_io::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void les_io::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}
