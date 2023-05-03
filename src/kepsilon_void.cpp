/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"kepsilon_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"
#include"convection.h"

kepsilon_void::kepsilon_void(lexer* p, fdm* a, ghostcell *pgc)
{
}

kepsilon_void::~kepsilon_void()
{
}

void kepsilon_void::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow, vrans* pvrans)
{
}
void kepsilon_void::isource(lexer* p, fdm* a)
{
	ULOOP
	a->F(i,j,k)=0.0;
}

void kepsilon_void::jsource(lexer *p,fdm* a)
{
	VLOOP
	a->G(i,j,k)=0.0;
}

void kepsilon_void::ksource(lexer *p,fdm* a)
{
	WLOOP
	a->H(i,j,k)=0.0;
}

void kepsilon_void::ktimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

void kepsilon_void::etimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

void kepsilon_void::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{

}

double kepsilon_void::kinval(int ii, int jj, int kk)
{
    double val=0.0;

    return val;
}

double kepsilon_void::epsval(int ii, int jj, int kk)
{
    double val=0.0;

    return val;
}

double kepsilon_void::ccipol_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double kepsilon_void::ccipol_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double kepsilon_void::ccipol_a_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double kepsilon_void::ccipol_a_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

void kepsilon_void::kinget(int ii, int jj, int kk,double val)
{

}

void kepsilon_void::epsget(int ii, int jj, int kk,double val)
{

}

void kepsilon_void::gcupdate(lexer *p, fdm *a, ghostcell *pgc)
{
}

void kepsilon_void::ini(lexer* p, fdm*a, ghostcell* pgc)
{
}

void kepsilon_void::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void kepsilon_void::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void kepsilon_void::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

