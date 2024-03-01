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

#include"nhflow_komega_void.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"strain.h"

nhflow_komega_void::nhflow_komega_void(lexer* p, fdm_nhf* d, ghostcell *pgc)
{
}

nhflow_komega_void::~nhflow_komega_void()
{
}

void nhflow_komega_void::start(fdm_nhf* d, lexer* p, nhflow_convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow, vrans* pvrans)
{
}
void nhflow_komega_void::isource(lexer* p, fdm_nhf* d)
{
	LOOP
	d->F[IJK]=0.0;
}

void nhflow_komega_void::jsource(lexer *p,fdm_nhf* d)
{
	LOOP
	d->G[IJK]=0.0;
}

void nhflow_komega_void::ksource(lexer *p,fdm_nhf* d)
{
	LOOP
	d->H[IJK]=0.0;
}

void nhflow_komega_void::ktimesave(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
}

void nhflow_komega_void::etimesave(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
}

void nhflow_komega_void::print_3D(lexer* p, fdm_nhf *d, ghostcell *pgc, ofstream &result)
{

}

double nhflow_komega_void::kinval(int ii, int jj, int kk)
{
    double val=0.0;

    return val;
}

double nhflow_komega_void::epsval(int ii, int jj, int kk)
{
    double val=0.0;

    return val;
}

double nhflow_komega_void::ccipol_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double nhflow_komega_void::ccipol_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double nhflow_komega_void::ccipol_a_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

double nhflow_komega_void::ccipol_a_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=0.0;

    return val;
}

void nhflow_komega_void::kinget(int ii, int jj, int kk,double val)
{

}

void nhflow_komega_void::epsget(int ii, int jj, int kk,double val)
{

}

void nhflow_komega_void::gcupdate(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void nhflow_komega_void::ini(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
}

void nhflow_komega_void::name_pvtu(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result)
{
}

void nhflow_komega_void::name_vtu(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void nhflow_komega_void::offset_vtu(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

