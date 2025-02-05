/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_les_io.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_les_io::nhflow_les_io(lexer *p, fdm_nhf *d) : nhflow_strain(p,d)
{

}

nhflow_les_io::~nhflow_les_io()
{
}

void nhflow_les_io::print_3D(lexer* p, fdm_nhf *d, ghostcell *pgc, ofstream &result)
{
    
    // eddyv
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
	if(p->j_dir==0)
    {
    jj=j;
    j=0;
	ffn=float(0.5*(d->EV[IJK]+d->EV[IJKp1]));
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.25*(d->EV[IJK]+d->EV[IJKp1]+d->EV[IJp1K]+d->EV[IJp1Kp1]));
        
        
	result.write((char*)&ffn, sizeof (float));
	}

}

double nhflow_les_io::ccipol_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4( kin, xp, yp, zp);

    return val;
}

double nhflow_les_io::ccipol_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4( eps, xp, yp, zp);

    return val;
}

double nhflow_les_io::ccipol_a_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4a( kin, xp, yp, zp);

    return val;
}

double nhflow_les_io::ccipol_a_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4a( eps, xp, yp, zp);

    return val;
}

double nhflow_les_io::kinval(int ii, int jj, int kk)
{
    double val=0.0;

    //val=kin(ii,jj,kk);

    return val;
}

double nhflow_les_io::epsval(int ii, int jj, int kk)
{
    double val=0.0;

    //val=eps(ii,jj,kk);

    return val;
}

void nhflow_les_io::kinget(int ii, int jj, int kk,double val)
{
}

void nhflow_les_io::epsget(int ii, int jj, int kk,double val)
{
}

void nhflow_les_io::gcupdate(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void nhflow_les_io::name_pvtu(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>"<<endl;
}

void nhflow_les_io::name_vtu(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void nhflow_les_io::offset_vtu(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}


void nhflow_les_io::ini(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
}


void nhflow_les_io::plain_wallfunc(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
}

void nhflow_les_io::inflow(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
}

void nhflow_les_io::tau_calc(fdm_nhf *d, lexer* p, double maxwdist)
{
}

void nhflow_les_io::isource(lexer *p, fdm_nhf *d)
{
	LOOP
	d->F[IJK]=0.0;
}

void nhflow_les_io::jsource(lexer *p, fdm_nhf *d)
{
	LOOP
	d->G[IJK]=0.0;
}

void nhflow_les_io::ksource(lexer *p, fdm_nhf *d)
{
	LOOP
	d->H[IJK]=0.0;
}