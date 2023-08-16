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

#include"rans_io.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

rans_io::rans_io(lexer *p, fdm *a) : strain(p,a), eps(p), kin(p), eddyv0(p), wallf(p),
									 ke_c_1e(1.44), ke_c_2e(1.92),ke_sigma_k(1.0),ke_sigma_e(1.3),
									 kw_alpha(5.0/9.0), kw_beta(3.0/40.0),kw_sigma_k(2.0),kw_sigma_w(2.0),
									 sst_alpha1(5.0/9.0), sst_alpha2(0.44), sst_beta1(3.0/40.0), sst_beta2(0.0828), 
									 sst_sigma_k1(0.85), sst_sigma_k2(1.0), sst_sigma_w1(0.5), sst_sigma_w2(0.856)
{
}

rans_io::~rans_io()
{
}

void rans_io::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
	ffn=float(p->ipol4_a(kin));
	result.write((char*)&ffn, sizeof (float));
	}
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));

	TPLOOP
	{
	ffn=float(p->ipol4_a(eps));
	result.write((char*)&ffn, sizeof (float));
	}

}

double rans_io::ccipol_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=p->ccipol4( kin, xp, yp, zp);

    return val;
}

double rans_io::ccipol_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=p->ccipol4( eps, xp, yp, zp);

    return val;
}

double rans_io::ccipol_a_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=p->ccipol4_kin( kin, xp, yp, zp);

    return val;
}

double rans_io::ccipol_a_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val;

    val=p->ccipol4_a( eps, xp, yp, zp);

    return val;
}

double rans_io::kinval(int ii, int jj, int kk)
{
    double val;

    val=kin(ii,jj,kk);

    return val;
}

double rans_io::epsval(int ii, int jj, int kk)
{
    double val;

    val=eps(ii,jj,kk);

    return val;
}

void rans_io::kinget(int ii, int jj, int kk,double val)
{
    kin(ii,jj,kk)=val;
}

void rans_io::epsget(int ii, int jj, int kk,double val)
{
    eps(ii,jj,kk)=val;
}

void rans_io::gcupdate(lexer *p, fdm *a, ghostcell *pgc)
{
    pgc->start4(p,kin,20);
	pgc->start4(p,eps,30);
}

void rans_io::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"kin\"/>"<<endl;
	
	if(p->T10==1||p->T10==11 || p->T10==21 ||p->T10==0 || p->T10>30)
	result<<"<PDataArray type=\"Float32\" Name=\"epsilon\"/>"<<endl;
	if(p->T10==2||p->T10==12 || p->T10==22||p->T10==3||p->T10==13)
    result<<"<PDataArray type=\"Float32\" Name=\"omega\"/>"<<endl;
}

void rans_io::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"kin\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	if(p->T10==1||p->T10==11 || p->T10==21 ||p->T10==0 || p->T10>30)
	result<<"<DataArray type=\"Float32\" Name=\"epsilon\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	if(p->T10==2||p->T10==12 || p->T10==22||p->T10==3||p->T10==13)
    result<<"<DataArray type=\"Float32\" Name=\"omega\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void rans_io::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}

