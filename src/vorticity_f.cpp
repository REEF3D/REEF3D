/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"vorticity_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

vorticity_f::vorticity_f(lexer *p, fdm *a) : strain(p,a), omega1(p), omega2(p), omega3(p)
{
}

vorticity_f::~vorticity_f()
{
}

void vorticity_f::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{
    double wx,wy,wz;

// xy plane

    LOOP
    {
     omega1(i,j,k) = qij(p,a,2,3);
     omega2(i,j,k) = qij(p,a,1,3);
     omega3(i,j,k) = qij(p,a,2,1);
    }

    pgc->start4(p,omega1,1);
	pgc->start4(p,omega2,1);
	pgc->start4(p,omega3,1);
    pgc->dgcpol(p,omega1,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,omega2,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,omega3,p->dgc4,p->dgc4_count,14);

    // --
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
    ffn=float(p->ipol4(omega1));

	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;n++)
	{
    ffn=float(p->ccipol4(omega1,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));

	result.write((char*)&ffn, sizeof (float));
	}
	
	// --
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
    ffn=float(p->ipol4(omega2));

	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;n++)
	{
    ffn=float(p->ccipol4(omega2,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));

	result.write((char*)&ffn, sizeof (float));
	}
	
	// --
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
    ffn=float(p->ipol4(omega3));

	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;n++)
	{
    ffn=float(p->ccipol4(omega3,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));

	result.write((char*)&ffn, sizeof (float));
	}


}

void vorticity_f::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"vorticity x\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"vorticity y\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"vorticity z\"/>"<<endl;
}

void vorticity_f::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"vorticity x\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"vorticity y\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"vorticity z\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void vorticity_f::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
}





/*
vorticity_f::vorticity_f(lexer *p, fdm *a)
{
}

vorticity_f::~vorticity_f()
{
}

void vorticity_f::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{
    // xy plane
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
	val = 0.25*(a->u(i,j,k) + a->u(i,j+1,k) + a->u(i,j,k+1) + a->u(i,j+1,k+1))*p->DXM
        - 0.25*(a->v(i,j,k) + a->v(i+1,j,k) + a->u(i,j,k+1) + a->u(i+1,j,k+1))*p->DXM;

    ffn=float(val);

	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;n++)
	{
	val = p->ccipol1(a->u,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2])*p->DXM
        - p->ccipol2(a->v,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2])*p->DXM;

    ffn=float(val);

	result.write((char*)&ffn, sizeof (float));
	}

	// xz plane
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
	val = -0.25*(a->u(i,j,k) + a->u(i,j+1,k) + a->u(i,j,k+1) + a->u(i,j+1,k+1))*p->DXM
        + 0.25*(a->w(i,j,k) + a->w(i+1,j,k) + a->w(i,j+1,k) + a->w(i+1,j+1,k))*p->DXM;

    ffn=float(val);

	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;n++)
	{
	val = -p->ccipol1(a->u,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2])*p->DXM
        + p->ccipol3(a->w,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2])*p->DXM;

    ffn=float(val);

	result.write((char*)&ffn, sizeof (float));
	}

	// yz plane
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));

    TPLOOP
	{
	val = -0.25*(a->v(i,j,k) + a->v(i+1,j,k) + a->u(i,j,k+1) + a->u(i+1,j,k+1))*p->DXM
        + 0.25*(a->w(i,j,k) + a->w(i,j+1,k) + a->w(i+1,j,k) + a->w(i+1,j+1,k))*p->DXM;

    ffn=float(val);

	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;n++)
	{
	val = -p->ccipol2(a->v,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2])*p->DXM
        + p->ccipol3(a->w,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2])*p->DXM;

    ffn=float(val);

	result.write((char*)&ffn, sizeof (float));
	}

}

void vorticity_f::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"vorticity Zn\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"vorticity Yn\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"vorticity Xn\"/>"<<endl;
}

void vorticity_f::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"vorticity Zn\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"vorticity Yn\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"vorticity Xn\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;

}

void vorticity_f::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
}


*/


