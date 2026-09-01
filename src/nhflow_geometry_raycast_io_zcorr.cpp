/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"nhflow_geometry.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

// Inside/outside flagging by vertical ray parity.
//
// The ray is vertical, so the 3D segment-triangle intersection reduces to a
// point-in-triangle test on the x-y projection of the triangle, plus a
// barycentric evaluation of the crossing height Rz. Triangles that are
// vertical project to zero area and cannot be crossed by a vertical ray -
// they are skipped by the area test.
//
// CL counts the crossings above a cell, CR the crossings below it. A cell is
// inside the entity when both counts are odd; requiring both to agree makes
// the test fail safe to "outside" on a leaky STL.

void nhflow_geometry::ray_cast_io_zcorr(lexer *p, fdm_nhf *d, ghostcell *pgc, int ts, int te)
{
	double xs,xe,ys,ye;
	double Ax,Ay,Az;
	double Bx,By,Bz;
	double Cx,Cy,Cz;
	double abx,aby,bcx,bcy,cax,cay;
	double e0,e1,e2;
	double area2,sgn,denom;
	double Px,Py,Rz;
	int is,ie,js,je;

    const double eps_area = 1.0e-20*DSM*DSM;
    const double margin   = 2.0*DSM;

    LOOP
	{
	CL[IJK]=0;
	CR[IJK]=0;
	}


	for(n=ts; n<te; ++n)
	{
	Ax = tri_x[n][0];
	Ay = tri_y[n][0];
	Az = tri_z[n][0];

	Bx = tri_x[n][1];
	By = tri_y[n][1];
	Bz = tri_z[n][1];

	Cx = tri_x[n][2];
	Cy = tri_y[n][2];
	Cz = tri_z[n][2];

    // triangle footprint in x-y
	xs = MIN3(Ax,Bx,Cx);
	xe = MAX3(Ax,Bx,Cx);

	ys = MIN3(Ay,By,Cy);
	ye = MAX3(Ay,By,Cy);

    // AABB reject: footprint does not overlap this subdomain
    if(xe < p->originx-margin || xs > p->endx+margin)
    continue;

    if(p->j_dir==1 && (ye < p->originy-margin || ys > p->endy+margin))
    continue;

    // projected edge vectors, loop invariant
    abx = Bx-Ax;
    aby = By-Ay;

    bcx = Cx-Bx;
    bcy = Cy-By;

    cax = Ax-Cx;
    cay = Ay-Cy;

    // twice the signed area of the projected triangle: (B-A) x (C-A)
    area2 = -abx*cay + aby*cax;

    // vertical triangle: cannot be crossed by a vertical ray
    if(fabs(area2) < eps_area)
    continue;

    // normalise the winding to counter-clockwise
    sgn   = (area2 > 0.0 ? 1.0 : -1.0);
    denom = 1.0/(sgn*area2);

	is = MAX(p->posc_i(xs)-1, 0);
	ie = MIN(p->posc_i(xe)+2, p->knox);

    js = 0;
    je = 1;

    if(p->j_dir==1)
    {
	js = MAX(p->posc_j(ys)-1, 0);
	je = MIN(p->posc_j(ye)+2, p->knoy);
    }


		for(i=is;i<ie;i++)
		for(j=js;j<je;j++)
		{
		Px = p->XP[IP];
		Py = p->YP[JP];

        // projected edge functions: e0 opposite C, e1 opposite A, e2 opposite B
        e0 = sgn*(abx*(Py-Ay) - aby*(Px-Ax));
        e1 = sgn*(bcx*(Py-By) - bcy*(Px-Bx));
        e2 = sgn*(cax*(Py-Cy) - cay*(Px-Cx));

        if(e0 < 0.0 || e1 < 0.0 || e2 < 0.0)
        continue;

        // barycentric interpolation of the crossing height
        Rz = (e1*Az + e2*Bz + e0*Cz)*denom;

            for(k=0;k<p->knoz;++k)
            {
				if(p->ZSP[IJK] < Rz)
				CL[IJK] += 1;

                else
				CR[IJK] += 1;
            }
		}
	}


    LOOP
	if((CL[IJK]+1)%2==0  && (CR[IJK]+1)%2==0)
	IO[IJK]=-1;
}