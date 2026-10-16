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

#include"6DOF_obj_2D.h"
#include"lexer.h"
#include"ghostcell.h"
#include<vector>
#include<cmath>

// The 2D body level set of SFLOW (fs: signed horizontal distance to the waterline of the hull, negative inside;
// draft: depth of the hull below the still water level) on the curvilinear grid (G 50). The Cartesian version
// casts rays along the grid lines x and y and re-initialises the distance on the grid; on the curvilinear grid
// the waterline (the section of the hull with the plane z = wd) is cut into segments and every cell centre
// gets
//   inside/outside: parity of the waterline segments crossed by the ray from the centre in +x
//   fs:             exact distance to the nearest segment (clamped to 10 cell widths), signed by inside/outside
//   draft:          the deepest crossing of the vertical ray of the centre with the hull below wd
// All ranks hold all triangles: the parity and the distance of a cell need no exchange.

void sixdof_obj_2D::ray_cast_2D_curv(lexer *p, ghostcell *pgc)
{
    if(curv_dsm<=0.0)
    {
        double ws=0.0;
        int nw=0;
        SLICELOOP4
        {
        ws += p->geo_width(i,j);
        ++nw;
        }
        ws = pgc->globalsum(ws);
        nw = pgc->globalisum(nw);
        curv_dsm = ws/double(MAX(nw,1));
    }

    const double dsm = curv_dsm;
    const double band = 10.0*dsm;
    const double psi = 1.0e-8*dsm;
    const double z0 = p->wd + psi;

    SLICELOOP4
	{
        fsio(i,j)=1;
	    fs(i,j)=1.0e9;
        Ls(i,j) = 1.0e-6;
        Bs(i,j) = 1.0e-6;
        draft(i,j) = 0.0;
        cl(i,j) = 0;
	}

    double bx0,bx1,by0,by1;
    p->geo_box(bx0,bx1,by0,by1);

    std::vector<int> cand;

	for(int n=0; n<tricount; ++n)
	{
        const double X[3] = {tri_x[n][0], tri_x[n][1], tri_x[n][2]};
        const double Y[3] = {tri_y[n][0], tri_y[n][1], tri_y[n][2]};
        const double Z[3] = {tri_z[n][0], tri_z[n][1], tri_z[n][2]};

        const double xs = MIN3(X[0],X[1],X[2]), xe = MAX3(X[0],X[1],X[2]);
        const double ys = MIN3(Y[0],Y[1],Y[2]), ye = MAX3(Y[0],Y[1],Y[2]);

        // ---- draft: vertical rays of the cell centres under the triangle
        const double area2 = (X[1]-X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0]);

        if(fabs(area2) > 1.0e-20*dsm*dsm && !(xe < bx0 || xs > bx1 || ye < by0 || ys > by1))
        {
            p->geo_cells(xs, xe, ys, ye, cand);

            for(size_t c=0; c<cand.size(); c+=2)
            {
                i = cand[c];
                j = cand[c+1];
                const double Px = p->XC2D[IJ] - psi, Py = p->YC2D[IJ] + psi;

                const double l0 = ((X[1]-Px)*(Y[2]-Py) - (X[2]-Px)*(Y[1]-Py))/area2;
                const double l1 = ((X[2]-Px)*(Y[0]-Py) - (X[0]-Px)*(Y[2]-Py))/area2;
                const double l2 = 1.0 - l0 - l1;

                if(l0<0.0 || l1<0.0 || l2<0.0)
                continue;

                const double Rz = l0*Z[0] + l1*Z[1] + l2*Z[2];

                if(p->wd-Rz>0.0)
                draft(i,j) = MAX(p->wd-Rz, draft(i,j));
            }
        }

        // ---- the waterline segment of the triangle (section with z = wd + psi)
        double sx[2], sy[2];
        int ns=0;

        for(int q=0; q<3 && ns<2; ++q)
        {
            const int r = (q+1)%3;
            const double dq = Z[q]-z0, dr = Z[r]-z0;

            if((dq<0.0 && dr>=0.0) || (dq>=0.0 && dr<0.0))
            {
                const double t = dq/(dq-dr);
                sx[ns] = X[q] + t*(X[r]-X[q]);
                sy[ns] = Y[q] + t*(Y[r]-Y[q]);
                ++ns;
            }
        }

        if(ns<2)
        continue;

        const double ax=sx[0], ay=sy[0], bx=sx[1], by=sy[1];
        const double ex=bx-ax, ey=by-ay;
        const double ll = ex*ex + ey*ey;

        // parity: the segment crosses the ray of the centre (yc in [min y, max y), half-open) to the right
        if(ay!=by)
        {
            const double ylo = MIN(ay,by), yhi = MAX(ay,by);

            if(!(yhi < by0 || ylo > by1) && MAX(ax,bx) >= bx0)
            {
                p->geo_cells(bx0-1.0, MAX(ax,bx), ylo, yhi, cand);

                for(size_t c=0; c<cand.size(); c+=2)
                {
                    i = cand[c];
                    j = cand[c+1];
                    const double xc = p->XC2D[IJ], yc = p->YC2D[IJ];

                    if(yc<ylo || yc>=yhi)
                    continue;

                    const double xi = ax + (yc-ay)*ex/ey;

                    if(xc < xi)
                    cl(i,j) += 1;
                }
            }
        }

        // distance to the segment in the band
        if(!(MAX(ax,bx)+band < bx0 || MIN(ax,bx)-band > bx1 || MAX(ay,by)+band < by0 || MIN(ay,by)-band > by1))
        {
            p->geo_cells(MIN(ax,bx)-band, MAX(ax,bx)+band, MIN(ay,by)-band, MAX(ay,by)+band, cand);

            for(size_t c=0; c<cand.size(); c+=2)
            {
                i = cand[c];
                j = cand[c+1];
                const double px = p->XC2D[IJ]-ax, py = p->YC2D[IJ]-ay;
                const double t = ll>0.0 ? MAX(0.0, MIN(1.0, (px*ex + py*ey)/ll)) : 0.0;
                const double dx = px - t*ex, dy = py - t*ey;
                const double d = sqrt(dx*dx + dy*dy);

                if(d < fs(i,j))
                fs(i,j) = d;
            }
        }
    }

	SLICELOOP4
    {
        if(cl(i,j)%2==1)
        fsio(i,j) = -1;

        fs(i,j) = MIN(fs(i,j), band);

        if(fsio(i,j)==-1)
        fs(i,j) = -fs(i,j);
    }

	pgc->gcsl_start4(p,fs,50);
    pgc->gcsl_start4(p,draft,50);
}
