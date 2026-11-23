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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#include"nhflow_shipbl.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<iostream>

// The hull triangles read from floating.stl are replaced by the displacement body of the hull (nhflow_shipbl): the
// solver then runs with the body as it does with a displacement body written by Cases/KCS/shipbl.py (single run,
// no first pass). Mass and inertia come from the control file (X 22 - X 24), as with the body file.
//
//   NHF_BL=1          on
//   NHF_BL_L          length between perpendiculars [m]                                   (needed)
//   NHF_BL_XFP        x of the forward perpendicular in the tank, hull at even keel [m]   (needed)
//   NHF_BL_U          speed [m/s]                         default NHF_SRC
//   NHF_BL_T, _B      draft, beam [m]                     default from the triangles
//   NHF_BL_ZW         still water level [m]               default F 60
//   NHF_BL_SINK       sinkage of the triangles [m, + down], NHF_BL_TRIM trim [deg, + bow up] about midship at the
//                     still water line: the boundary layer is computed for the hull at even keel (as shipbl.py does)
//                     and the body is given the attitude of the hull
//   NHF_BL_XTR        trip at s = x/L                     default 0.05
//   NHF_BL_FSW=0      without the surface wake            default 1
//   NHF_BL_DUMP=1     writes REEF3D_NHFLOW_6DOF/NHF_BL_dispbody.stl
//
// Rank 0 computes (seconds; the dense matrix of the double body takes 8 N^2 bytes, N about 8000), the triangles
// are broadcast.
void nhflow_shipbl_apply(lexer *p, ghostcell *pgc, double **&tri_x, double **&tri_y, double **&tri_z,
                         double **&tri_x0, double **&tri_y0, double **&tri_z0, int ts, int &tricount)
{
    auto envd = [](const char *name, double def)
    {
        return getenv(name)!=nullptr ? atof(getenv(name)) : def;
    };

    if(int(envd("NHF_BL",0.0))==0)
    return;

    const double L = envd("NHF_BL_L",0.0);
    const double xfp = envd("NHF_BL_XFP",-1.0e20);
    const double U = envd("NHF_BL_U",envd("NHF_SRC",0.0));

    if(L<=0.0 || xfp<-1.0e19 || U<=0.0)
    {
        if(p->mpirank==0)
        std::cout<<"NHF_BL: NHF_BL_L, NHF_BL_XFP and the speed (NHF_BL_U or NHF_SRC) are needed"<<std::endl;
        pgc->final(true);
    }

    const double zw = envd("NHF_BL_ZW",p->F60);
    const double sink = envd("NHF_BL_SINK",0.0);
    const double trim = envd("NHF_BL_TRIM",0.0)*3.14159265358979323846/180.0;
    const double xm = xfp + 0.5*L;
    const double ct = cos(trim), st = sin(trim);

    int nnew=0;
    std::vector<double> buf;

    if(p->mpirank==0)
    {
        // hull to even keel:  x' = xm + xr cos + zr sin,  z' = zw - xr sin + zr cos - sink
        std::vector<nhflow_shipbl::tri3> hull(tricount-ts), body;
        double ylo=1.0e30, yhi=-1.0e30;
        for(int n=ts; n<tricount; ++n)
        for(int q=0; q<3; ++q)
        {
            const double X = tri_x[n][q]-xm, Z = tri_z[n][q]-zw+sink;
            hull[n-ts][3*q]   = xm + X*ct - Z*st;
            hull[n-ts][3*q+1] = tri_y[n][q];
            hull[n-ts][3*q+2] = zw + X*st + Z*ct;
            ylo = std::min(ylo,tri_y[n][q]);
            yhi = std::max(yhi,tri_y[n][q]);
        }

        nhflow_shipbl bl;
        bl.L = L;
        bl.T = envd("NHF_BL_T",0.0);
        bl.B = envd("NHF_BL_B",0.0);
        bl.U = U;
        bl.nu = p->W2;
        bl.x_fp = xfp;
        bl.y_c = envd("NHF_BL_YC",0.5*(ylo+yhi));
        bl.z_w = zw;
        bl.xtr = envd("NHF_BL_XTR",0.05);
        bl.fsw = int(envd("NHF_BL_FSW",1.0))!=0;
        bl.run(hull,body);

        nnew = int(body.size());
        buf.resize(9*nnew);
        for(int n=0; n<nnew; ++n)
        for(int q=0; q<3; ++q)
        {
            const double xr = body[n][3*q]-xm, zr = body[n][3*q+2]-zw;
            buf[9*n+3*q]   = xm + xr*ct + zr*st;
            buf[9*n+3*q+1] = body[n][3*q+1];
            buf[9*n+3*q+2] = zw - xr*st + zr*ct - sink;
        }

        std::cout<<"NHF_BL: the hull ("<<tricount-ts<<" triangles) is replaced by its displacement body, U "<<U<<" m/s, nu "<<p->W2
                 <<", sinkage "<<sink<<" m, trim "<<envd("NHF_BL_TRIM",0.0)<<" deg"<<std::endl<<bl.log<<std::endl;

        if(int(envd("NHF_BL_DUMP",0.0))!=0)
        {
            FILE *fp = fopen("./REEF3D_NHFLOW_6DOF/NHF_BL_dispbody.stl","w");
            if(fp!=nullptr)
            {
                fprintf(fp,"solid body\n");
                for(int n=0; n<nnew; ++n)
                {
                    fprintf(fp," facet normal 0 0 0\n  outer loop\n");
                    for(int q=0; q<3; ++q)
                    fprintf(fp,"   vertex %.17g %.17g %.17g\n",buf[9*n+3*q],buf[9*n+3*q+1],buf[9*n+3*q+2]);      // exact: a run with this file equals the run that wrote it
                    fprintf(fp,"  endloop\n endfacet\n");
                }
                fprintf(fp,"endsolid body\n");
                fclose(fp);
            }
        }
    }

    pgc->bcast_int(&nnew,1);
    if(p->mpirank!=0)
    buf.resize(9*nnew);
    pgc->bcast_double(buf.data(),9*nnew);

    p->Dresize(tri_x,tricount,ts+nnew,3,3);
    p->Dresize(tri_y,tricount,ts+nnew,3,3);
    p->Dresize(tri_z,tricount,ts+nnew,3,3);
    p->Dresize(tri_x0,tricount,ts+nnew,3,3);
    p->Dresize(tri_y0,tricount,ts+nnew,3,3);
    p->Dresize(tri_z0,tricount,ts+nnew,3,3);

    for(int n=0; n<nnew; ++n)
    for(int q=0; q<3; ++q)
    {
        tri_x[ts+n][q] = buf[9*n+3*q];
        tri_y[ts+n][q] = buf[9*n+3*q+1];
        tri_z[ts+n][q] = buf[9*n+3*q+2];
    }
    tricount = ts+nnew;
}
