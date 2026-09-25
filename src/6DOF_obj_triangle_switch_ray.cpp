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

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fieldint.h"
#include<vector>
#include<cmath>
#include<cstddef>

namespace
{
// Uniform 2D bin grid over the projection plane (a,b) of one ray direction.
// Every triangle is stored in all bins overlapped by its psi-expanded bounding
// box. A point lookup therefore returns a superset of the triangles that pass
// the bounding-box test of the ray cast, so the crossing count is identical to
// testing against all triangles, at O(candidates) instead of O(tricount).
struct tri_bin2d
{
    double amin=0.0, bmin=0.0, ha=1.0, hb=1.0;
    int na=1, nb=1;
    std::vector<std::size_t> start;
    std::vector<int> id;

    inline int ia(double a) const
    {
        int i = int((a-amin)/ha);
        return i<0 ? 0 : (i>=na ? na-1 : i);
    }

    inline int ib(double b) const
    {
        int i = int((b-bmin)/hb);
        return i<0 ? 0 : (i>=nb ? nb-1 : i);
    }

    inline int cell(double a, double b) const
    {
        return ia(a)*nb + ib(b);
    }

    void build(const std::vector<double> &alo, const std::vector<double> &ahi,
               const std::vector<double> &blo, const std::vector<double> &bhi, int num)
    {
        if(num<=0)
        {
            start.assign(2,0);
            return;
        }

        amin = alo[0];
        bmin = blo[0];
        double amax = ahi[0];
        double bmax = bhi[0];

        for(int n=1; n<num; ++n)
        {
            amin = alo[n]<amin ? alo[n] : amin;
            bmin = blo[n]<bmin ? blo[n] : bmin;
            amax = ahi[n]>amax ? ahi[n] : amax;
            bmax = bhi[n]>bmax ? bhi[n] : bmax;
        }

        double La = amax-amin;
        double Lb = bmax-bmin;

        // about one bin per triangle, split according to the aspect ratio
        const int nmax = 4096;
        double target = double(num);

        if(La>0.0 && Lb>0.0)
        {
            na = int(std::sqrt(target*La/Lb)+0.5);
            na = na<1 ? 1 : (na>nmax ? nmax : na);
            nb = int(target/double(na)+0.5);
            nb = nb<1 ? 1 : (nb>nmax ? nmax : nb);
        }
        else
        {
            na = La>0.0 ? (num<nmax ? num : nmax) : 1;
            nb = Lb>0.0 ? (num<nmax ? num : nmax) : 1;
        }

        ha = La>0.0 ? La/double(na) : 1.0;
        hb = Lb>0.0 ? Lb/double(nb) : 1.0;

        const std::size_t ncell = std::size_t(na)*std::size_t(nb);
        start.assign(ncell+1,0);

        // pass 1: count
        for(int n=0; n<num; ++n)
        {
            int ia0=ia(alo[n]), ia1=ia(ahi[n]);
            int ib0=ib(blo[n]), ib1=ib(bhi[n]);

            for(int ii=ia0; ii<=ia1; ++ii)
            for(int jj=ib0; jj<=ib1; ++jj)
            ++start[std::size_t(ii)*nb + jj + 1];
        }

        for(std::size_t c=0; c<ncell; ++c)
        start[c+1] += start[c];

        // pass 2: fill
        id.resize(start[ncell]);
        std::vector<std::size_t> pos(start.begin(), start.end()-1);

        for(int n=0; n<num; ++n)
        {
            int ia0=ia(alo[n]), ia1=ia(ahi[n]);
            int ib0=ib(blo[n]), ib1=ib(bhi[n]);

            for(int ii=ia0; ii<=ia1; ++ii)
            for(int jj=ib0; jj<=ib1; ++jj)
            id[pos[std::size_t(ii)*nb + jj]++] = n;
        }
    }
};
}

void sixdof_obj::triangle_switch_ray(lexer *p, ghostcell *pgc)
{
	double Px,Py,Pz;
	double Qx,Qy,Qz;
    double Sx=0.0,Sy=0.0,Sz=0.0;
    double Tx=0.0,Ty=0.0,Tz=0.0;
	double Rx,Ry,Rz;
	double Ax,Ay,Az;
	double Bx,By,Bz;
	double Cx,Cy,Cz;
	double PQx,PQy,PQz;
	double Mx,My,Mz;
	double u,v,w;
	double denom;
	double psi = 1.0e-8*p->DXM;

    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double xc,yc,zc;
	double nx,ny,nz,norm;
    int tricount_local_max,sum;
    int cutnum;
    int domdir;

    if(p->mpirank==0)
	cout<<"Triangle Switch "<<endl;

    // allocate
    p->Iarray(tricount_local_list,p->M10+1);
    p->Iarray(tricount_local_displ,p->M10+1);

    p->Iarray(tri_switch,tricount);

    for (n=0;n<tricount;++n)
    tri_switch[n]=0;

    // divide triangles to local processors
    tricount_local = int(tricount/p->M10);

    sum = 0;
    for(q=0; q<p->M10-1; ++q)
    {
    tricount_local_list[q]=tricount_local;
    sum += tricount_local;
    }

    tricount_local_list[p->M10-1] = tricount - sum;


    // displacement
    tricount_local_displ[0]=0;

    for(q=1;q<p->M10+1;++q)
    tricount_local_displ[q] = tricount_local_displ[q-1] + tricount_local_list[q-1];

    tricount_local_max = 1;
    for(q=0; q<p->M10; ++q)
    tricount_local_max = MAX(tricount_local_list[q],tricount_local_max);

    p->Iarray(tri_switch_local,tricount_local_max);


    // triangle bounding boxes (psi-expanded), computed once
    std::vector<double> bxs(tricount),bxe(tricount);
    std::vector<double> bys(tricount),bye(tricount);
    std::vector<double> bzs(tricount),bze(tricount);

    for(n=0; n<tricount; ++n)
    {
    bxs[n] = MIN3(tri_x[n][0],tri_x[n][1],tri_x[n][2]) - psi;
    bxe[n] = MAX3(tri_x[n][0],tri_x[n][1],tri_x[n][2]) + psi;

    bys[n] = MIN3(tri_y[n][0],tri_y[n][1],tri_y[n][2]) - psi;
    bye[n] = MAX3(tri_y[n][0],tri_y[n][1],tri_y[n][2]) + psi;

    bzs[n] = MIN3(tri_z[n][0],tri_z[n][1],tri_z[n][2]) - psi;
    bze[n] = MAX3(tri_z[n][0],tri_z[n][1],tri_z[n][2]) + psi;
    }

    // bin grids for the three ray directions
    tri_bin2d grid_x, grid_y, grid_z;

    grid_x.build(bys,bye,bzs,bze,tricount);    // x-ray: (y,z)-plane

    if(p->y_dir==1)
    grid_y.build(bxs,bxe,bzs,bze,tricount);    // y-ray: (x,z)-plane

    grid_z.build(bxs,bxe,bys,bye,tricount);    // z-ray: (x,y)-plane


    // ray cast
	for(q=tricount_local_displ[p->mpirank];q<tricount_local_displ[p->mpirank+1];++q)
	{
        // triangle points
        x0 = tri_x[q][0];
        y0 = tri_y[q][0];
        z0 = tri_z[q][0];

        x1 = tri_x[q][1];
        y1 = tri_y[q][1];
        z1 = tri_z[q][1];

        x2 = tri_x[q][2];
        y2 = tri_y[q][2];
        z2 = tri_z[q][2];

        // normals
        nx = (y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0);
        ny = (x2 - x0)*(z1 - z0) - (x1 - x0)*(z2 - z0);
        nz = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

        norm = sqrt(nx*nx + ny*ny + nz*nz);

        nx /= norm>1.0e-20?norm:1.0e20;
        ny /= norm>1.0e-20?norm:1.0e20;
        nz /= norm>1.0e-20?norm:1.0e20;

        domdir = 0;

        if(fabs(nx)>=fabs(ny))
        {
        domdir=1;

            if(fabs(nx)<fabs(nz))
            domdir=3;
        }

        if(fabs(nx)<fabs(ny))
        {
        domdir=2;

            if(fabs(ny)<fabs(nz))
            domdir=3;
        }

        // Center of triangle
		xc = (x0 + x1 + x2)/3.0;
		yc = (y0 + y1 + y2)/3.0;
		zc = (z0 + z1 + z2)/3.0;

        cutnum=0;

        // no y-rays in 2D
        if(domdir==2 && p->y_dir==0)
        domdir=0;

        // degenerate triangle (zero normal): never counted a crossing before either
        if((domdir==1 && nx==0.0) || (domdir==2 && ny==0.0) || (domdir==3 && nz==0.0))
        domdir=0;

        if(domdir>0)
        {
        // ray and half-ray bounds, only depend on q
        Px = Qx = xc;
        Py = Qy = yc;
        Pz = Qz = zc;

        if(domdir==1)
        {
        Px = p->global_xmin-10.0*p->DXM;
        Qx = p->global_xmax+10.0*p->DXM;

        Sx = nx>0.0 ? xc : p->global_xmin-10.0*p->DXM;
        Tx = nx>0.0 ? p->global_xmax+10.0*p->DXM : xc;
        }

        if(domdir==2)
        {
        Py = p->global_ymin-1000.0*p->DXM;
        Qy = p->global_ymax+1000.0*p->DXM;

        Sy = ny>0.0 ? yc : p->global_ymin-1000.0*p->DXM;
        Ty = ny>0.0 ? p->global_ymax+1000.0*p->DXM : yc;
        }

        if(domdir==3)
        {
        Pz = p->global_zmin-1000.0*p->DXM;
        Qz = p->global_zmax+1000.0*p->DXM;

        Sz = nz>0.0 ? zc : p->global_zmin-1000.0*p->DXM;
        Tz = nz>0.0 ? p->global_zmax+1000.0*p->DXM : zc;
        }

        PQx = Qx-Px;
        PQy = Qy-Py;
        PQz = Qz-Pz;

        Mx = PQy*Pz - PQz*Py;
        My = PQz*Px - PQx*Pz;
        Mz = PQx*Py - PQy*Px;

        // candidate triangles from the bin containing the ray
        const tri_bin2d *grid = &grid_z;
        double pa = xc;
        double pb = yc;

        if(domdir==1)
        {
        grid = &grid_x;
        pa = yc;
        pb = zc;
        }

        if(domdir==2)
        {
        grid = &grid_y;
        pa = xc;
        pb = zc;
        }

        const int c = grid->cell(pa,pb);

            // ray cast loop
            for(std::size_t e=grid->start[c]; e<grid->start[c+1]; ++e)
            {
            n = grid->id[e];

            if(n==q)
            continue;

            // bounding box test
            if(!((domdir==1 && yc>=bys[n] && yc<=bye[n] && zc>=bzs[n] && zc<=bze[n])
              || (domdir==2 && xc>=bxs[n] && xc<=bxe[n] && zc>=bzs[n] && zc<=bze[n])
              || (domdir==3 && xc>=bxs[n] && xc<=bxe[n] && yc>=bys[n] && yc<=bye[n])))
            continue;

            Ax = tri_x[n][0];
            Ay = tri_y[n][0];
            Az = tri_z[n][0];

            Bx = tri_x[n][1];
            By = tri_y[n][1];
            Bz = tri_z[n][1];

            Cx = tri_x[n][2];
            Cy = tri_y[n][2];
            Cz = tri_z[n][2];

                // uvw
                u = PQx*(Cy*Bz - Cz*By) + PQy*(Cz*Bx - Cx*Bz) + PQz*(Cx*By - Cy*Bx)
                  + Mx*(Cx-Bx) + My*(Cy-By) + Mz*(Cz-Bz);

                v = PQx*(Ay*Cz - Az*Cy) + PQy*(Az*Cx - Ax*Cz) + PQz*(Ax*Cy - Ay*Cx)
                  + Mx*(Ax-Cx) + My*(Ay-Cy) + Mz*(Az-Cz);

                w = PQx*(By*Az - Bz*Ay) + PQy*(Bz*Ax - Bx*Az) + PQz*(Bx*Ay - By*Ax)
                  + Mx*(Bx-Ax) + My*(By-Ay) + Mz*(Bz-Az);


                int check=1;
                if(fabs(u)<1.0e-15 && fabs(v)<1.0e-15 && fabs(w)<1.0e-15)
                check = 0;

                    // note: && binds before ||, so check only guards the negative branch (as before)
                    if((u>=0.0 && v>=0.0 && w>=0.0) || (u<0.0 && v<0.0 && w<0.0) && check==1)
                    {
                    denom = 1.0/(u+v+w);
                    u *= denom;
                    v *= denom;
                    w *= denom;

                        if(domdir==1)
                        {
                        Rx = u*Ax + v*Bx + w*Cx;

                        if(Rx>=Sx && Rx<=Tx)
                        ++cutnum;
                        }

                        if(domdir==2 && p->j_dir==1)
                        {
                        Ry = u*Ay + v*By + w*Cy;

                        if(Ry>=Sy && Ry<=Ty)
                        ++cutnum;
                        }

                        if(domdir==3)
                        {
                        Rz = u*Az + v*Bz + w*Cz;

                        if(Rz>=Sz && Rz<=Tz)
                        ++cutnum;
                        }
                    }
            } // ray cast loop end
        }

        // odd number of crossings -> inside pointing -> switch
        tri_switch_local[q-tricount_local_displ[p->mpirank]] = (cutnum%2!=0) ? 1 : 0;
	}

    pgc->allgatherv_int(tri_switch_local, tricount_local_list[p->mpirank], tri_switch, tricount_local_list, tricount_local_displ);

    // start loop part 3: global switch
        tricount_switch_total=0;
        for (int n=0;n<tricount;++n)
        if(tri_switch[n]==1)
        {
        x1 = tri_x[n][1];
        y1 = tri_y[n][1];
        z1 = tri_z[n][1];

        x2 = tri_x[n][2];
        y2 = tri_y[n][2];
        z2 = tri_z[n][2];


        tri_x[n][1] = x2;
        tri_y[n][1] = y2;
        tri_z[n][1] = z2;

        tri_x[n][2] = x1;
        tri_y[n][2] = y1;
        tri_z[n][2] = z1;

        ++tricount_switch_total;
        }

        if(p->mpirank==0)
        cout<<"6DOF STL triangle switch count: "<<tricount_switch_total<<endl;

        // free allocated arrays
    p->del_Iarray(tri_switch,tricount);
    p->del_Iarray(tricount_local_list,p->M10+1);
    p->del_Iarray(tricount_local_displ,p->M10+1);
    p->del_Iarray(tri_switch_local,tricount_local_max);
}
