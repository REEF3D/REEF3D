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

#include"geo_raycast.h"
#include"lexer.h"

// Cartesian ray casting, taken from the 6DOF floating body geometry
// (6DOF_obj_ray3D_io_x/ycorr/zcorr, 6DOF_obj_ray3D_x/y/z, 6DOF_obj_ray3D_direct).
//
// A ray along the grid line through the cell centres of (j,k) for dir 0, (i,k) for
// dir 1 and (i,j) for dir 2 is intersected with every triangle (Plücker coordinates
// u,v,w: all of one sign = crossing). The ray is tilted by psi against rays through
// edges and vertices.

namespace
{
    // crossing of the ray PQ with ABC; returns 1 and the barycentric weights for a crossing
    inline int cross(const geo_ray &r, double tol,
                     double Ax, double Ay, double Az,
                     double Bx, double By, double Bz,
                     double Cx, double Cy, double Cz,
                     double &u, double &v, double &w)
    {
        const double PQx = r.Qx-r.Px;
        const double PQy = r.Qy-r.Py;
        const double PQz = r.Qz-r.Pz;

        const double Mx = PQy*r.Pz - PQz*r.Py;
        const double My = PQz*r.Px - PQx*r.Pz;
        const double Mz = PQx*r.Py - PQy*r.Px;

        u = PQx*(Cy*Bz - Cz*By) + PQy*(Cz*Bx - Cx*Bz) + PQz*(Cx*By - Cy*Bx)
          + Mx*(Cx-Bx) + My*(Cy-By) + Mz*(Cz-Bz);

        v = PQx*(Ay*Cz - Az*Cy) + PQy*(Az*Cx - Ax*Cz) + PQz*(Ax*Cy - Ay*Cx)
          + Mx*(Ax-Cx) + My*(Ay-Cy) + Mz*(Az-Cz);

        w = PQx*(By*Az - Bz*Ay) + PQy*(Bz*Ax - Bx*Az) + PQz*(Bx*Ay - By*Ax)
          + Mx*(Bx-Ax) + My*(By-Ay) + Mz*(Bz-Az);

        int check=1;

        if(tol==0.0)
        {
            if(u==0.0 && v==0.0 && w==0.0)
            check = 0;

            if(((u>0.0 && v>0.0 && w>0.0) || (u<0.0 && v<0.0 && w<0.0)) && check==1)
            {
                const double denom = 1.0/(u+v+w);
                u *= denom;
                v *= denom;
                w *= denom;
                return 1;
            }
        }
        else
        {
            if(fabs(u)<=tol && fabs(v)<=tol && fabs(w)<=tol)
            check = 0;

            if(((u>tol && v>tol && w>tol) || (u<-tol && v<-tol && w<-tol)) && check==1)
            {
                const double denom = 1.0/(u+v+w);
                u *= denom;
                v *= denom;
                w *= denom;
                return 1;
            }
        }

        return 0;
    }
}

// ray along dir through the centre of the cell (i,j,k)
void geo_raycast::setray(lexer *p, int dir, int i, int j, int k, const geo_cart &g, geo_ray &r)
{
    const double psi = g.psi;

    if(g.raystyle==1)
    {
        // shifted ray (DIVEMesh)
        if(dir==0)
        {
            r.Px = p->global_xmin-g.ext;
            r.Py = p->YP[JP]+psi;
            r.Pz = p->ZP[KP]+psi;
            r.Qx = p->global_xmax+g.ext;
            r.Qy = p->YP[JP]+psi;
            r.Qz = p->ZP[KP]+psi;
        }
        else if(dir==1)
        {
            r.Px = p->XP[IP]+psi;
            r.Py = p->global_ymin-g.ext;
            r.Pz = p->ZP[KP]+psi;
            r.Qx = p->XP[IP]+psi;
            r.Qy = p->global_ymax+g.ext;
            r.Qz = p->ZP[KP]+psi;
        }
        else
        {
            r.Px = p->XP[IP]+psi;
            r.Py = p->YP[JP]+psi;
            r.Pz = p->global_zmin-g.ext;
            r.Qx = p->XP[IP]+psi;
            r.Qy = p->YP[JP]+psi;
            r.Qz = p->global_zmax+g.ext;
        }

        return;
    }

    // tilted ray (6DOF)
    if(dir==0)
    {
        r.Px = p->global_xmin-g.ext;
        r.Py = p->YP[JP]-psi;
        r.Pz = p->ZP[KP]+psi;
        r.Qx = p->global_xmax+g.ext;
        r.Qy = p->YP[JP]+psi;
        r.Qz = p->ZP[KP]-psi;
    }
    else if(dir==1)
    {
        r.Px = p->XP[IP]+psi;
        r.Py = p->global_ymin-g.ext;
        r.Pz = p->ZP[KP]-psi;
        r.Qx = p->XP[IP]-psi;
        r.Qy = p->global_ymax+g.ext;
        r.Qz = p->ZP[KP]+psi;
    }
    else
    {
        r.Px = p->XP[IP]-psi;
        r.Py = p->YP[JP]+psi;
        r.Pz = p->global_zmin-g.ext;
        r.Qx = p->XP[IP]+psi;
        r.Qy = p->YP[JP]-psi;
        r.Qz = p->global_zmax+g.ext;
    }
}

void geo_raycast::cart_io(lexer *p, int dir, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                          const geo_cart &g, int *IO, int *CL, int *CR)
{
    int i,j,k,n;
    int as,ae,bs,be;
    double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
    double u,v,w,R;
    double s1,e1,s2,e2;
    geo_ray r;

    // crossing threshold: exact, DIVEMesh vertical rays 1.0e-20
    const double tol = (g.raystyle==1 && dir==2) ? 1.0e-20 : 0.0;

    // transverse directions of the ray
    const int d1 = (dir==0) ? 1 : 0;
    const int d2 = (dir==2) ? 1 : 2;

    for(i=g.is; i<g.ie; ++i)
    for(j=g.js; j<g.je; ++j)
    for(k=g.ks; k<g.ke; ++k)
    if(!g.interior || p->flag4[IJK]>0)
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

        if(g.clip && !checkin(p,dir,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz))
        continue;

        const double c1[3] = {(d1==0)?Ax:Ay, (d1==0)?Bx:By, (d1==0)?Cx:Cy};
        const double c2[3] = {(d2==1)?Ay:Az, (d2==1)?By:Bz, (d2==1)?Cy:Cz};

        s1 = MIN3(c1[0],c1[1],c1[2]);
        e1 = MAX3(c1[0],c1[1],c1[2]);
        s2 = MIN3(c2[0],c2[1],c2[2]);
        e2 = MAX3(c2[0],c2[1],c2[2]);

        bracket(p,d1,s1,e1,g,as,ae);
        bracket(p,d2,s2,e2,g,bs,be);

        for(int a=as; a<ae; ++a)
        for(int b=bs; b<be; ++b)
        {
            if(dir==0)
            {j=a; k=b;}
            else if(dir==1)
            {i=a; k=b;}
            else
            {i=a; j=b;}
            
            setray(p,dir,i,j,k,g,r);

            if(cross(r,tol,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,u,v,w)==0)
            continue;

            if(dir==0)
            {
                R = u*Ax + v*Bx + w*Cx;

                for(i=g.is; i<g.ie; ++i)
                {
                    if(p->XP[IP]<R)
                    CR[IJK] += 1;

                    if(p->XP[IP]>=R)
                    CL[IJK] += 1;
                }
            }
            else if(dir==1)
            {
                R = u*Ay + v*By + w*Cy;

                for(j=g.js; j<g.je; ++j)
                {
                    if(p->YP[JP]<R)
                    CR[IJK] += 1;

                    if(p->YP[JP]>=R)
                    CL[IJK] += 1;
                }
            }
            else
            {
                R = u*Az + v*Bz + w*Cz;

                for(k=g.ks; k<g.ke; ++k)
                {
                    if(p->ZP[KP]<R)
                    CR[IJK] += 1;

                    if(p->ZP[KP]>=R)
                    CL[IJK] += 1;
                }
            }
        }
    }

    for(i=g.is; i<g.ie; ++i)
    for(j=g.js; j<g.je; ++j)
    for(k=g.ks; k<g.ke; ++k)
    if(!g.interior || p->flag4[IJK]>0)
    {
        if(g.raymode==2)
        {
            if(CL[IJK]%2==0 && CR[IJK]%2==0)
            IO[IJK]=-1;
        }
        else
        {
            if((CL[IJK]+1)%2==0 && (CR[IJK]+1)%2==0)
            IO[IJK]=-1;
        }
    }
}

void geo_raycast::cart_dist(lexer *p, int dir, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                            const geo_cart &g, const int *IO, double *LS, double *bed)
{
    int i,j,k,n;
    int as,ae,bs,be;
    int ii,jj,kk;
    double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
    double u,v,w,R;
    double s1,e1,s2,e2;
    geo_ray r;

    // crossing threshold: 6DOF 1.0e-20 for x and z, exact for y; DIVEMesh 1.0e-20 for z only
    const double tol = (g.raystyle==1) ? ((dir==2) ? 1.0e-20 : 0.0) : ((dir==1) ? 0.0 : 1.0e-20);

    // the cell of a crossing is checked against its neighbours where IO is known
    const int off = g.interior ? 0 : 1;

    const int d1 = (dir==0) ? 1 : 0;
    const int d2 = (dir==2) ? 1 : 2;

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

        if(g.clip && !checkin(p,dir,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz))
        continue;

        const double c1[3] = {(d1==0)?Ax:Ay, (d1==0)?Bx:By, (d1==0)?Cx:Cy};
        const double c2[3] = {(d2==1)?Ay:Az, (d2==1)?By:Bz, (d2==1)?Cy:Cz};

        s1 = MIN3(c1[0],c1[1],c1[2]);
        e1 = MAX3(c1[0],c1[1],c1[2]);
        s2 = MIN3(c2[0],c2[1],c2[2]);
        e2 = MAX3(c2[0],c2[1],c2[2]);

        bracket(p,d1,s1,e1,g,as,ae);
        bracket(p,d2,s2,e2,g,bs,be);

        for(int a=as; a<ae; ++a)
        for(int b=bs; b<be; ++b)
        {
            if(dir==0)
            {j=a; k=b;}
            else if(dir==1)
            {i=a; k=b;}
            else
            {i=a; j=b;}
            
            setray(p,dir,i,j,k,g,r);
            if(cross(r,tol,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,u,v,w)==0)
            continue;

            // crossings between two solid cells are no interface: no distance
            // (grid solids: only crossings in cells of the global domain are checked, as in DIVEMesh)
            int distcheck=1;

            if(dir==0)
            {
                R = u*Ax + v*Bx + w*Cx;

                i = cell(p,0,R,g);

                const bool valid = (i>=g.is+off && i<g.ie-off) && (g.interior || (i+p->origin_i>=0 && i+p->origin_i<p->gknox));

                if(R<p->XP[IP])
                if(valid)
                if(IO[IJK]<0 && IO[Im1JK]<0)
                distcheck=0;

                if(R>=p->XP[IP])
                if(valid)
                if(IO[IJK]<0 && IO[Ip1JK]<0)
                distcheck=0;

                if(distcheck==1)
                for(i=g.is; i<g.ie; ++i)
                LS[IJK]=MIN(fabs(R-p->XP[IP]),fabs(LS[IJK]));
            }
            else if(dir==1)
            {
                R = u*Ay + v*By + w*Cy;

                j = cell(p,1,R,g);

                const bool valid = (j>=g.js+off && j<g.je-off) && (g.interior || (j+p->origin_j>=0 && j+p->origin_j<p->gknoy));

                if(R<p->YP[JP])
                if(valid)
                if(IO[IJK]<0 && IO[IJm1K]<0)
                distcheck=0;

                if(R>=p->YP[JP])
                if(valid)
                if(IO[IJK]<0 && IO[IJp1K]<0)
                distcheck=0;

                if(distcheck==1)
                for(j=g.js; j<g.je; ++j)
                LS[IJK]=MIN(fabs(R-p->YP[JP]),fabs(LS[IJK]));
            }
            else
            {
                R = u*Az + v*Bz + w*Cz;

                k = cell(p,2,R,g);

                const bool valid = (k>=g.ks+off && k<g.ke-off) && (g.interior || (k+p->origin_k>=0 && k+p->origin_k<p->gknoz));

                if(R<p->ZP[KP])
                if(valid)
                if(IO[IJK]<0 && IO[IJKm1]<0)
                distcheck=0;

                if(R>=p->ZP[KP])
                if(valid)
                if(IO[IJK]<0 && IO[IJKp1]<0)
                distcheck=0;

                if(distcheck==1)
                for(k=g.ks; k<g.ke; ++k)
                LS[IJK]=MIN(fabs(R-p->ZP[KP]),LS[IJK]);

                // bed level of the column: highest crossing inside the domain
                if(bed!=nullptr && i>=0 && i<p->knox && j>=0 && j<p->knoy)
                if(R>p->global_zmin && R<p->global_zmax)
                if(R<p->global_zmax-1.0e-10)
                bed[IJ] = MAX(bed[IJ],R);
            }
        }
    }
}

void geo_raycast::cart_vertexdist(lexer *p, double **tri_x, double **tri_y, double **tri_z, int ts, int te, double *LS)
{
    int i,j,k,n;
    int is,ie,js,je,ks,ke;
    double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
    double xs,xe,ys,ye,zs,ze;
    double xc,yc,zc,dist;

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

        if(!checkin(p,1,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz))
        continue;

        xs = MIN3(Ax,Bx,Cx);
        xe = MAX3(Ax,Bx,Cx);

        ys = MIN3(Ay,By,Cy);
        ye = MAX3(Ay,By,Cy);

        zs = MIN3(Az,Bz,Cz);
        ze = MAX3(Az,Bz,Cz);

        is = p->posc_i(xs)-2;
        ie = p->posc_i(xe)+2;

        js = p->posc_j(ys)-2;
        je = p->posc_j(ye)+2;

        ks = p->posc_k(zs)-2;
        ke = p->posc_k(ze)+2;

        is = MAX(is,0);
        ie = MIN(ie,p->knox);

        js = MAX(js,0);
        je = MIN(je,p->knoy);

        ks = MAX(ks,0);
        ke = MIN(ke,p->knoz);

        for(i=is;i<ie;i++)
        for(j=js;j<je;j++)
        for(k=ks;k<ke;k++)
        {
            xc = p->XP[IP];
            yc = p->YP[JP];
            zc = p->ZP[KP];

            dist = sqrt(pow(xc-Ax,2.0) + pow(yc-Ay,2.0) + pow(zc-Az,2.0));
            LS[IJK]=MIN(dist,fabs(LS[IJK]));

            dist = sqrt(pow(xc-Bx,2.0) + pow(yc-By,2.0) + pow(zc-Bz,2.0));
            LS[IJK]=MIN(dist,fabs(LS[IJK]));

            dist = sqrt(pow(xc-Cx,2.0) + pow(yc-Cy,2.0) + pow(zc-Cz,2.0));
            LS[IJK]=MIN(dist,fabs(LS[IJK]));
        }
    }
}

void geo_raycast::bracket(lexer *p, int d, double s, double e, const geo_cart &g, int &is, int &ie)
{
    // cells of the triangle footprint along d, extended by epsi cells (6DOF kernels)
    const double *DP = (d==0) ? p->DXP : ((d==1) ? p->DYP : p->DZP);
    const int lo = (d==0) ? g.is : ((d==1) ? g.js : g.ks);
    const int hi = (d==0) ? g.ie : ((d==1) ? g.je : g.ke);
    const int kno = (d==0) ? p->knox : ((d==1) ? p->knoy : p->knoz);

    is = cell(p,d,s,g);
    ie = cell(p,d,e,g);

    // spacing lookup inside the node arrays
    const int qs = g.interior ? is : MAX(MIN(is,kno+marge-1),-marge);
    const int qe = g.interior ? ie : MAX(MIN(ie,kno+marge-1),-marge);

    s = s - epsi*DP[qs + marge];
    e = e + epsi*DP[qe + marge];

    is = cell(p,d,s,g);
    ie = cell(p,d,e,g);

    is = MAX(is,lo);
    ie = MIN(ie,hi);
}
