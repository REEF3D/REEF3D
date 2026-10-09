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

// Porous layers (B 202): depth below the exposed surface of the porous structures.
//
// The layer depth is the exact point-triangle distance to the exposed surface
// triangles, evaluated in a band of the total layer thickness. A triangle is
// exposed when the point just outside it (outward normal) lies
//   - inside the global domain (faces on or beyond the domain boundary are skipped),
//   - outside all other porous entities (internal faces of the union are skipped),
//   - and, with B 203 1, not below the face: faces whose outward normal points
//     more than 45 deg downwards rest on the bed and are skipped.
// In 2D the y-end caps are skipped, as in geo_raycast::sigma_band.
//
// The triangles are static, so the classification is done once (layer_faces).
// The sigma nodes move with the free surface, so the distance is recomputed
// with the level set in every update (layer_dist).

// inside/outside of entity qn by vertical ray parity, same test and ray offset as
// geo_raycast::sigma_io, incl. raymode and inversion of the entity
bool nhflow_geometry::point_inside(lexer *p, int qn, double Px, double Py, double Pz)
{
    double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
    double abx,aby,bcx,bcy,cax,cay;
    double e0,e1,e2,area2,sgn,denom,Rz;
    int above=0, below=0;

    const double eps_area = 1.0e-20*DSM*DSM;
    const double psi      = 1.0e-8*DSM;

    Px -= psi;
    Py += psi;

    for(int n=tstart[qn]; n<tend[qn]; ++n)
    {
    Ax = tri_x[n][0];  Ay = tri_y[n][0];  Az = tri_z[n][0];
    Bx = tri_x[n][1];  By = tri_y[n][1];  Bz = tri_z[n][1];
    Cx = tri_x[n][2];  Cy = tri_y[n][2];  Cz = tri_z[n][2];

    // footprint reject
    if(Px < MIN3(Ax,Bx,Cx) || Px > MAX3(Ax,Bx,Cx))
    continue;

    if(Py < MIN3(Ay,By,Cy) || Py > MAX3(Ay,By,Cy))
    continue;

    abx = Bx-Ax;  aby = By-Ay;
    bcx = Cx-Bx;  bcy = Cy-By;
    cax = Ax-Cx;  cay = Ay-Cy;

    area2 = -abx*cay + aby*cax;

    // vertical triangle: not crossed by a vertical ray
    if(fabs(area2) < eps_area)
    continue;

    sgn   = (area2 > 0.0 ? 1.0 : -1.0);
    denom = 1.0/(sgn*area2);

    e0 = sgn*(abx*(Py-Ay) - aby*(Px-Ax));
    e1 = sgn*(bcx*(Py-By) - bcy*(Px-Bx));
    e2 = sgn*(cax*(Py-Cy) - cay*(Px-Cx));

    if(e0 < 0.0 || e1 < 0.0 || e2 < 0.0)
    continue;

    Rz = (e1*Az + e2*Bz + e0*Cz)*denom;

    if(Pz < Rz)
    ++above;

    else
    ++below;
    }

    const int mode = (qn<int(ent_raymode.size())) ? ent_raymode[qn] : 1;

    bool inside;

    if(mode==2)
    inside = (above%2==0 && below%2==0);

    else
    inside = (above%2==1 && below%2==1);

    if(qn<int(ent_invert.size()) && ent_invert[qn]==1)
    inside = !inside;

    return inside;
}

void nhflow_geometry::layer_faces(lexer *p, ghostcell *pgc)
{
    int qn,qe,n;
    int exposed=0, ambiguous=0;

    tri_layer.assign(tricount,0);

    // largest horizontal cell size: the band must cover the Heaviside width B209*dx
    layer_dsmax = 0.0;

    SLICELOOP4
    {
    layer_dsmax = MAX(layer_dsmax, p->DXN[IP]);

    if(p->j_dir==1)
    layer_dsmax = MAX(layer_dsmax, p->DYN[JP]);
    }

    layer_dsmax = pgc->globalmax(layer_dsmax);

    // offset of the test points from the face
    const double h = 1.0e-3*DSM;

    for(qn=0; qn<entity_sum; ++qn)
    for(n=tstart[qn]; n<tend[qn]; ++n)
    {
        double Ax=tri_x[n][0], Ay=tri_y[n][0], Az=tri_z[n][0];
        double Bx=tri_x[n][1], By=tri_y[n][1], Bz=tri_z[n][1];
        double Cx=tri_x[n][2], Cy=tri_y[n][2], Cz=tri_z[n][2];

        double nx = (By-Ay)*(Cz-Az) - (Bz-Az)*(Cy-Ay);
        double ny = (Bz-Az)*(Cx-Ax) - (Bx-Ax)*(Cz-Az);
        double nz = (Bx-Ax)*(Cy-Ay) - (By-Ay)*(Cx-Ax);
        double nm = sqrt(nx*nx + ny*ny + nz*nz);

        // degenerate triangle
        if(nm < 1.0e-20)
        continue;

        nx/=nm;  ny/=nm;  nz/=nm;

        // 2D: y-end caps
        if(p->j_dir==0 && fabs(ny)>0.9)
        continue;

        double xc = (Ax+Bx+Cx)/3.0;
        double yc = (Ay+By+Cy)/3.0;
        double zc = (Az+Bz+Cz)/3.0;

        bool in_p = point_inside(p,qn, xc+h*nx, yc+h*ny, zc+h*nz);
        bool in_m = point_inside(p,qn, xc-h*nx, yc-h*ny, zc-h*nz);

        // orientation undecided (open or degenerate surface): keep the face
        if(in_p==in_m)
        {
        tri_layer[n]=1;
        ++exposed;
        ++ambiguous;
        continue;
        }

        // outward normal
        double ox = in_m ? nx : -nx;
        double oy = in_m ? ny : -ny;
        double oz = in_m ? nz : -nz;

        double qx = xc + h*ox;
        double qy = yc + h*oy;
        double qz = zc + h*oz;

        // face on or beyond the domain boundary
        if(qx < p->global_xmin || qx > p->global_xmax)
        continue;

        if(p->j_dir==1 && (qy < p->global_ymin || qy > p->global_ymax))
        continue;

        // face resting on the bed
        if(p->B203==1 && oz < -0.7071)
        continue;

        // internal face of the union of the porous entities
        bool internal=false;

        for(qe=0; qe<entity_sum; ++qe)
        if(qe!=qn && point_inside(p,qe,qx,qy,qz))
        {
        internal=true;
        break;
        }

        if(internal)
        continue;

        tri_layer[n]=1;
        ++exposed;
    }

    if(p->mpirank==0)
    {
    cout<<"VRANS porous layers: "<<p->B202<<" layers, "<<exposed<<" of "<<tricount<<" surface triangles define the layer depth"<<endl;

    if(ambiguous>0)
    cout<<"VRANS porous layers: "<<ambiguous<<" triangles with undecided orientation are kept, check the STL for open edges"<<endl;

    if(exposed==0)
    cout<<"VRANS porous layers: !!! no exposed surface triangles, all porous cells get the core properties B 201 !!!"<<endl;
    }
}

// LD = distance to the exposed surface, capped at the band t_total + (2 B209 + 1) dx_max,
// evaluated where LS < band (inside the structures and in the Heaviside band around them)
void nhflow_geometry::layer_dist(lexer *p, fdm_nhf *d, double t_total, const double *LS, double *LD)
{
    int i,j,k,n;

    const double band  = t_total + (2.0*p->B209 + 1.0)*layer_dsmax;
    const double band2 = band*band;

    LOOP
    LD[IJK] = band;

    for(n=0; n<tricount; ++n)
    {
        if(tri_layer[n]==0)
        continue;

        double Ax=tri_x[n][0], Ay=tri_y[n][0], Az=tri_z[n][0];
        double Bx=tri_x[n][1], By=tri_y[n][1], Bz=tri_z[n][1];
        double Cx=tri_x[n][2], Cy=tri_y[n][2], Cz=tri_z[n][2];

        double txs=MIN3(Ax,Bx,Cx)-band, txe=MAX3(Ax,Bx,Cx)+band;
        double tys=MIN3(Ay,By,Cy)-band, tye=MAX3(Ay,By,Cy)+band;
        double tzs=MIN3(Az,Bz,Cz)-band, tze=MAX3(Az,Bz,Cz)+band;

        if(txe<p->originx || txs>p->endx) continue;
        if(p->j_dir==1 && (tye<p->originy || tys>p->endy)) continue;

        int is = MAX(p->posc_i(txs)-1, 0);
        int ie = MIN(p->posc_i(txe)+2, p->knox);
        int js = 0, je = 1;

        if(p->j_dir==1)
        {
            js = MAX(p->posc_j(tys)-1, 0);
            je = MIN(p->posc_j(tye)+2, p->knoy);
        }

        for(i=is; i<ie; ++i)
        for(j=js; j<je; ++j)
        {
            int ks = MAX(p->posc_sig(i,j,tzs)-1, 0);
            int ke = MIN(p->posc_sig(i,j,tze)+2, p->knoz);

            for(k=ks; k<ke; ++k)
            {
                if(LS[IJK] >= band) continue;

                if(p->XP[IP] < txs || p->XP[IP] > txe) continue;
                if(p->j_dir==1 && (p->YP[JP] < tys || p->YP[JP] > tye)) continue;
                if(p->ZSP[IJK] < tzs || p->ZSP[IJK] > tze) continue;

                double dd = geo_raycast::dist2_tri(p->XP[IP], p->YP[JP], p->ZSP[IJK],
                                                   Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz);

                if(dd < band2 && dd < LD[IJK]*LD[IJK])
                LD[IJK] = sqrt(dd);
            }
        }
    }
}
