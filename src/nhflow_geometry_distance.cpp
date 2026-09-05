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
#define WLVL (fabs(d->WL(i,j))>0.00005?d->WL(i,j):1.0e20)

// squared distance from P to triangle ABC
double nhflow_geometry::dist2_tri(const double Px,const double Py,const double Pz,
                                  const double Ax,const double Ay,const double Az,
                                  const double Bx,const double By,const double Bz,
                                  const double Cx,const double Cy,const double Cz)
{
    double abx=Bx-Ax, aby=By-Ay, abz=Bz-Az;
    double acx=Cx-Ax, acy=Cy-Ay, acz=Cz-Az;
    double apx=Px-Ax, apy=Py-Ay, apz=Pz-Az;

    double qx,qy,qz;   // closest point

    double d1 = abx*apx + aby*apy + abz*apz;
    double d2 = acx*apx + acy*apy + acz*apz;

    if(d1<=0.0 && d2<=0.0)                       // vertex A
    {qx=Ax; qy=Ay; qz=Az; goto done;}

    {
    double bpx=Px-Bx, bpy=Py-By, bpz=Pz-Bz;
    double d3 = abx*bpx + aby*bpy + abz*bpz;
    double d4 = acx*bpx + acy*bpy + acz*bpz;

    if(d3>=0.0 && d4<=d3)                        // vertex B
    {qx=Bx; qy=By; qz=Bz; goto done;}

    double vc = d1*d4 - d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0)            // edge AB
    {
        double v = d1/((d1-d3)!=0.0 ? (d1-d3) : 1.0e-20);
        qx=Ax+v*abx; qy=Ay+v*aby; qz=Az+v*abz; goto done;
    }

    double cpx=Px-Cx, cpy=Py-Cy, cpz=Pz-Cz;
    double d5 = abx*cpx + aby*cpy + abz*cpz;
    double d6 = acx*cpx + acy*cpy + acz*cpz;

    if(d6>=0.0 && d5<=d6)                        // vertex C
    {qx=Cx; qy=Cy; qz=Cz; goto done;}

    double vb = d5*d2 - d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0)            // edge AC
    {
        double w = d2/((d2-d6)!=0.0 ? (d2-d6) : 1.0e-20);
        qx=Ax+w*acx; qy=Ay+w*acy; qz=Az+w*acz; goto done;
    }

    double va = d3*d6 - d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0)  // edge BC
    {
        double den = (d4-d3)+(d5-d6);
        double w = (d4-d3)/(den!=0.0 ? den : 1.0e-20);
        qx=Bx+w*(Cx-Bx); qy=By+w*(Cy-By); qz=Bz+w*(Cz-Bz); goto done;
    }

    {                                            // face interior
        double den = va+vb+vc;
        den = 1.0/(fabs(den)>1.0e-20 ? den : 1.0e-20);
        double v = vb*den;
        double w = vc*den;
        qx=Ax+abx*v+acx*w; qy=Ay+aby*v+acy*w; qz=Az+abz*v+acz*w;
    }
    }

done:
    return (Px-qx)*(Px-qx) + (Py-qy)*(Py-qy) + (Pz-qz)*(Pz-qz);
}


void nhflow_geometry::band_distance(lexer *p, fdm_nhf *d, ghostcell *pgc, double *LS, int ts, int te)
{
    const double band = NB*DSM;      // NB ~ 4 cells; must exceed A526 by 2-3
    const double band2 = band*band;

    for(n=ts; n<te; ++n)
    {
        double Ax=tri_x[n][0], Ay=tri_y[n][0], Az=tri_z[n][0];
        double Bx=tri_x[n][1], By=tri_y[n][1], Bz=tri_z[n][1];
        double Cx=tri_x[n][2], Cy=tri_y[n][2], Cz=tri_z[n][2];
        
        // <<<< INSERT HERE >>>>
        // 2D: skip the y-normal end caps. In a one-cell-wide domain they sit
        // DYN/2 from every interior cell and win every distance comparison,
        // flattening the interior to a constant. Matches ray_cast_z and
        // force_calc_stl, which both drop |ny|>0.9 facets when j_dir==0.
        if(p->j_dir==0)
        {
            double enx = (By-Ay)*(Cz-Az) - (Cy-Ay)*(Bz-Az);
            double eny = (Cx-Ax)*(Bz-Az) - (Bx-Ax)*(Cz-Az);
            double enz = (Bx-Ax)*(Cy-Ay) - (Cx-Ax)*(By-Ay);
            double enm = sqrt(enx*enx + eny*eny + enz*enz);

            if(enm > 1.0e-20 && fabs(eny)/enm > 0.9)
            continue;
        }

        double txs=MIN3(Ax,Bx,Cx)-band, txe=MAX3(Ax,Bx,Cx)+band;
        double tys=MIN3(Ay,By,Cy)-band, tye=MAX3(Ay,By,Cy)+band;
        double tzs=MIN3(Az,Bz,Cz)-band, tze=MAX3(Az,Bz,Cz)+band;

        // cheap reject: triangle band outside this subdomain
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
                // AABB pre-test in physical space, avoids the full kernel
                if(p->XP[IP] < txs || p->XP[IP] > txe) continue;
                if(p->j_dir==1 && (p->YP[JP] < tys || p->YP[JP] > tye)) continue;
                if(p->ZSP[IJK] < tzs || p->ZSP[IJK] > tze) continue;

                double dd = dist2_tri(p->XP[IP], p->YP[JP], p->ZSP[IJK],
                                      Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz);

                if(dd < band2 && dd < LS[IJK]*LS[IJK])
                LS[IJK] = sqrt(dd);
            }
        }
    }
}