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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"met_forcing.h"
#include"met_provider.h"
#include"met_frame.h"
#include"met_drag.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<algorithm>
#include<iostream>

met_forcing::met_forcing(lexer *p, ghostcell *pgc) : U10(p), V10(p), Pa(p), Tx(p), Ty(p), Gx(p), Gy(p), Fc(p)
{
    wind = p->W200>0;
    pressure = wind && p->W208==1;
    coriolis = p->W211>0;

    if(wind)
    {
    pm = met_provider::get(p,pgc);

    // box of the cell centres of the grid incl. its ghost cells: the part of the met file the rank reads
    double xs = 1.0e300, xe = -1.0e300, ys = 1.0e300, ye = -1.0e300;

        IMALOOP
        {
        xs = std::min(xs,p->XP[IP]);
        xe = std::max(xe,p->XP[IP]);
        }

        JMALOOP
        {
        ys = std::min(ys,p->YP[JP]);
        ye = std::max(ye,p->YP[JP]);
        }

    pm->window(xs,xe,ys,ye);
    }

    // Coriolis parameter f = 2 Omega sin(latitude)
    if(coriolis)
    {
    const double Omega = 7.2921159e-5;
    const double rad = 3.14159265358979323846/180.0;
    met_frame fr;

        if(p->W205==1)
        fr = met_frame(p->W205_lon,p->W205_lat,p->W205_rot,p->W206_x,p->W206_y);

    const double f0 = p->W212!=0.0 ? p->W212 : 2.0*Omega*std::sin(p->W205_lat*rad);

        IMALOOP
        JMALOOP
        {
            if(p->W211==1)
            Fc(i,j) = f0;

            if(p->W211==2)
            {
            double lon, lat;
            // curvilinear grid (G 50): the physical cell centre (XP, YP are the logical coordinates there)
            if(p->geo_curv==1)
            fr.to_geo(p->XC2D[IJ],p->YC2D[IJ],lon,lat);
            else
            fr.to_geo(p->XP[IP],p->YP[JP],lon,lat);
            Fc(i,j) = 2.0*Omega*std::sin(lat*rad);
            }
        }
    }
}

void met_forcing::update(lexer *p, ghostcell *pgc)
{
    if(!wind || (has_t && p->simtime==tlast))
    return;

    pm->update(p,pgc,p->simtime);

    const double r = pm->ramp(p,p->simtime);
    const double ra = r*p->W3/p->W1;
    const int law = met_drag_law(p->W207,p->W200);
    double u, v, pa, speed;

    IMALOOP
    JMALOOP
    {
    pm->at(p->XP[IP],p->YP[JP],u,v,pa);

    speed = std::sqrt(u*u + v*v);

    U10(i,j) = u;
    V10(i,j) = v;
    Pa(i,j)  = pa;

    Tx(i,j) = ra*met_drag(law,speed,p->W213,p->W214)*speed*u;
    Ty(i,j) = ra*met_drag(law,speed,p->W213,p->W214)*speed*v;
    }

    // kinematic air pressure gradient, central differences (the ghost cells hold the source's values)
    if(pressure)
    {
    const double rp = r/p->W1;

        SLICELOOP4
        {
        Gx(i,j) = rp*(Pa(i+1,j) - Pa(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
        Gy(i,j) = p->j_dir==1 ? rp*(Pa(i,j+1) - Pa(i,j-1))/(p->DYP[JP] + p->DYP[JM1]) : 0.0;
        }
    }

    tlast = p->simtime;
    has_t = true;
}

void met_forcing::source_sflow(lexer *p, ghostcell *pgc, slice &F, int dir, slice &WL, slice &UH, slice &VH)
{
    update(p,pgc);

    if(dir==1 && p->j_dir==0)
    return;

    slice4 &T = (dir==0) ? Tx : Ty;
    slice4 &G = (dir==0) ? Gx : Gy;

    if(wind)
    SLICELOOP4
    WETDRY
    F(i,j) += T(i,j);

    if(pressure)
    SLICELOOP4
    WETDRY
    F(i,j) -= WL(i,j)*G(i,j);

    if(coriolis)
    {
        if(dir==0)
        SLICELOOP4
        WETDRY
        F(i,j) += Fc(i,j)*VH(i,j);

        if(dir==1)
        SLICELOOP4
        WETDRY
        F(i,j) -= Fc(i,j)*UH(i,j);
    }
}

void met_forcing::source_nhflow(lexer *p, ghostcell *pgc, double *F, int dir, slice &WL, double *UH, double *VH)
{
    update(p,pgc);

    if(dir==1 && p->j_dir==0)
    return;

    slice4 &T = (dir==0) ? Tx : Ty;
    slice4 &G = (dir==0) ? Gx : Gy;

    // surface stress in the top layer
    if(wind)
    {
    k = p->knoz-1;

        SLICELOOP4
        WETDRY
        F[IJK] += T(i,j)/p->DZN[KP];
    }

    // air pressure gradient and Coriolis force in every layer
    if(pressure)
    LOOP
    WETDRY
    F[IJK] -= WL(i,j)*G(i,j);

    if(coriolis)
    {
        if(dir==0)
        LOOP
        WETDRY
        F[IJK] += Fc(i,j)*VH[IJK];

        if(dir==1)
        LOOP
        WETDRY
        F[IJK] -= Fc(i,j)*UH[IJK];
    }
}

void met_forcing::offset_vtp(lexer *p, int *offset, int &n)
{
    if(!wind)
    return;

    // wind10 (3 components), air_pressure, wind_stress (3 components)
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)*3+sizeof(int);
    ++n;
}

void met_forcing::name_vtp(ofstream &result, int *offset, int &n)
{
    if(!wind)
    return;

    result<<"<DataArray type=\"Float32\" Name=\"wind10\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"air_pressure\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"wind_stress\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void met_forcing::name_pvtp(ofstream &result)
{
    if(!wind)
    return;

    result<<"<PDataArray type=\"Float32\" Name=\"wind10\" NumberOfComponents=\"3\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"air_pressure\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"wind_stress\" NumberOfComponents=\"3\"/>\n";
}

void met_forcing::print_vtp(lexer *p, ofstream &result)
{
    if(!wind)
    return;

    int iin;
    float ffn;

    // wind10
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(U10));
    result.write((char*)&ffn, sizeof(float));

    ffn=float(p->sl_ipol4(V10));
    result.write((char*)&ffn, sizeof(float));

    ffn=0.0;
    result.write((char*)&ffn, sizeof(float));
    }

    // air_pressure
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(Pa));
    result.write((char*)&ffn, sizeof(float));
    }

    // wind_stress [N/m^2]
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->W1*p->sl_ipol4(Tx));
    result.write((char*)&ffn, sizeof(float));

    ffn=float(p->W1*p->sl_ipol4(Ty));
    result.write((char*)&ffn, sizeof(float));

    ffn=0.0;
    result.write((char*)&ffn, sizeof(float));
    }
}
