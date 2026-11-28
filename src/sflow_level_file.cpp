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


#include"sflow_level_file.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<fstream>
#include<iomanip>
#include<sys/stat.h>
#include<vector>
#include<cmath>

void level_grid::read(lexer *p, ghostcell *pgc, const char *file, const char *key)
{
    std::ifstream f(file);
    
    if(!f)
    {
    if(p->mpirank==0)
    cout<<key<<": "<<file<<" not found"<<endl;
    pgc->final(EXIT_FAILURE);
    }
    
    f>>nx>>ny>>x0>>y0>>dx>>dy;
    
    if(nx>0 && ny>0)
    L.assign((size_t)nx*ny,-9999.0);
    
    for(size_t n=0; n<L.size(); ++n)
    f>>L[n];
    
    if(!f || nx<2 || ny<2 || dx<=0.0 || dy<=0.0)
    {
    if(p->mpirank==0)
    cout<<key<<": "<<file<<" could not be read (header nx ny x0 y0 dx dy, then nx*ny levels)"<<endl;
    pgc->final(EXIT_FAILURE);
    }
}

double level_grid::value(double x, double y) const
{
    const double fi = (x - x0)/dx, fj = (y - y0)/dy;
    
    if(fi<-1.0 || fj<-1.0 || fi>nx || fj>ny)
    return -9999.0;
    
    auto val = [&](int ii, int jj) { return L[(size_t)jj*nx + ii]; };
    
    const int i0 = std::max(0,std::min(nx-2,(int)floor(fi)));
    const int j0 = std::max(0,std::min(ny-2,(int)floor(fj)));
    const double s = std::max(0.0,std::min(1.0,fi-i0)), t = std::max(0.0,std::min(1.0,fj-j0));
    
    const double v00=val(i0,j0), v10=val(i0+1,j0), v01=val(i0,j0+1), v11=val(i0+1,j0+1);
    
    if(v00>-9000.0 && v10>-9000.0 && v01>-9000.0 && v11>-9000.0)
    return (1.0-s)*(1.0-t)*v00 + s*(1.0-t)*v10 + (1.0-s)*t*v01 + s*t*v11;
    
    const int in = std::max(0,std::min(nx-1,(int)lround(fi)));
    const int jn = std::max(0,std::min(ny-1,(int)lround(fj)));
    return val(in,jn);
}

double level_grid::nearest(double x, double y) const
{
    const double fi = (x - x0)/dx, fj = (y - y0)/dy;
    
    if(fi<-1.0 || fj<-1.0 || fi>nx || fj>ny)
    return -9999.0;
    
    const int in = std::max(0,std::min(nx-1,(int)lround(fi)));
    const int jn = std::max(0,std::min(ny-1,(int)lround(fj)));
    return L[(size_t)jn*nx + in];
}

void sflow_level_file::start(lexer *p, fdm2D *b, ghostcell *pgc, sflow_release *rel)
{
    level_grid G;
    G.read(p,pgc,"sflow-level.dat","F 73");
    
    int count = 0;
    double vol = 0.0;
    
    SLICELOOP4
    {
        const double x = p->geo_curv==1 ? p->XC2D[IJ] : p->XP[IP];
        const double y = p->geo_curv==1 ? p->YC2D[IJ] : p->YP[JP];
        
        const double level = G.value(x,y);
        
        if(level<-9000.0 || level<=b->bed(i,j))
        continue;
        
        // F 74: released later, part of the bed until then
        if(rel!=nullptr && rel->time(x,y)>0.0)
        {
        rel->lock(p,b,i,j,level);
        continue;
        }
        
        b->eta(i,j) = level - p->wd;
        ++count;
        vol += (level - b->bed(i,j))*((p->geo_curv==1) ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP]);
    }
    
    count = pgc->globalisum(count);
    vol = pgc->globalsum(vol);
    
    if(p->mpirank==0)
    cout<<"F 73: initial level from sflow-level.dat ("<<G.nx<<" x "<<G.ny<<", spacing "<<G.dx<<" x "<<G.dy<<" m) in "<<count<<" cells, volume above the bed "<<vol<<" m3"<<(rel!=nullptr ? " (released at the start)" : "")<<endl;
    
    if(rel!=nullptr)
    {
    pgc->gcsl_start4(p,b->bed,50);
    rel->report_ini(p,pgc);
    }
}

// ---------------------------------------------------------------------------------------------
// F 74: staged release
// ---------------------------------------------------------------------------------------------

sflow_release::sflow_release(lexer *p, ghostcell *pgc) : tr(p), thick(p)
{
    R.read(p,pgc,"sflow-release.dat","F 74");
    
    SLICELOOP4
    {
    tr(i,j) = 0.0;
    thick(i,j) = 0.0;
    }
}

double sflow_release::time(double x, double y) const
{
    const double t = R.nearest(x,y);
    return t<-9000.0 ? 0.0 : t;
}

void sflow_release::lock(lexer *p, fdm2D *b, int i, int j, double level)
{
    const double x = p->geo_curv==1 ? p->XC2D[IJ] : p->XP[IP];
    const double y = p->geo_curv==1 ? p->YC2D[IJ] : p->YP[JP];
    
    tr(i,j) = time(x,y);
    thick(i,j) = level - b->bed(i,j);
    b->bed(i,j) = level;
}

void sflow_release::report_ini(lexer *p, ghostcell *pgc)
{
    int count=0;
    double vol=0.0, tmin=1.0e20, tmax=-1.0e20;
    
    SLICELOOP4
    if(thick(i,j)>0.0)
    {
    ++count;
    vol += thick(i,j)*((p->geo_curv==1) ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP]);
    tmin = MIN(tmin,tr(i,j));
    tmax = MAX(tmax,tr(i,j));
    }
    
    count = pgc->globalisum(count);
    vol = pgc->globalsum(vol);
    tmin = pgc->globalmin(tmin);
    tmax = pgc->globalmax(tmax);
    
    if(p->mpirank==0)
    cout<<"F 74: staged release from sflow-release.dat: "<<count<<" cells, "<<vol<<" m3 part of the bed until their time, released between t = "<<tmin<<" and "<<tmax<<" s"<<endl;
}

void sflow_release::step(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the cells whose time falls within this step: at the end of the step
    const double tend = p->simtime + p->dt;
    
    int count=0;
    double vol=0.0;
    
    SLICELOOP4
    if(thick(i,j)>0.0 && tr(i,j)<=tend+1.0e-12)
    {
        const double T = thick(i,j);
        const bool dry = b->WL(i,j)<p->A244;
        
        b->bed(i,j) -= T;
        b->WL(i,j) += T;
        thick(i,j) = 0.0;
        
        // wet now: the wetting of SFLOW (A 243 2) wets a dry cell only from a higher wet neighbour
        p->wet[IJ] = 1;
        
        // a dry cell starts at rest
        if(dry)
        {
        b->P(i,j) = b->P(i-1,j) = 0.0;
        b->Q(i,j) = b->Q(i,j-1) = 0.0;
        }
        
        ++count;
        vol += T*((p->geo_curv==1) ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP]);
    }
    
    count = pgc->globalisum(count);
    
    if(count>0)
    {
    vol = pgc->globalsum(vol);
    
    pgc->gcsl_start4(p,b->bed,50);
    pgc->gcsl_start4(p,b->WL,50);
    pgc->gcsl_start4Vint(p,p->wet,50);
    pgc->gcsl_start1(p,b->P,10);
    pgc->gcsl_start2(p,b->Q,11);
    
    ++nevents;
    
    if(p->mpirank==0)
    cout<<"F 74: released "<<count<<" cells, "<<vol<<" m3 at t = "<<tend<<" s"<<endl;
    }
    
    // the log: the volume of the flow and the volume still in the bed (their sum is the volume of F 73)
    double vflow=0.0, vbed=0.0;
    
    SLICELOOP4
    {
    const double A = (p->geo_curv==1) ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP];
    vflow += b->WL(i,j)*A;
    vbed += thick(i,j)*A;
    }
    
    vflow = pgc->globalsum(vflow);
    vbed = pgc->globalsum(vbed);
    
    if(p->mpirank==0)
    {
        if(!logout.is_open())
        {
        mkdir("./REEF3D_SFLOW_Log",0777);
        logout.open("./REEF3D_SFLOW_Log/REEF3D_SFLOW_release.dat");
        logout<<"# time  V_flow  V_bed (not yet released)  V_flow + V_bed  (m3; F 74)"<<endl;
        }
        
        logout<<setprecision(8)<<tend<<" "<<vflow<<" "<<vbed<<" "<<vflow+vbed<<endl;
    }
}
