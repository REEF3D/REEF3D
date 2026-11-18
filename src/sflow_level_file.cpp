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
#include<vector>
#include<cmath>

void sflow_level_file::start(lexer *p, fdm2D *b, ghostcell *pgc)
{
    std::ifstream f("sflow-level.dat");
    
    if(!f)
    {
    if(p->mpirank==0)
    cout<<"F 73: sflow-level.dat not found"<<endl;
    pgc->final(EXIT_FAILURE);
    }
    
    int nx=0, ny=0;
    double x0=0.0, y0=0.0, dx=1.0, dy=1.0;
    f>>nx>>ny>>x0>>y0>>dx>>dy;
    
    std::vector<double> L((size_t)nx*ny,-9999.0);
    
    for(size_t n=0; n<L.size(); ++n)
    f>>L[n];
    
    if(!f || nx<2 || ny<2 || dx<=0.0 || dy<=0.0)
    {
    if(p->mpirank==0)
    cout<<"F 73: sflow-level.dat could not be read (header nx ny x0 y0 dx dy, then nx*ny levels)"<<endl;
    pgc->final(EXIT_FAILURE);
    }
    
    auto val = [&](int ii, int jj) { return L[(size_t)jj*nx + ii]; };
    
    int count = 0;
    double vol = 0.0;
    
    SLICELOOP4
    {
        const double x = p->geo_curv==1 ? p->XC2D[IJ] : p->XP[IP];
        const double y = p->geo_curv==1 ? p->YC2D[IJ] : p->YP[JP];
        
        const double fi = (x - x0)/dx, fj = (y - y0)/dy;
        
        if(fi<-1.0 || fj<-1.0 || fi>nx || fj>ny)
        continue;
        
        const int i0 = std::max(0,std::min(nx-2,(int)floor(fi)));
        const int j0 = std::max(0,std::min(ny-2,(int)floor(fj)));
        const double s = std::max(0.0,std::min(1.0,fi-i0)), t = std::max(0.0,std::min(1.0,fj-j0));
        
        const double v00=val(i0,j0), v10=val(i0+1,j0), v01=val(i0,j0+1), v11=val(i0+1,j0+1);
        double level;
        
        if(v00>-9000.0 && v10>-9000.0 && v01>-9000.0 && v11>-9000.0)
        level = (1.0-s)*(1.0-t)*v00 + s*(1.0-t)*v10 + (1.0-s)*t*v01 + s*t*v11;
        
        else
        {
            const int in = std::max(0,std::min(nx-1,(int)lround(fi)));
            const int jn = std::max(0,std::min(ny-1,(int)lround(fj)));
            level = val(in,jn);
        }
        
        if(level<-9000.0 || level<=b->bed(i,j))
        continue;
        
        b->eta(i,j) = level - p->wd;
        ++count;
        vol += (level - b->bed(i,j))*((p->geo_curv==1) ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP]);
    }
    
    count = pgc->globalisum(count);
    vol = pgc->globalsum(vol);
    
    if(p->mpirank==0)
    cout<<"F 73: initial level from sflow-level.dat ("<<nx<<" x "<<ny<<", spacing "<<dx<<" x "<<dy<<" m) in "<<count<<" cells, volume above the bed "<<vol<<" m3"<<endl;
}
