/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres::seed_particles(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    double spacing;
    double particles_per_dim;
    
    // estimate number particles
    int count=0;
    ALOOP
    if(a->topo(i,j,k)<=0)
        ++count;

    // safety
    count += 100;

    count *= p->Q24;

    P.resize(p,count);
    
    const int irand(100000);
    const double drand(100000.0);
    
    // seed equal spacing
    if(p->Q29==1)
    {
        n=0;
        ALOOP
        if(a->topo(i,j,k)<=0)
        {
        int particles_per_cell = p->Q24;  
        int particles_per_dim = int(pow(double(particles_per_cell), 1.0/3.0) + 0.5);
        
        spacing = (1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]) / double(particles_per_dim);

            for(int ii = 0; ii < particles_per_dim; ++ii)
            for(int jj = 0; jj < particles_per_dim; ++jj)
            for(int kk = 0; kk < particles_per_dim; ++kk)
            {
                n=P.Empty[P.index_empty];
                P.X[n] = p->XN[IP] + double(ii) * spacing ;
                P.Y[n] = p->YN[JP] + double(jj) * spacing ;
                P.Z[n] = p->ZN[KP] + double(kk) * spacing ;

                P.D[n] = p->S20;
                P.RO[n] = p->S22;

                P.Flag[n] = ACTIVE;
                --P.index_empty;
            }
        }
    }
    
    // seed equal spacing with jitter
    if(p->Q29==2)
    {
        n=0;
        ALOOP
        if(a->topo(i,j,k)<=0)
        {
        int particles_per_cell = p->Q24; 
        int particles_per_dim = int(pow(double(particles_per_cell), 1.0/3.0) + 0.5);        
        
        spacing = (1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]) / double(particles_per_dim);
        
        //cout<<"particles_per_dim: "<<particles_per_dim<<" spacing: "<<spacing<<endl;

            for(int ii = 0; ii < particles_per_dim; ++ii)
            for(int jj = 0; jj < particles_per_dim; ++jj)
            for(int kk = 0; kk < particles_per_dim; ++kk)
            {
                n=P.Empty[P.index_empty];
                P.X[n] = p->XN[IP] + double(ii) * spacing + 0.1*p->DXN[IP]*double(rand() % irand)/drand;
                P.Y[n] = p->YN[JP] + double(jj) * spacing + 0.1*p->DYN[JP]*double(rand() % irand)/drand;
                P.Z[n] = p->ZN[KP] + double(kk) * spacing + 0.1*p->DZN[KP]*double(rand() % irand)/drand;

                P.D[n] = p->S20;
                P.RO[n] = p->S22;

                P.Flag[n] = ACTIVE;
                --P.index_empty;
            }
        }
    }
    
    // seed fully random
    if(p->Q29==3)
    {
        n=0;
        ALOOP
        if(a->topo(i,j,k)<=0)
        {
            for(int qn=0;qn<p->Q24;++qn)
            {
                n=P.Empty[P.index_empty];
                P.X[n] = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                P.Y[n] = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                P.Z[n] = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;

                P.D[n] = p->S20;
                P.RO[n] = p->S22;

                P.Flag[n] = ACTIVE;
                --P.index_empty;
            }
        }
    }

    // remove above bed
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        double topoval  = p->ccipol4_b(a->topo,P.X[n],P.Y[n],P.Z[n]);
        double solidval = p->ccipol4_b(a->solid,P.X[n],P.Y[n],P.Z[n]);

        if(topoval>0.0)
            P.remove(n);

        if(solidval<0.0)
            P.remove(n);
    }

    cellSum_update(p,pgc,s,P.X,P.Y,P.Z);
    
    wallbc(p,pgc,s);
}
