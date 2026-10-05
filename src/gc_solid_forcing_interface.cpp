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
Authors: Hans Bihs, Tobias Martin, Ahmet Soydan
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

double ghostcell::Hsolidface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi, H, phival_sf,dirac;
    
    //cout<<"p->topoforcing: "<<p->topoforcing<<" p->toporead: "<<p->toporead<<endl;
    
    if (p->j_dir==0)
    psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
	
    if (p->j_dir==1)
    psi = p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    // Construct solid heaviside function

    if(p->topoforcing>0 && p->solidread>0)
    phival_sf = MIN(0.5*(a->solid(i,j,k) + a->solid(i+aa,j+bb,k+cc)), 0.5*(a->topo(i,j,k) + a->topo(i+aa,j+bb,k+cc))); 
    
    if(p->topoforcing>0 && p->solidread==0)
    phival_sf = 0.5*(a->topo(i,j,k) + a->topo(i+aa,j+bb,k+cc)); 
    
    if(p->topoforcing==0 && p->solidread>0)
    phival_sf = 0.5*(a->solid(i,j,k) + a->solid(i+aa,j+bb,k+cc));
    
	
    if (-phival_sf > psi)
    H = 1.0;

    else if (-phival_sf < -psi)
    H = 0.0;

    else
    H = 0.5*(1.0 + -phival_sf/psi + (1.0/PI)*sin((PI*-phival_sf)/psi));

    // particle bed (CPM, S 10 1): the bed level is the iso-surface of the solid fraction of the parcels and
    // moves freely inside the top cell of the bed; the forcing must not reach into the fluid above it
    // (it would damp the first fluid cell whenever the bed level comes within psi of its centre):
    // one-sided transition over the depth psi inside the bed; solid bodies keep the symmetric transition
    if(p->S10==1 && p->Q10>0 && p->topoforcing>0)
    {
        double xs = -0.5*(a->topo(i,j,k) + a->topo(i+aa,j+bb,k+cc));

        if(xs<=0.0)
        H = 0.0;

        else if(xs>=psi)
        H = 1.0;

        else
        {
            double sx = 2.0*xs/psi - 1.0;
            H = 0.5*(1.0 + sx + (1.0/PI)*sin(PI*sx));
        }
        
        if(p->solidread>0)
        {
            double ps = -0.5*(a->solid(i,j,k) + a->solid(i+aa,j+bb,k+cc));
            double Hs;
            
            if(ps > psi)
            Hs = 1.0;
            
            else if(ps < -psi)
            Hs = 0.0;
            
            else
            Hs = 0.5*(1.0 + ps/psi + (1.0/PI)*sin((PI*ps)/psi));
            
            H = MAX(H,Hs);
        }
    }
    
    if(p->topoforcing==0 && p->solidread==0)
    H = 0.0;
    
    return H;
}

double ghostcell::Hsolidface_t(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi, H, phival_sf,dirac;
	
    psi = 0.5*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if (p->knoy==1)
    {
        psi = 0.5*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
    }

    // Construct solid heaviside function
    phival_sf = MIN(0.5*(a->solid(i,j,k) + a->solid(i+aa,j+bb,k+cc)), 0.5*(a->topo(i,j,k) + a->topo(i+aa,j+bb,k+cc))); 
	
    if (-phival_sf > psi)
    {
        H = 1.0;
    }
    
    else if (-phival_sf < -psi)
    {
        H = 0.0;
    }
    
    else
    {
        H = 0.5*(1.0 + -phival_sf/psi + (1.0/PI)*sin((PI*-phival_sf)/psi));
    }
	
    return H;
}



