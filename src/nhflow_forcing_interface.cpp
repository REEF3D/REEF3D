/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_reinidisc_fsf.h"

double nhflow_forcing::Hsolidface(lexer *p, fdm_nhf *d, int aa, int bb, int cc)
{
    double psi, H, phival_sf,dirac;
    
    if (p->j_dir==0)
    psi = p->X41*(1.0/1.0)*(p->DXN[IP]);
	
    if (p->j_dir==1)
    psi = p->X41*(1.0/2.0)*(p->DXN[IP]+p->DYN[JP]);


    // Construct solid heaviside function
    phival_sf = d->SOLID[IJK];
    
	
    if(-phival_sf > psi)
    H = 1.0;

    if(-phival_sf < -psi)
    H = 0.0;

    if(fabs(phival_sf)<=psi)
    H = 0.5*(1.0 + -phival_sf/psi + (1.0/PI)*sin((PI*-phival_sf)/psi));
    
    /*
    H = 0.0;
    if(fabs(phival_sf)<psi)
    H = (0.5/psi)*(1.0 + cos((PI*(phival_sf))/psi));

    
    H=MIN(H,1.0);*/
    
    
    /*
    H = 0.0;
    if(fabs(phival_sf)<10.0*psi)
    H = psi/(PI*(phival_sf*phival_sf + psi*psi));
    */
    
    //H=MIN(H,1.0);
    

    return H;
}