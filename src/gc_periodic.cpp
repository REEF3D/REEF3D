/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"field.h"
#include"vec.h"
#include"fdm.h"

void ghostcell::gc_periodic(lexer *p, field& f, int gcv, int cs)
{
    double val1,val2,val3;
    
    if(cs==1)
    JLOOP
    KLOOP
    PCHECK
    {
    // 4 to 1 coupling
    i=p->knox-1;
    
    val1 = f(i,j,k);
    val2 = f(i-1,j,k);
    val3 = f(i-2,j,k);
    

    i=0;
    
	f(i-1,j,k) = val1;
    f(i-2,j,k) = val2;
    f(i-3,j,k) = val3;
    
    
    // 1 to 4 coupling
    i=0;
    
    val1 = f(i,j,k);
    val2 = f(i+1,j,k);
    val3 = f(i+2,j,k);
    
    i=p->knox-1;
    
	f(i+1,j,k) = val1;
    f(i+2,j,k) = val2;
    f(i+3,j,k) = val3;
    }
    
    
    if(cs==2)
    ILOOP
    KLOOP
    PCHECK
    {
    // 2 to 3 coupling
    j=p->knoy-1;
    
    val1 = f(i,j,k);
    val2 = f(i,j-1,k);
    val3 = f(i,j-2,k);
    
    j=0;
    
	f(i,j-1,k) = val1;
    f(i,j-2,k) = val2;
    f(i,j-3,k) = val3;
    
    // 3 to 2 coupling
    j=0;
    
    val1 = f(i,j,k);
    val2 = f(i,j+1,k);
    val3 = f(i,j+2,k);
    
    j=p->knoy-1;
    
	f(i,j+1,k) = val1;
    f(i,j+2,k) = val2;
    f(i,j+3,k) = val3;
    }
    
    
    if(cs==3)
    ILOOP
    JLOOP
    PCHECK
    {
    // 6 to 5 coupling
    k=p->knoz-1;
    
    val1 = f(i,j,k);
    val2 = f(i,j,k-1);
    val3 = f(i,j,k-2);
    
    k=0;
    
	f(i,j,k-1) = val1;
    f(i,j,k-2) = val2;
    f(i,j,k-3) = val3;
    
    // 5 to 6 coupling
    k=0;
    
    val1 = f(i,j,k);
    val2 = f(i,j,k+1);
    val3 = f(i,j,k+2);
    
    k=p->knoz-1;

	f(i,j,k+1) = val1;
    f(i,j,k+2) = val2;
    f(i,j,k+3) = val3;
    }
}

void ghostcell::gcV_periodic(lexer *p, vec &x, int gcv, int cs)
{
    double val1,val2,val3;
    
    if(cs==1)
    {
        for(q=0;q<p->gc4periodic_count[0];++q)
        {
        // 4 to 1 coupling
        n = p->gc4periodic[3][q];
        val1 = x.V[I_J_K_4];
        val2 = x.V[Im1_J_K_4];
        val3 = x.V[Im2_J_K_4]; 
        
        n = p->gc4periodic[0][q];
        x.V[Im1_J_K_4] = val1;
        x.V[Im2_J_K_4] = val2;
        x.V[Im3_J_K_4] = val3; 
        
        // 1 to 4 coupling
        n = p->gc4periodic[0][q];
        val1 = x.V[I_J_K_4];
        val2 = x.V[Ip1_J_K_4];
        val3 = x.V[Ip2_J_K_4]; 
        
        n = p->gc4periodic[3][q];
        x.V[Ip1_J_K_4] = val1;
        x.V[Ip2_J_K_4] = val2;
        x.V[Ip3_J_K_4] = val3;         
        }
    }
    
    if(cs==2)
    {
        for(q=0;q<p->gc4periodic_count[2];++q)
        {
        // 2 to 3 coupling
        n = p->gc4periodic[1][q];
        val1 = x.V[I_J_K_4];
        val2 = x.V[I_Jm1_K_4];
        val3 = x.V[I_Jm2_K_4]; 
        
        n = p->gc4periodic[2][q];
        x.V[I_Jm1_K_4] = val1;
        x.V[I_Jm2_K_4] = val2;
        x.V[I_Jm3_K_4] = val3; 
        
        // 3 to 2 coupling
        n = p->gc4periodic[2][q];
        val1 = x.V[I_J_K_4];
        val2 = x.V[I_Jp1_K_4];
        val3 = x.V[I_Jp2_K_4]; 
        
        n = p->gc4periodic[1][q];
        x.V[I_Jp1_K_4] = val1;
        x.V[I_Jp2_K_4] = val2;
        x.V[I_Jp3_K_4] = val3;         
        }
    }
    
    if(cs==3)
    {
        for(q=0;q<p->gc4periodic_count[4];++q)
        {
        // 6 to 5 coupling
        n = p->gc4periodic[5][q];
        val1 = x.V[I_J_K_4];
        val2 = x.V[I_J_Km1_4];
        val3 = x.V[I_J_Km2_4]; 
        
        n = p->gc4periodic[4][q];
        x.V[I_J_Km1_4] = val1;
        x.V[I_J_Km2_4] = val2;
        x.V[I_J_Km3_4] = val3; 
        
        // 5 to 6 coupling
        n = p->gc4periodic[4][q];
        val1 = x.V[I_J_K_4];
        val2 = x.V[I_J_Kp1_4];
        val3 = x.V[I_J_Kp2_4]; 
        
        n = p->gc4periodic[5][q];
        x.V[I_J_Kp1_4] = val1;
        x.V[I_J_Kp2_4] = val2;
        x.V[I_J_Kp3_4] = val3;         
        }
    }
}

void ghostcell::gcV_periodic_all(lexer *p, vec &x, int gcv, int cs)
{
    double val1,val2,val3;
    
    if(cs==1)
    {
        for(q=0;q<p->gc4periodic_count[0];++q)
        {
        // 4 to 1 coupling
        n = p->gc4aperiodic[3][q];
        val1 = x.V[I_J_K_4a];
        val2 = x.V[Im1_J_K_4a];
        val3 = x.V[Im2_J_K_4a]; 
        
        n = p->gc4aperiodic[0][q];
        x.V[Im1_J_K_4a] = val1;
        x.V[Im2_J_K_4a] = val2;
        x.V[Im3_J_K_4a] = val3; 
        
        // 1 to 4 coupling
        n = p->gc4aperiodic[0][q];
        val1 = x.V[I_J_K_4a];
        val2 = x.V[Ip1_J_K_4a];
        val3 = x.V[Ip2_J_K_4a]; 
        
        n = p->gc4aperiodic[3][q];
        x.V[Ip1_J_K_4a] = val1;
        x.V[Ip2_J_K_4a] = val2;
        x.V[Ip3_J_K_4a] = val3;         
        }
    }
}
