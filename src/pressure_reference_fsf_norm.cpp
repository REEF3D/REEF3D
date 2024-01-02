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

#include"pressure_reference.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void pressure_reference::fsf_normalize(lexer*p, fdm* a, ghostcell *pgc)
{
    double epsi;
	double dirac;
    double pressval;
    double dirac_sum;

    
    // epsi
    if(p->j_dir==0)        
    epsi = 2.1*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    epsi = 2.1*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);


    // pressval dirac

    pressval=0.0;
    dirac_sum=0.0;
	LOOP
	{
        if(fabs(a->phi(i,j,k))<epsi)
        dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
            
        if(fabs(a->phi(i,j,k))>=epsi)
        dirac=0.0;
        
        if(dirac>1.0e-10 && a->phi(i,j,k)<0.0)
        {
        pressval += dirac*a->press(i,j,k);
        dirac_sum += dirac;
        }
    }
    
    pressval = pgc->globalsum(pressval);
    
    dirac_sum = pgc->globalsum(dirac_sum);
    
    if(dirac_sum>0)
    pressval = pressval/dirac_sum;
    
    
    if(p->B33==1)
    p->pressgage=pressval;
    
    if(p->B33==2)
    LOOP
    a->press(i,j,k) -= pressval;
    
    
    // pressval orig
    /*
    pressval=0.0;
    count=0;
	LOOP
	{
        if(fabs(a->phi(i,j,k))<epsi)
        dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
            
        if(fabs(a->phi(i,j,k))>=epsi)
        dirac=0.0;
        
        if(dirac>1.0e-10 && a->phi(i,j,k)<0.0)
        {
        pressval += a->press(i,j,k);
        ++count;
        }
	}
    
    pressval = pgc->globalsum(pressval);
    
    count = pgc->globalisum(count);
    
    if(count>0)
    pressval = pressval/double(count);
    
    LOOP
    a->press(i,j,k) -= pressval;
    */

}




