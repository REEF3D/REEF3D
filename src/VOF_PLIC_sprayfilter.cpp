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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"lexer.h"
#include"fdm.h"

void VOF_PLIC::sprayfilter(fdm* a, lexer* p)
{
    int voltest = 0;
    int jrange;
    if(p->j_dir>0)
        jrange=1;
    else 
        jrange=0;
        
   
        for(int irun=-1; irun<1; irun++)
        {
            for(int krun=-1; krun<1; krun++)
            {
                for(int jrun=-jrange; jrun<jrange; jrun++)
                {
                    if(vofstep(i+irun,j+jrun,k+krun)>=p->F91)
                        voltest=1;
                }
            }
        }
        
        if(voltest==0)
            vofstep(i,j,k)=0.0;
}