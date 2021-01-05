/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"patchBC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"


void patchBC::patchBC_IDcount(lexer *p, ghostcell *pgc)
{    
    int check;
    
    geo_count = p->B221 + p->B222;
    
    p->Iarray(ID_array,geo_count);
    
    
    // ini ID array
    count=0;
    check=0;
    
    if(p->B221>0 && check==0)
    {
    ID_array[0] = p->B221_ID[0];
    count=1;
    check=1;
    }
    
    if(p->B222>0 && check==0)
    {
    ID_array[0] = p->B222_ID[0];
    count=1;
    check=1;
    }
    
    // fill ID array
    for(n=0; n<p->B221;++n)
    {
        check=1;
        for(qn=0;qn<count;++qn)
        if(ID_array[qn] == p->B221_ID[n]);
        check=0;
        
        if(check==1)
        {
        ID_array[count] = p->B221_ID[n];
        ++count;       
        }
    }
    
    for(n=0; n<p->B222;++n)
    {
        check=1;
        for(qn=0;qn<count;++qn)
        if(ID_array[qn] == p->B221_ID[n]);
        check=0;
        
        if(check==1)
        {
        ID_array[count] = p->B221_ID[n];
        ++count;       
        }
    }
    
    obj_count=count;

} 

