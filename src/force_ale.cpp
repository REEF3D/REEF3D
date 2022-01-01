/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"force_ale.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

force_ale::force_ale(lexer* p, fdm *a, ghostcell *pgc, int qn) : ID(qn)
{
	
    force_aleprintcount=0;
        
    // open files
    print_ini(p,a,pgc);
    
    is = p->posc_i(p->P85_x[ID]);
    js = p->posc_j(p->P85_y[ID]);

}

force_ale::~force_ale()
{
}
void force_ale::ini(lexer *p, fdm *a, ghostcell *pgc)
{

} 

void force_ale::start(lexer *p, fdm *a, ghostcell *pgc)
{
} 

