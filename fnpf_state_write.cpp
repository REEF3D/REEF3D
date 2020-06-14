/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"fnpf_state.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void fnpf_state::write(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // mainheader file
    if(p->mpirank==0)
    mainheader(p,c,pgc);
    
    
    // Open File
	int num=0;

    if(p->P15>=1)
    num = printcount;

    
    // result file
    filename(p,c,pgc,num);
	 
	ofstream result;
	result.open(name, ios::binary);

    
    SLICELOOP4
    {
    ffn=float(c->eta(i,j));
    result.write((char*)&ffn, sizeof (float));
    } 
    
    FLOOP
    {
    ffn=float(c->U[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 

	FLOOP
    {
    ffn=float(c->V[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 

	FLOOP
    {
    ffn=float(c->W[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	
	result.close();
	
	++printcount;
}