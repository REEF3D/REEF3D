/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"lexer.h"


lexer::lexer() : cmu(0.09), position(this), interpolation(this), grid_sigma(this)
{
    sigT=0.9;
    
	ini_default();
    mpirank=0;
}

void lexer::lexer_read()
{

    if(mpirank==0)
	read_control();

	ctrlsize=12500;
	
    Iarray(ictrl,ctrlsize);
    Darray(dctrl,ctrlsize);

	read_grid();
	
	lexer_ini();
    
    if(G2==1)
    sigma_coord_ini(this);
    
    lexer_gridspacing();
	parse();	
}

lexer::~lexer()
{
}

int lexer::xmax,lexer::ymax,lexer::zmax;
int lexer::knox,lexer::knoy,lexer::knoz;
int lexer::margin;








