/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"fnpf_sg_v.h"

fnpf_sg_void::fnpf_sg_void()
{
}

fnpf_sg_void::~fnpf_sg_void()
{
}

void fnpf_sg_void::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	
	
}

void fnpf_sg_void::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc,ioflow *pflow, reini *preini, convection *pconvec)
{	
}

void fnpf_sg_void::inidisc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    
}
   