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

#include"sflow_diffusion_void.h"
#include"lexer.h"
#include"fdm2D.h"

sflow_diffusion_void::sflow_diffusion_void(lexer* p)
{
}

sflow_diffusion_void::~sflow_diffusion_void()
{
}

void sflow_diffusion_void::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{


}

void sflow_diffusion_void::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{


}

void sflow_diffusion_void::diff_w(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, slice &w, double alpha)
{


}

void sflow_diffusion_void::diff_scalar(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &f, double sig, double alpha)
{
    
}
