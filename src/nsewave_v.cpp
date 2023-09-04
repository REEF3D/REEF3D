/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nsewave_v.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heat.h"
#include"concentration.h"

nsewave_v::nsewave_v(lexer *p, fdm *a, ghostcell *pgc, heat *&pheat, concentration *&pconc)
{
}

nsewave_v::~nsewave_v()
{
}

void nsewave_v::start(lexer* p, fdm* a, ghostcell* pgc, momentum *pmom, diffusion *pdiff, turbulence *pturb,
                      convection* pconvec, pressure *ppress, poisson *ppois, solver *ppoissonsolv, solver *psolv, 
                      ioflow* pflow, vrans* pvrans, sixdof_df_base *p6dof_df, vector<net*>& pnet)
{
}

void nsewave_v::update(lexer *p, fdm *a, ghostcell *pgc, slice &f)
{
}

void nsewave_v::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{

}
