/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4a.h"
#include"sediment_fdm.h"
#include"turbulence.h"

partres::partres(lexer *p) : drho(p->W1/p->S22) ,invKinVis(p->W1/p->W2), 
                                            Ps(p->Q14),beta(p->Q15),epsilon(p->Q16),theta_crit(p->Q17),bedChange(p)
{
        p->Darray(stressTensor,p->imax*p->jmax*p->kmax);
        p->Darray(cellSum,p->imax*p->jmax*p->kmax);
        p->Darray(cellSumTopo,p->imax*p->jmax*p->kmax);
        p->Darray(columnSum,p->imax*p->jmax);
        dx=p->global_xmax-p->global_xmin;
        LOOP
        {
            dx = min(dx,MIN3(p->DXN[IP],p->DYN[JP],p->DZN[KP]));
        }
        
        relax_ini(p);
}

partres::~partres()
{
        delete[] stressTensor;
        stressTensor = nullptr;
        delete[] cellSum;
        cellSum = nullptr;
        delete[] cellSumTopo;
        cellSumTopo = nullptr;
        delete[] columnSum;
        columnSum = nullptr;
}








