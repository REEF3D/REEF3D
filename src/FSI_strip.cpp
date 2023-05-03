/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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
--------------------------------------------------------------------*/

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fsi_strip::fsi_strip(int num):nstrip(num),beam(num){}
    
fsi_strip::~fsi_strip(){}

void fsi_strip::start(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
	// Set mooring time step
	t_strip_n = t_strip;
	t_strip += alpha*p->dt;

    // Integrate from t_mooring_n to t_mooring
    Integrate(t_strip_n,t_strip);
};

void fsi_strip::update_points()
{
    getTransPos(x_el);
    getTransVel(xdot_el);
    getRotPos(q_el);
    getRotVel(qdot_el);

    for (int eI = 0; eI < Ne; eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            Xil[eI].col(pI) = rotVec(Xil_0[eI].col(pI),q_el.col(eI+1));
            lagrangePoints[eI].col(pI) = (x_el.col(eI+1) + x_el.col(eI))/2.0 + Xil[eI].col(pI);
        }
    }
}

void fsi_strip::store_variables(lexer *p)
{
    P_el_n = P_el;
    I_el_n = I_el;
}
