/*------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs
​
This file is part of REEF3D.
​
REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
​
This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.
​
You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Arun Kamath​
-----------*/
#include"force_fit.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include <math.h>

double force_fit::du_V_dt(lexer *p, fdm_fnpf *c, ghostcell *pgc) // calculate d(u*V)/dt four FK-force
{
    double du_V_dt = (u_V-u_V_n)/p->dt;
    return du_V_dt;
}

double force_fit::dv_V_dt(lexer *p, fdm_fnpf *c, ghostcell *pgc) // calculate d(v*V)/dt four FK-force
{
    double dv_V_dt = (v_V-v_V_n)/p->dt;
    return dv_V_dt;
}

double force_fit::dw_V_dt(lexer *p, fdm_fnpf *c, ghostcell *pgc) // calculate d(w*V)/dt four FK-force
{
    double dw_V_dt = (w_V-w_V_n)/p->dt;
    return dw_V_dt;
}

 void force_fit::force_fit_force(lexer* p, fdm_fnpf *c, ghostcell *pgc)
 {
    Fx= Fy= Fz=0.0;
    Fz_buoy=Fx_FK=Fy_FK=Fz_FK=0.0;
    u_V=v_V=w_V=0.0;
    
    for(k=0; k<p->knoz; ++k)
    {
        Fz_buoy+= (-1)*p->W1*p->W22*PI*rc*rc*p->DZN[KP]*(1/p->sigz[IJ]);
        u_V+=c->U[FIJK]*PI*rc*rc*p->DZN[KP]*(1/p->sigz[IJ]);
        v_V+=c->V[FIJK]*PI*rc*rc*p->DZN[KP]*(1/p->sigz[IJ]);
        w_V+=c->W[FIJK]*PI*rc*rc*p->DZN[KP]*(1/p->sigz[IJ]);
    }
    
    Fx_FK= (-1)*p->W1*du_V_dt(p,c,pgc);
    Fy_FK= (-1)*p->W1*dv_V_dt(p,c,pgc);
    Fz_FK= (-1)*p->W1*dw_V_dt(p,c,pgc);
    
    Fx=Fx_FK;
    Fy=Fy_FK;
    Fz=Fz_buoy+Fz_FK;
    
    u_V_n=u_V;
    v_V_n=v_V;
    w_V_n=v_V;
}