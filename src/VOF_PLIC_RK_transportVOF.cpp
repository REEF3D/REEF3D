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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

// in 2D scheme one function is used for Fields
void VOF_PLIC::vof_transport_COSMIC2D_RK
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    int nSweep,
    int sweep,
    field& uvel,
    field& vvel,
    field& wvel
)
{
    if(nSweep<2 && sweep==0)
    {   
        LOOP
        {
            Flux_x(i,j,k)=(Vn_p(i,j,k)-Vn_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

            F_x(i,j,k)=F_n(i,j,k)
                        -Flux_x(i,j,k)
                        +p->dt*F_n(i,j,k)*(uvel(i,j,k)-uvel(i-1,j,k))/p->DXN[IP];
        }
        pgc->start4(p,F_x,1);
        if(nSweep==1)
            LOOP
            {
                Flux_xz(i,j,k)=(Vz_p(i,j,k)-Vz_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
        
    }
    
    if(nSweep<2 && sweep==2)
    {
        LOOP
        {
            Flux_z(i,j,k)=(Vn_p(i,j,k)-Vn_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
        
            F_z(i,j,k)=F_n(i,j,k)
                        -Flux_z(i,j,k)
                        +p->dt*F_n(i,j,k)*(wvel(i,j,k)-wvel(i,j,k-1))/p->DZN[KP];
        }
        pgc->start4(p,F_z,1);
        if(nSweep==1)
            LOOP
            {
                Flux_zx(i,j,k)=(Vx_p(i,j,k)-Vx_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
    }
    
    if(nSweep==2)
    {
        LOOP
        {
            if(sweep==0)
            {
                Flux_zx(i,j,k)=(Vx_p(i,j,k)-Vx_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
            else if(sweep==2)
            {
                Flux_xz(i,j,k)=(Vz_p(i,j,k)-Vz_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
        
        }               
    }
    
}

//in 3D too many Fields are present -> individual functions for individual fields and fluxes

void VOF_PLIC::get_Fx_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel
)
{
    LOOP
    {
        Flux_x(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_x(i,j,k)=F_n(i,j,k)
                    -Flux_x(i,j,k)
                    +p->dt*F_n(i,j,k)*(uvel(i,j,k)-uvel(i-1,j,k))/p->DXN[IP];
    }
    pgc->start4(p,F_x,1);
}

void VOF_PLIC::get_Fy_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& vvel
)
{
    LOOP
    {
        Flux_y(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_y(i,j,k)=F_n(i,j,k)
                    -Flux_y(i,j,k)
                    +p->dt*F_n(i,j,k)*(vvel(i,j,k)-vvel(i,j-1,k))/p->DYN[JP];
    }
    pgc->start4(p,F_y,1);
}

void VOF_PLIC::get_Fz_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& wvel
)
{
    LOOP
    {
        Flux_z(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_z(i,j,k)=F_n(i,j,k)
                    -Flux_z(i,j,k)
                    +p->dt*F_n(i,j,k)*(wvel(i,j,k)-wvel(i,j,k-1))/p->DZN[KP];
    }
    pgc->start4(p,F_z,1);
}

void VOF_PLIC::get_Fxy_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel
)
{
    LOOP
    {
        Flux_xy(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_xy(i,j,k)=F_y(i,j,k)
                    -Flux_xy(i,j,k)
                    +p->dt*F_y(i,j,k)*(uvel(i,j,k)-uvel(i-1,j,k))/p->DXN[IP];
    }
    pgc->start4(p,F_xy,1);
}

void VOF_PLIC::get_Fxz_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel
)
{
    LOOP
    {
        Flux_xz(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_xz(i,j,k)=F_z(i,j,k)
                    -Flux_xz(i,j,k)
                    +p->dt*F_z(i,j,k)*(uvel(i,j,k)-uvel(i-1,j,k))/p->DXN[IP];
    }
    pgc->start4(p,F_xz,1);
}

void VOF_PLIC::get_Fyx_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& vvel
)
{
    LOOP
    {
        Flux_yx(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_yx(i,j,k)=F_x(i,j,k)
                    -Flux_yx(i,j,k)
                    +p->dt*F_x(i,j,k)*(vvel(i,j,k)-vvel(i,j-1,k))/p->DYN[JP];
    }
    pgc->start4(p,F_yx,1);
}

void VOF_PLIC::get_Fyz_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& vvel
)
{
    LOOP
    {
        Flux_yz(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_yz(i,j,k)=F_z(i,j,k)
                    -Flux_yz(i,j,k)
                    +p->dt*F_z(i,j,k)*(vvel(i,j,k)-vvel(i,j-1,k))/p->DYN[JP];
    }
    pgc->start4(p,F_yz,1);
}

void VOF_PLIC::get_Fzx_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& wvel
)
{
    LOOP
    {
        Flux_zx(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_zx(i,j,k)=F_x(i,j,k)
                    -Flux_zx(i,j,k)
                    +p->dt*F_x(i,j,k)*(wvel(i,j,k)-wvel(i,j,k-1))/p->DZN[KP];
    }
    pgc->start4(p,F_zx,1);
}

void VOF_PLIC::get_Fzy_and_Flux
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& wvel
)
{
    LOOP
    {
        Flux_zy(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

        F_zy(i,j,k)=F_y(i,j,k)
                    -Flux_zy(i,j,k)
                    +p->dt*F_y(i,j,k)*(wvel(i,j,k)-wvel(i,j,k-1))/p->DZN[KP];
    }
    pgc->start4(p,F_zy,1);
}

void VOF_PLIC::get_Flux_xzy
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel
)
{
    LOOP
        Flux_xzy(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
}

void VOF_PLIC::get_Flux_xyz
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel
)
{
    LOOP
        Flux_xyz(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
}

void VOF_PLIC::get_Flux_yxz
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& vvel
)
{
    LOOP
        Flux_yxz(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
}

void VOF_PLIC::get_Flux_yzx
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& vvel
)
{
    LOOP
        Flux_yzx(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
}

void VOF_PLIC::get_Flux_zxy
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& wvel
)
{
    LOOP
        Flux_zxy(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
}

void VOF_PLIC::get_Flux_zyx
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& wvel
)
{
    LOOP
        Flux_zyx(i,j,k)=(V_p(i,j,k)-V_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
}