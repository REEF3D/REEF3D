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
#include"freesurface_header.h"


void VOF_PLIC::fieldloop_xy
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_y(i,j,k)>=a_thres && F_y(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_y);
                advectPlane_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
            }
            else if(F_y(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i-1,j,k)>0.0)
                V_m(i,j,k)=V_p(i-1,j,k);
                    
            if(V_m(i+1,j,k)<0.0)
                V_p(i,j,k)=V_m(i+1,j,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Fxy_and_Flux(a,p,pgc,uvel);
        swtch_xy=1;
        //F_xy and Flux_xy done
}
    
void VOF_PLIC::fieldloop_xz
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_z(i,j,k)>=a_thres && F_z(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_z);
                advectPlane_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
            }
            else if(F_z(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i-1,j,k)>0.0)
                V_m(i,j,k)=V_p(i-1,j,k);
                    
            if(V_m(i+1,j,k)<0.0)
                V_p(i,j,k)=V_m(i+1,j,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Fxz_and_Flux(a,p,pgc,uvel);
        swtch_xz=1;
        //F_xz and Flux_xz done
}

void VOF_PLIC::fieldloop_yx
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_x(i,j,k)>=a_thres && F_x(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_x);
                advectPlane_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
            }
            else if(F_x(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j-1,k)>0.0)
                V_m(i,j,k)=V_p(i,j-1,k);
                    
            if(V_m(i,j+1,k)<0.0)
                V_p(i,j,k)=V_m(i,j+1,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Fyx_and_Flux(a,p,pgc,vvel);
        swtch_yx=1;
        //F_yx and Flux_yx done
}

void VOF_PLIC::fieldloop_yz
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_z(i,j,k)>=a_thres && F_z(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_z);
                advectPlane_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
            }
            else if(F_z(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j-1,k)>0.0)
                V_m(i,j,k)=V_p(i,j-1,k);
                    
            if(V_m(i,j+1,k)<0.0)
                V_p(i,j,k)=V_m(i,j+1,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Fyz_and_Flux(a,p,pgc,vvel);
        swtch_yz=1;
        //F_yz and Flux_yz done
}

void VOF_PLIC::fieldloop_zx
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_x(i,j,k)>=a_thres && F_x(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_x);
                advectPlane_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
            }
            else if(F_x(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j,k-1)>0.0)
                V_m(i,j,k)=V_p(i,j,k-1);
                    
            if(V_m(i,j,k+1)<0.0)
                V_p(i,j,k)=V_m(i,j,k+1);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Fzx_and_Flux(a,p,pgc,wvel);
        swtch_zx=1;
        //F_zx and Flux_zx done
}

void VOF_PLIC::fieldloop_zy
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_y(i,j,k)>=a_thres && F_y(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_y);
                advectPlane_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
            }
            else if(F_y(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j,k-1)>0.0)
                V_m(i,j,k)=V_p(i,j,k-1);
                    
            if(V_m(i,j,k+1)<0.0)
                V_p(i,j,k)=V_m(i,j,k+1);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Fzy_and_Flux(a,p,pgc,wvel);
        swtch_zy=1;
        //F_zy and Flux_zy done
}
    
void VOF_PLIC::fieldloop_xzy
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_zy(i,j,k)>=a_thres && F_zy(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_zy);
                advectPlane_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
            }
            else if(F_zy(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i-1,j,k)>0.0)
                V_m(i,j,k)=V_p(i-1,j,k);
                    
            if(V_m(i+1,j,k)<0.0)
                V_p(i,j,k)=V_m(i+1,j,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Flux_xzy(a,p,pgc,uvel);
        swtch_xzy=1;
        //Flux_xzy done
}

void VOF_PLIC::fieldloop_xyz
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_yz(i,j,k)>=a_thres && F_yz(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_yz);
                advectPlane_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
            }
            else if(F_yz(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,0,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i-1,j,k)>0.0)
                V_m(i,j,k)=V_p(i-1,j,k);
                    
            if(V_m(i+1,j,k)<0.0)
                V_p(i,j,k)=V_m(i+1,j,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Flux_xyz(a,p,pgc,uvel);
        swtch_xyz=1;
        //Flux_xyz done
}

void VOF_PLIC::fieldloop_yxz
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_xz(i,j,k)>=a_thres && F_xz(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_xz);
                advectPlane_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
            }
            else if(F_xz(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j-1,k)>0.0)
                V_m(i,j,k)=V_p(i,j-1,k);
                    
            if(V_m(i,j+1,k)<0.0)
                V_p(i,j,k)=V_m(i,j+1,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Flux_yxz(a,p,pgc,vvel);
        swtch_yxz=1;
        //Flux_yxz done
}

void VOF_PLIC::fieldloop_yzx
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_zx(i,j,k)>=a_thres && F_zx(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_zx);
                advectPlane_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
            }
            else if(F_zx(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,1,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j-1,k)>0.0)
                V_m(i,j,k)=V_p(i,j-1,k);
                    
            if(V_m(i,j+1,k)<0.0)
                V_p(i,j,k)=V_m(i,j+1,k);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Flux_yzx(a,p,pgc,vvel);
        swtch_yzx=1;
        //Flux_yzx done
}

void VOF_PLIC::fieldloop_zxy
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_xy(i,j,k)>=a_thres && F_xy(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_xy);
                advectPlane_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
            }
            else if(F_xy(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j,k-1)>0.0)
                V_m(i,j,k)=V_p(i,j,k-1);
                    
            if(V_m(i,j,k+1)<0.0)
                V_p(i,j,k)=V_m(i,j,k+1);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Flux_zxy(a,p,pgc,wvel);
        swtch_zxy=1;
        //Flux_zxy done
}

void VOF_PLIC::fieldloop_zyx
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
        LOOP
        {
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_yx(i,j,k)>=a_thres && F_yx(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_yx);
                advectPlane_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
            }
            else if(F_yx(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,2,uvel,vvel,wvel);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        LOOP
        {
            if(V_p(i,j,k-1)>0.0)
                V_m(i,j,k)=V_p(i,j,k-1);
                    
            if(V_m(i,j,k+1)<0.0)
                V_p(i,j,k)=V_m(i,j,k+1);
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        get_Flux_zyx(a,p,pgc,wvel);
        swtch_zyx=1;
        //Flux_zyx done
}