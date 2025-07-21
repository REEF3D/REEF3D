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
#include"initialize.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"
#include"interpolation.h"

void VOF_PLIC::symmetric_scheme2D_FCRK3
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
        F_n(i,j,k)=a->vof(i,j,k);
    }
    pgc->start4(p,F_n,1);
    
    for(int nSweep=0; nSweep<Sweepdim; nSweep++)
    {
        sweep=S_2D[sSweep][nSweep];
        
        LOOP
        {   
            Vn_p(i,j,k)=0.0;
            Vn_m(i,j,k)=0.0;
            Vx_p(i,j,k)=0.0;
            Vx_m(i,j,k)=0.0;
            Vz_p(i,j,k)=0.0;
            Vz_m(i,j,k)=0.0;
        }
        pgc->start4(p,Vn_p,1);
        pgc->start4(p,Vn_m,1);
        pgc->start4(p,Vx_p,1);
        pgc->start4(p,Vx_m,1);
        pgc->start4(p,Vz_p,1);
        pgc->start4(p,Vz_m,1);
        
        LOOP
        {
            if(F_n(i,j,k)>=a_thres && F_n(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_n);
                advectPlane_forCOSMIC2D_RK(a,p,sweep,-1,uvel,vvel,wvel);
            }
            else if(F_n(i,j,k)>w_thres)
                advectWater_forCOSMIC2D_RK(a,p,sweep,-1,uvel,vvel,wvel);
                
            if(nSweep==1 && sweep==0)
            {
                if(F_z(i,j,k)>=a_thres && F_z(i,j,k)<=w_thres)
                {
                    reconstructPlane_alt(a,p,F_z);
                    advectPlane_forCOSMIC2D_RK(a,p,sweep,2,uvel,vvel,wvel);
                }
                else if(F_z(i,j,k)>w_thres)
                    advectWater_forCOSMIC2D_RK(a,p,sweep,2,uvel,vvel,wvel);
            }
            else if(nSweep==1 && sweep==2)
            {
                if(F_x(i,j,k)>=a_thres && F_x(i,j,k)<=w_thres)
                {
                    reconstructPlane_alt(a,p,F_x);
                    advectPlane_forCOSMIC2D_RK(a,p,sweep,0,uvel,vvel,wvel);
                }
                else if(F_x(i,j,k)>w_thres)
                    advectWater_forCOSMIC2D_RK(a,p,sweep,0,uvel,vvel,wvel);
            }
        }
        
        pgc->start4(p,Vn_p,1);
        pgc->start4(p,Vn_m,1);
        pgc->start4(p,Vx_p,1);
        pgc->start4(p,Vx_m,1);
        pgc->start4(p,Vz_p,1);
        pgc->start4(p,Vz_m,1);
        
        if(sweep==0)
        {
            LOOP
            {
                if(Vn_p(i-1,j,k)>0.0)
                    Vn_m(i,j,k)=Vn_p(i-1,j,k);
                    
                if(Vn_m(i+1,j,k)<0.0)
                    Vn_p(i,j,k)=Vn_m(i+1,j,k);
                    
                if(Vz_p(i-1,j,k)>0.0)
                    Vz_m(i,j,k)=Vz_p(i-1,j,k);
                    
                if(Vz_m(i+1,j,k)<0.0)
                    Vz_p(i,j,k)=Vz_m(i+1,j,k);
            }
            pgc->start4(p,Vn_p,1);
            pgc->start4(p,Vn_m,1);
            pgc->start4(p,Vz_p,1);
            pgc->start4(p,Vz_m,1);
        }
        else if(sweep==2)
        {
            LOOP
            {
                if(Vn_p(i,j,k-1)>0.0)
                    Vn_m(i,j,k)=Vn_p(i,j,k-1);
                    
                if(Vn_m(i,j,k+1)<0.0)
                    Vn_p(i,j,k)=Vn_m(i,j,k+1);
            
                if(Vx_p(i,j,k-1)>0.0)
                    Vx_m(i,j,k)=Vx_p(i,j,k-1);
                    
                if(Vx_m(i,j,k+1)<0.0)
                    Vx_p(i,j,k)=Vx_m(i,j,k+1);
            }        
            pgc->start4(p,Vn_p,1);
            pgc->start4(p,Vn_m,1);
            pgc->start4(p,Vx_p,1);
            pgc->start4(p,Vx_m,1);
        }
                    
        vof_transport_COSMIC2D_RK(a,p,pgc,nSweep,sweep,uvel,vvel,wvel);
    }
    
    if(sweep==0)
    {
            LOOP
            {
                if(F_x(i,j,k)>=a_thres && F_x(i,j,k)<=w_thres)
                {
                    reconstructPlane_alt(a,p,F_x);
                    advectPlane_forCOSMIC2D_RK(a,p,2,0,uvel,vvel,wvel);
                }
                else if(F_x(i,j,k)>w_thres)
                    advectWater_forCOSMIC2D_RK(a,p,2,0,uvel,vvel,wvel);
            }
            
            pgc->start4(p,Vx_p,1);
            pgc->start4(p,Vx_m,1);
            LOOP
            {
                if(Vx_p(i,j,k-1)>0.0)
                    Vx_m(i,j,k)=Vx_p(i,j,k-1);
                    
                if(Vx_m(i,j,k+1)<0.0)
                    Vx_p(i,j,k)=Vx_m(i,j,k+1);
            }
            pgc->start4(p,Vx_p,1);
            pgc->start4(p,Vx_m,1);
         
    }
    else if(sweep==2)
    {
            LOOP
            {
                if(F_z(i,j,k)>=a_thres && F_z(i,j,k)<=w_thres)
                {
                    reconstructPlane_alt(a,p,F_z);
                    advectPlane_forCOSMIC2D_RK(a,p,0,2,uvel,vvel,wvel);
                }
                else if(F_z(i,j,k)>w_thres)
                    advectWater_forCOSMIC2D_RK(a,p,0,2,uvel,vvel,wvel);
            }
            pgc->start4(p,Vz_p,1);
            pgc->start4(p,Vz_m,1);
            LOOP
            {
                if(Vz_p(i-1,j,k)>0.0)
                    Vz_m(i,j,k)=Vz_p(i-1,j,k);
                    
                if(Vz_m(i+1,j,k)<0.0)
                    Vz_p(i,j,k)=Vz_m(i+1,j,k);
            }
            pgc->start4(p,Vz_p,1);
            pgc->start4(p,Vz_m,1);
            
           
    }
    
    vof_transport_COSMIC2D_RK(a,p,pgc,2,sweep,uvel,vvel,wvel);
   // cout<<"F_n:"<<F_n(5,0,5)<<" ;F_x:"<<F_x(5,0,5)<<" ;F_z:"<<F_z(5,0,5)<<endl;
   // cout<<"Flux_x:"<<Flux_x(5,0,5)<<" ;Flux_z:"<<Flux_z(5,0,5)<<" ;Flux_xz:"<<Flux_xz(5,0,5)<<" ;Flux_zx:"<<Flux_zx(5,0,5)<<endl;
   // cout<<"F_new:"<<F_new(5,0,5)<<endl;
    
    LOOP
    {
        a->L(i,j,k)=-0.5*(Flux_x(i,j,k)+Flux_zx(i,j,k)+Flux_z(i,j,k)+Flux_xz(i,j,k));
    }
    pgc->start4(p,a->L,1);
   // cout<<"vofstep:"<<vofstep(5,0,5)<<endl;
}

void VOF_PLIC::symmetric_scheme3D_FCRK3
(
    fdm* a,
    lexer* p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
    swtch_x=0;
    swtch_y=0;
    swtch_z=0;
    swtch_xy=0;
    swtch_xz=0;
    swtch_yx=0;
    swtch_yz=0;
    swtch_zx=0;
    swtch_zy=0;
    swtch_xyz=0;
    swtch_xzy=0;
    swtch_yxz=0;
    swtch_yzx=0;
    swtch_zxy=0;
    swtch_zyx=0;
    
    LOOP
        F_n(i,j,k)=a->vof(i,j,k);
    pgc->start4(p,F_n,1);
        
    
    for(int nSweep=0; nSweep<Sweepdim; nSweep++)
    {
        
        sweep=S_S[sSweep][nSweep];
        LOOP
        {   
            V_p(i,j,k)=0.0;
            V_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        LOOP
        {
            if(F_n(i,j,k)>=a_thres && F_n(i,j,k)<=w_thres)
            {
                reconstructPlane_alt(a,p,F_n);
                advectPlane_forCOSMIC3D_RK(a,p,sweep,uvel,vvel,wvel);
            }
            else if(F_n(i,j,k)>w_thres)
                advectWater_forCOSMIC3D_RK(a,p,sweep,uvel,vvel,wvel);
        }
        
        pgc->start4(p,V_p,1);
        pgc->start4(p,V_m,1);
        
        if(sweep==0)
        {
            LOOP
            {
                if(V_p(i-1,j,k)>0.0)
                    V_m(i,j,k)=V_p(i-1,j,k);
                    
                if(V_m(i+1,j,k)<0.0)
                    V_p(i,j,k)=V_m(i+1,j,k);
            }
            pgc->start4(p,V_p,1);
            pgc->start4(p,V_m,1);
            get_Fx_and_Flux(a,p,pgc,uvel);
            swtch_x=1;
            //F_x and Flux_x done
        }
        else if(sweep==1)
        {
            LOOP
            {
                if(V_p(i,j-1,k)>0.0)
                    V_m(i,j,k)=V_p(i,j-1,k);
                    
                if(V_m(i,j+1,k)<0.0)
                    V_p(i,j,k)=V_m(i,j+1,k);
            }
            pgc->start4(p,V_p,1);
            pgc->start4(p,V_m,1);
            get_Fy_and_Flux(a,p,pgc,vvel);
            swtch_y=1;
            //F_y and Flux_y done
        }
        else if(sweep==2)
        {
            LOOP
            {
                if(V_p(i,j,k-1)>0.0)
                    V_m(i,j,k)=V_p(i,j,k-1);
                    
                if(V_m(i,j,k+1)<0.0)
                    V_p(i,j,k)=V_m(i,j,k+1);
            }        
            pgc->start4(p,V_p,1);
            pgc->start4(p,V_m,1);
            get_Fz_and_Flux(a,p,pgc,wvel);
            swtch_z=1;
            //F_z and Flux_z done
        }
        
        if(sweep==0)
        {
            if(swtch_y==1)
                fieldloop_xy(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_z==1)
                fieldloop_xz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_yz==1)
                fieldloop_xyz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_zy==1)
                fieldloop_xzy(a,p,pgc,uvel,vvel,wvel);
        }
        else if(sweep==1)
        {
            if(swtch_x==1)
                fieldloop_yx(a,p,pgc,uvel,vvel,wvel);
            
            if(swtch_z==1)
                fieldloop_yz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_xz==1)
                fieldloop_yxz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_zx==1)
                fieldloop_yzx(a,p,pgc,uvel,vvel,wvel);
        }
        else if(sweep==2)
        {
            if(swtch_x==1)
                fieldloop_zx(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_y==1)
                fieldloop_zy(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_xy==1)
                fieldloop_zxy(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_yx==1)
                fieldloop_zyx(a,p,pgc,uvel,vvel,wvel);
        }
        
    }
    for(int nSweep=0; nSweep<Sweepdim; nSweep++)
    {
        sweep=S_S[sSweep][nSweep];
        
        if(sweep==0)
        {
            if(swtch_xy==0)
                fieldloop_xy(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_xz==0)
                fieldloop_xz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_yz==1)
                fieldloop_xyz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_zy==1)
                fieldloop_xzy(a,p,pgc,uvel,vvel,wvel);
        }
        else if (sweep==1)
        {
            if(swtch_yx==0)
                fieldloop_yx(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_yz==0)
                fieldloop_yz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_xz==1)
                fieldloop_yxz(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_zx==1)
                fieldloop_yzx(a,p,pgc,uvel,vvel,wvel);
        }
        else if (sweep==2)
        {
            if(swtch_zx==0)
                fieldloop_zx(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_zy==0)
                fieldloop_zy(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_xy==1)
                fieldloop_zxy(a,p,pgc,uvel,vvel,wvel);
                
            if(swtch_yx==1)
                fieldloop_zyx(a,p,pgc,uvel,vvel,wvel);
        }
    }
    
    sweep=S_S[sSweep][0];
    
    if(sweep==0)
    {
        if(swtch_xyz==0)
            fieldloop_xyz(a,p,pgc,uvel,vvel,wvel);
            
        if(swtch_xzy==0)
            fieldloop_xzy(a,p,pgc,uvel,vvel,wvel);
    }
    else if(sweep==1)
    {
        if(swtch_yxz==0)
            fieldloop_yxz(a,p,pgc,uvel,vvel,wvel);
            
        if(swtch_yzx==0)
            fieldloop_yzx(a,p,pgc,uvel,vvel,wvel);
    }
    else if(sweep==2)
    {
        if(swtch_zxy==0)
            fieldloop_zxy(a,p,pgc,uvel,vvel,wvel);
            
        if(swtch_zyx==0)
            fieldloop_zyx(a,p,pgc,uvel,vvel,wvel);
    }
    
    int swtchsum;
    swtchsum = swtch_x+swtch_y+swtch_z+swtch_xy+swtch_xz+swtch_yx+swtch_yz+swtch_zx+swtch_zy+swtch_zx+swtch_zy
                +swtch_xyz+swtch_xzy+swtch_yxz+swtch_yzx+swtch_zxy+swtch_zyx;
                
    if(swtchsum < 15)
        cout<<"not all switches activated!!!"<<endl;
                    
        
   
    LOOP
    {
        a->L(i,j,k)=-1.0/6.0*(2.0*Flux_x(i,j,k)+Flux_xy(i,j,k)+Flux_xz(i,j,k)+Flux_xyz(i,j,k)+Flux_xzy(i,j,k)
                             +2.0*Flux_y(i,j,k)+Flux_yx(i,j,k)+Flux_yz(i,j,k)+Flux_yxz(i,j,k)+Flux_yzx(i,j,k)
                             +2.0*Flux_z(i,j,k)+Flux_zx(i,j,k)+Flux_zy(i,j,k)+Flux_zxy(i,j,k)+Flux_zyx(i,j,k));
    }
    pgc->start4(p,a->L,1);
   // cout<<"vofstep:"<<vofstep(5,0,5)<<endl;
}

