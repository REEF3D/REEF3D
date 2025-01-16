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

void VOF_PLIC::stepwise_scheme
(
    fdm* a,
    lexer* p,
    ghostcell* pgc
)
{
    for(int nSweep=0 ; nSweep<Sweepdim; nSweep++)
    {
        LOOP
        {
            V_w_p(i,j,k)=0.0;
            V_w_m(i,j,k)=0.0;
        }
        pgc->start4(p,V_w_p,1);
        pgc->start4(p,V_w_m,1);
        
        if(p->j_dir>0)
            sweep=S_S[sSweep][nSweep];
        else
            sweep=S_2D[sSweep][nSweep];
        
        if(nSweep==0)
        {
            pgc->start4(p,a->phi,1);
            pgc->start4(p,a->vof,1);
            LOOP
            {
                phistep(i,j,k)=a->phi(i,j,k);
                vofstep(i,j,k)=a->vof(i,j,k);
            }
            pgc->start4(p,phistep,1);
            pgc->start4(p,vofstep,1);
        }
        else if(nSweep==1)
        {
            pgc->start4(p,phiS0,1);
            pgc->start4(p,vofS0,1);
            LOOP
            {
                phistep(i,j,k)=phiS0(i,j,k);
                vofstep(i,j,k)=vofS0(i,j,k);
            }
            pgc->start4(p,phistep,1);
            pgc->start4(p,vofstep,1);
        }
        else
        {
            pgc->start4(p,phiS1,1);
            pgc->start4(p,vofS1,1);
            LOOP
            {
                phistep(i,j,k)=phiS1(i,j,k);
                vofstep(i,j,k)=vofS1(i,j,k);
            }
            pgc->start4(p,phistep,1);
            pgc->start4(p,vofstep,1);
        }
        
        LOOP
        {
            switch(p->F89)
                {
                    case 0:
                        break;
                    case 1:
                        transportPhi_Bonn(a,p,nSweep,sweep);
                        break;
                    case 2:
                        break;
                }
            
            
            if((vofstep(i,j,k)>0.0001 && vofstep(i,j,k)<0.9999))
            {
                reconstructPlane_alt(a,p,vofstep);
                
                switch(p->F90)
                {
                    case 0:
                        break;
                    case 1:
                        advectPlane_forBonnScheme(a,p,sweep);
                        break;
                    case 3:
                        advectPlane_NewWang(a,p,sweep);
                        break;
                }
                
            }
            else if( vofstep(i,j,k)>0.9999)
            {   
                switch(p->F90)
                {
                    case 0:
                        break;
                    case 1:
                        advectWater_forBonnScheme(a,p,sweep);
                        break;
                    case 3:
                        advectPlane_NewWang(a,p,sweep);
                        break;
                }
                
            }
            else
            {
                switch(p->F90)
                {
                    case 0:
                        break;
                    case 1:
                        break;
                    case 3:
                        advectPlane_NewWang(a,p,sweep);
                        break;
                }
            }
        }
        
        pgc->start4(p,nx,1);
        pgc->start4(p,ny,1);
        pgc->start4(p,nz,1);
        pgc->start4(p,alpha,1);
        pgc->start4(p,V_w_p,1);
        pgc->start4(p,V_w_m,1);
        pgc->start4(p,phiS1,1);
        pgc->start4(p,vofS1,1);
        pgc->start4(p,phiS0,1);
        pgc->start4(p,vofS0,1);
        
        switch(p->F90)
        {
            case 0:
                break;
            case 1:
                transportVOF_Bonn(a,p,nSweep,sweep);
                break;
            case 3:
                transportVOF_NewWang(a,p,nSweep);
                break;
        }
        
        pgc->start4(p,phiS1,1);
        pgc->start4(p,vofS1,1);
        pgc->start4(p,phiS0,1);
        pgc->start4(p,vofS0,1);
        pgc->start4(p,vofS2,1);
        
        
        if(p->j_dir>0)
        {
             if(nSweep==2)
             {
                 pgc->start4(p,phiS2,1);
                 pgc->start4(p,vofS2,1);
                 LOOP
                 {
                     phistep(i,j,k)=phiS2(i,j,k);
                     vofstep(i,j,k)=vofS2(i,j,k);
                 }
                 pgc->start4(p,phistep,1);
                 pgc->start4(p,vofstep,1);
             }
        }
        else
        {
            if(nSweep==1)
            {
                pgc->start4(p,phiS1,1);
                pgc->start4(p,vofS1,1);
                LOOP
                {
                    phistep(i,j,k)=phiS1(i,j,k);
                    vofstep(i,j,k)=vofS1(i,j,k);
                }
                pgc->start4(p,phistep,1);
                pgc->start4(p,vofstep,1);
            }
        }
    }
    
    pgc->start4(p,phistep,1);
    pgc->start4(p,vofstep,1);
}

void VOF_PLIC::symmetric_scheme2D
(
    fdm* a,
    lexer* p,
    ghostcell* pgc
)
{
    LOOP
    {
        F_n(i,j,k)=a->vof(i,j,k);
        F_x(i,j,k)=0.0;
        F_z(i,j,k)=0.0;
        F_new(i,j,k)=0.0;
        Crossflux_xz(i,j,k)=0.0;
        Crossflux_zx(i,j,k)=0.0;
    }
    pgc->start4(p,F_n,1);
    pgc->start4(p,F_x,1);
    pgc->start4(p,F_z,1);
    pgc->start4(p,F_new,1);
    pgc->start4(p,Crossflux_xz,1);
    pgc->start4(p,Crossflux_zx,1);
    
    for(int nSweep=0; nSweep<Sweepdim; nSweep++)
    {
        if(p->j_dir>0)
            sweep=S_S[sSweep][nSweep];
        else
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
            if(F_n(i,j,k)>=0.001 && F_n(i,j,k)<=0.999)
            {
                reconstructPlane_alt(a,p,F_n);
                advectPlane_forCOSMIC2D_simple(a,p,sweep,-1);
            }
            else if(F_n(i,j,k)>0.999)
                advectWater_forCOSMIC2D_simple(a,p,sweep,-1);
                
            if(nSweep==1 && sweep==0)
            {
                if(F_z(i,j,k)>=0.001 && F_z(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_z);
                    advectPlane_forCOSMIC2D_simple(a,p,sweep,2);
                }
                else if(F_z(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,sweep,2);
            }
            else if(nSweep==1 && sweep==2)
            {
                if(F_x(i,j,k)>=0.001 && F_x(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_x);
                    advectPlane_forCOSMIC2D_simple(a,p,sweep,0);
                }
                else if(F_x(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,sweep,0);
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
                    
        vof_transport_COSMIC2D(a,p,nSweep,sweep);
    }
    
    if(sweep==0)
    {
            LOOP
            {
                if(F_x(i,j,k)>=0.001 && F_x(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_x);
                    advectPlane_forCOSMIC2D_simple(a,p,2,0);
                }
                else if(F_x(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,2,0);
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
                if(F_z(i,j,k)>=0.001 && F_z(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_z);
                    advectPlane_forCOSMIC2D_simple(a,p,0,2);
                }
                else if(F_z(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,0,2);
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
    
    vof_transport_COSMIC2D(a,p,2,sweep);
   // cout<<"F_n:"<<F_n(5,0,5)<<" ;F_x:"<<F_x(5,0,5)<<" ;F_z:"<<F_z(5,0,5)<<endl;
   // cout<<"Flux_x:"<<Flux_x(5,0,5)<<" ;Flux_z:"<<Flux_z(5,0,5)<<" ;Crossflux_xz:"<<Crossflux_xz(5,0,5)<<" ;Crossflux_zx:"<<Crossflux_zx(5,0,5)<<endl;
   // cout<<"F_new:"<<F_new(5,0,5)<<endl;
    
    LOOP
    {
        vofstep(i,j,k)=F_new(i,j,k);
    }
   // cout<<"vofstep:"<<vofstep(5,0,5)<<endl;
}

void VOF_PLIC::symmetric_scheme2D_FCRK3
(
    fdm* a,
    lexer* p,
    ghostcell* pgc
)
{
    LOOP
    {
        F_n(i,j,k)=a->vof(i,j,k);
        F_x(i,j,k)=0.0;
        F_z(i,j,k)=0.0;
        F_new(i,j,k)=0.0;
        Crossflux_xz(i,j,k)=0.0;
        Crossflux_zx(i,j,k)=0.0;
    }
    pgc->start4(p,F_n,1);
    pgc->start4(p,F_x,1);
    pgc->start4(p,F_z,1);
    pgc->start4(p,F_new,1);
    pgc->start4(p,Crossflux_xz,1);
    pgc->start4(p,Crossflux_zx,1);
    
    for(int nSweep=0; nSweep<Sweepdim; nSweep++)
    {
        if(p->j_dir>0)
            sweep=S_S[sSweep][nSweep];
        else
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
            if(F_n(i,j,k)>=0.001 && F_n(i,j,k)<=0.999)
            {
                reconstructPlane_alt(a,p,F_n);
                advectPlane_forCOSMIC2D_simple(a,p,sweep,-1);
            }
            else if(F_n(i,j,k)>0.999)
                advectWater_forCOSMIC2D_simple(a,p,sweep,-1);
                
            if(nSweep==1 && sweep==0)
            {
                if(F_z(i,j,k)>=0.001 && F_z(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_z);
                    advectPlane_forCOSMIC2D_simple(a,p,sweep,2);
                }
                else if(F_z(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,sweep,2);
            }
            else if(nSweep==1 && sweep==2)
            {
                if(F_x(i,j,k)>=0.001 && F_x(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_x);
                    advectPlane_forCOSMIC2D_simple(a,p,sweep,0);
                }
                else if(F_x(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,sweep,0);
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
                    
        vof_transport_COSMIC2D(a,p,nSweep,sweep);
    }
    
    if(sweep==0)
    {
            LOOP
            {
                if(F_x(i,j,k)>=0.001 && F_x(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_x);
                    advectPlane_forCOSMIC2D_simple(a,p,2,0);
                }
                else if(F_x(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,2,0);
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
                if(F_z(i,j,k)>=0.001 && F_z(i,j,k)<=0.999)
                {
                    reconstructPlane_alt(a,p,F_z);
                    advectPlane_forCOSMIC2D_simple(a,p,0,2);
                }
                else if(F_z(i,j,k)>0.999)
                    advectWater_forCOSMIC2D_simple(a,p,0,2);
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
    
    vof_transport_COSMIC2D(a,p,2,sweep);
   // cout<<"F_n:"<<F_n(5,0,5)<<" ;F_x:"<<F_x(5,0,5)<<" ;F_z:"<<F_z(5,0,5)<<endl;
   // cout<<"Flux_x:"<<Flux_x(5,0,5)<<" ;Flux_z:"<<Flux_z(5,0,5)<<" ;Crossflux_xz:"<<Crossflux_xz(5,0,5)<<" ;Crossflux_zx:"<<Crossflux_zx(5,0,5)<<endl;
   // cout<<"F_new:"<<F_new(5,0,5)<<endl;
    
    LOOP
    {
        a->L(i,j,k)=-0.5*(Flux_x(i,j,k)+Crossflux_zx(i,j,k)+Flux_z(i,j,k)+Crossflux_xz(i,j,k));
    }
    pgc->start4(p,a->L,1);
   // cout<<"vofstep:"<<vofstep(5,0,5)<<endl;
}