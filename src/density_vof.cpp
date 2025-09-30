/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"initialize.h"
#include"density_vof.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

density_vof::density_vof(lexer* p) : epsi(p->F45*p->DXM), eps(2.1*p->DXM)
{
    
    
    double psim;
    int count;
    
    if(p->j_dir==0)        
    p->psi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    p->psi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
    
    
    p->psi0=p->psi;
}

density_vof::~density_vof()
{
}

double density_vof::roface(lexer *p, fdm *a, int d_i, int d_j, int d_k)
{  
    if(p->F92==1)   
    {
        double phival, psiro;
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+d_i,j+d_j,k+d_k));
        psiro = p->psi;
        
        if(phival>psiro)
            H=1.0;

        if(phival<-psiro)
            H=0.0;

        if(fabs(phival)<=psiro)
            H=0.5*(1.0 + phival/(psiro) + (1.0/PI)*sin((PI*phival)/(psiro)));
    
        roval = p->W1*H + p->W3*(1.0-H);
    }
    
    else if(p->F92==2||p->F92==22)
    {
            if(a->vof(i,j,k)>p->F94 && a->vof(i+d_i,j+d_j,k+d_k)>p->F94)
                roval=p->W1;
            else if(a->vof(i,j,k)<p->F93 && a->vof(i+d_i,j+d_j,k+d_k)<p->F93)
                roval=p->W3;
            else
            {
                if(d_i==1)
                    roval=(0.5*p->DXN[IP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                    +0.5*p->DXN[IP1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                    )/p->DXP[IP];
                else if(d_i==-1)
                    roval=(0.5*p->DXN[IP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                    +0.5*p->DXN[IM1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                    )/p->DXP[IM1];
            else if(d_k==1)
                roval=(0.5*p->DZN[KP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                +0.5*p->DZN[KP1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                )/p->DZP[KP];
            else if(d_k==-1)
                roval=(0.5*p->DZN[KP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                +0.5*p->DZN[KM1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                )/p->DZP[KM1];
            else if(d_j==1)
                roval=(0.5*p->DYN[JP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                +0.5*p->DYN[JP1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                )/p->DYP[JP];
            else if(d_j==-1)
                roval=(0.5*p->DYN[JP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                +0.5*p->DYN[JM1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                )/p->DYP[JM1];
            else
                cout<<"density case missing, d_i:"<<d_i<<" d_j:"<<d_j<<" d_k:"<<d_k<<endl;
            }
    }
    else if(p->F92==3||p->F92==32)
    {
        double H;
        int excepcheck;
        excepcheck=0;
        if(p->F96==1)
        {
            if(p->pos_x()>=p->F96_xs && p->pos_x()<=p->F96_xe)
            {
                excepcheck=1;
                
                if(a->vof(i,j,k)>p->F94 && a->vof(i+d_i,j+d_j,k+d_k)>p->F94)
                    roval=p->W1;
                else if(a->vof(i,j,k)<p->F93 && a->vof(i+d_i,j+d_j,k+d_k)<p->F93)
                    roval=p->W3;
                else
                {
                    if(d_i==1)
                        roval=(0.5*p->DXN[IP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                        +0.5*p->DXN[IP1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                        )/p->DXP[IP];
                    else if(d_i==-1)
                        roval=(0.5*p->DXN[IP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                        +0.5*p->DXN[IM1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                        )/p->DXP[IM1];
                    else if(d_k==1)
                        roval=(0.5*p->DZN[KP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                        +0.5*p->DZN[KP1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                        )/p->DZP[KP];
                    else if(d_k==-1)
                        roval=(0.5*p->DZN[KP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                        +0.5*p->DZN[KM1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                        )/p->DZP[KM1];
                    else if(d_j==1)
                        roval=(0.5*p->DYN[JP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                        +0.5*p->DYN[JP1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                        )/p->DYP[JP];
                    else if(d_j==-1)
                        roval=(0.5*p->DYN[JP]*(a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                        +0.5*p->DYN[JM1]*(a->vof(i+d_i,j+d_j,k+d_k)*p->W1+(1.0-a->vof(i+d_i,j+d_j,k+d_k))*p->W3)
                        )/p->DYP[JM1];
                    else
                        cout<<"density case missing d_i:"<<d_i<<" d_j:"<<d_j<<" d_k:"<<d_k<<" exepcheck 1"<<endl;
                }
            }
        }
                

        if(excepcheck==0)
        {
            if(p->j_dir<1)
            {
                if(d_i==1)
                    H=(0.25*p->DXN[IP]*a->vof_nt(i,j,k)+0.25*p->DXN[IP]*a->vof_nb(i,j,k)+0.25*p->DXN[IP1]*a->vof_st(i+1,j,k)+0.25*p->DXN[IP1]*a->vof_sb(i+1,j,k))/p->DXP[IP];
                else if(d_i==-1)
                    H=(0.25*p->DXN[IP]*a->vof_st(i,j,k)+0.25*p->DXN[IP]*a->vof_sb(i,j,k)+0.25*p->DXN[IM1]*a->vof_nt(i-1,j,k)+0.25*p->DXN[IM1]*a->vof_nb(i-1,j,k))/p->DXP[IM1];
                else if(d_j==1)
                    H=a->vof(i,j,k);
                else if(d_j==-1)
                    H=a->vof(i,j,k);
                else if(d_k==1)
                    H=(0.25*p->DZN[KP]*a->vof_nt(i,j,k)+0.25*p->DZN[KP]*a->vof_st(i,j,k)+0.25*p->DZN[KP1]*a->vof_nb(i,j,k+1)+0.25*p->DZN[KP1]*a->vof_sb(i,j,k+1))/p->DZP[KP];
                else if(d_k==-1)
                    H=(0.25*p->DZN[KP]*a->vof_nb(i,j,k)+0.25*p->DZN[KP]*a->vof_sb(i,j,k)+0.25*p->DZN[KM1]*a->vof_nt(i,j,k-1)+0.25*p->DZN[KM1]*a->vof_st(i,j,k-1))/p->DZP[KM1];
                else
                    cout<<"density case missing 2D, d_i:"<<d_i<<" d_j:"<<d_j<<" d_k:"<<d_k<<" exepcheck 0"<<endl;;
            }
            else
            {
                if(d_i==1)
                {
                    H=(0.5*p->DXN[IP]*(0.25*a->vof_nte(i,j,k)+0.25*a->vof_ntw(i,j,k)+0.25*a->vof_nbe(i,j,k)+0.25*a->vof_nbw(i,j,k))
                        +0.5*p->DXN[IP1]*(0.25*a->vof_ste(i+1,j,k)+0.25*a->vof_stw(i+1,j,k)+0.25*a->vof_sbe(i+1,j,k)+0.25*a->vof_sbw(i+1,j,k))
                        )/p->DXP[IP];
                }
                else if(d_i==-1)
                {
                    H=(0.5*p->DXN[IP]*(0.25*a->vof_ste(i,j,k)+0.25*a->vof_stw(i,j,k)+0.25*a->vof_sbe(i,j,k)+0.25*a->vof_sbw(i,j,k))
                        +0.5*p->DXN[IM1]*(0.25*a->vof_nte(i-1,j,k)+0.25*a->vof_ntw(i-1,j,k)+0.25*a->vof_nbe(i-1,j,k)+0.25*a->vof_nbw(i-1,j,k))
                        )/p->DXP[IM1];
                }
                else if(d_j==1)
                {
                    H=(0.5*p->DYN[JP]*(0.25*a->vof_ntw(i,j,k)+0.25*a->vof_nbw(i,j,k)+0.25*a->vof_stw(i,j,k)+0.25*a->vof_sbw(i,j,k))
                        +0.5*p->DYN[JP1]*(0.25*a->vof_nte(i,j+1,k)+0.25*a->vof_nbe(i,j+1,k)+0.25*a->vof_ste(i,j+1,k)+0.25*a->vof_sbe(i,j+1,k))
                        )/p->DYP[JP];
                }
                else if(d_j==-1)
                {
                    H=(0.5*p->DYN[JP]*(0.25*a->vof_nte(i,j,k)+0.25*a->vof_nbe(i,j,k)+0.25*a->vof_ste(i,j,k)+0.25*a->vof_sbe(i,j,k))
                        +0.5*p->DYN[JM1]*(0.25*a->vof_ntw(i,j-1,k)+0.25*a->vof_nbw(i,j-1,k)+0.25*a->vof_stw(i,j-1,k)+0.25*a->vof_sbw(i,j-1,k))
                        )/p->DYP[JM1];
                }
                else if(d_k==1)
                {
                    H=(0.5*p->DZN[KP]*(0.25*a->vof_nte(i,j,k)+0.25*a->vof_ntw(i,j,k)+0.25*a->vof_ste(i,j,k)+0.25*a->vof_stw(i,j,k))
                        +0.5*p->DZN[KP1]*(0.25*a->vof_nbe(i,j,k+1)+0.25*a->vof_nbw(i,j,k+1)+0.25*a->vof_sbe(i,j,k+1)+0.25*a->vof_sbw(i,j,k+1))
                        )/p->DZP[KP];
                }
                else if(d_k==-1)
                {
                    H=(0.5*p->DZN[KP]*(0.25*a->vof_nbe(i,j,k)+0.25*a->vof_nbw(i,j,k)+0.25*a->vof_sbe(i,j,k)+0.25*a->vof_sbw(i,j,k))
                        +0.5*p->DZN[KM1]*(0.25*a->vof_nte(i,j,k-1)+0.25*a->vof_ntw(i,j,k-1)+0.25*a->vof_ste(i,j,k-1)+0.25*a->vof_stw(i,j,k-1))
                        )/p->DZP[KM1];
                }
                else
                    cout<<"density case missing 3D d_i:"<<d_i<<" d_j:"<<d_j<<" d_k:"<<d_k<<" exepcheck 0"<<endl;
            }
        }
        roval=roval = p->W1*H + p->W3*(1.0-H);
    }
    
    return roval;
}




