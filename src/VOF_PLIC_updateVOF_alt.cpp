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

void VOF_PLIC::updateVOF_alt(fdm* a, lexer* p,int sweep)
{
    /*double V_w_celltotal, V_a_celltotal;
	LOOP
    {
        V_w_celltotal=V_w_old(i,j,k)+V_w_update(i,j,k)+Watersafe(i,j,k);
        if(V_w_celltotal<0.0)
        {
          //  cout<<"water vol loss"<<endl;
            a->vof(i,j,k)=0.0;
        }
        else
        {
            if(V_w_celltotal>(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]))
            {
                Watersafe(i,j,k)=V_w_celltotal-(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                V_w_celltotal -=Watersafe(i,j,k);
            }
            
            a->vof(i,j,k)=V_w_celltotal/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            
            if(a->vof(i,j,k)>1.0)
            {
               // cout<<"water overflow"<<endl;
                a->vof(i,j,k)=1.0;
            }
        }
    }*/
    /*
    LOOP
    {
        if(a->vof(i,j,k)>0.001 && a->vof(i,j,k)<0.999)
        {
            int checko = 0;
            for(int isch=-1; isch<2;isch++)
            {
                for(int jsch=-1; jsch<2;jsch++)
                {
                    for(int ksch=-1; ksch<2;ksch++)
                    {
                        if(a->vof(i+isch,j+jsch,k+ksch)>0.999)
                            checko++;
                    }
                }
            }
            
            if(checko<1)
            {
                if(sweep==0)
                {
                    if(a->vof(i-2,j,k)>0.999)
                    {
                        a->vof(i-1,j,k)+=a->vof(i,j,k);
                        a->vof(i,j,k)=0.0;
                        if(a->vof(i-1,j,k)>1.0)
                        {
                            Watersafe(i-1,j,k)=(a->vof(i-1,j,k)-1.0)*p->DXN[IM1]*p->DYN[IM1]*p->DZN[IM1];
                            a->vof(i-1,j,k)=1.0;
                        }
                    }
                    else if(a->vof(i+2,j,k)>0.999)
                    {
                        a->vof(i+1,j,k)+=a->vof(i,j,k);
                        a->vof(i,j,k)=0.0;
                        if(a->vof(i+1,j,k)>1.0)
                        {
                            Watersafe(i+1,j,k)=(a->vof(i+1,j,k)-1.0)*p->DXN[IP1]*p->DYN[IP1]*p->DZN[IP1];
                            a->vof(i+1,j,k)=1.0;
                        }
                    }
                    else
                        cout<<"water lost along sweep:"<<sweep<<endl;
                }
                else if(sweep==1)
                {
                    if(a->vof(i,j-2,k)>0.999)
                    {
                        a->vof(i,j-1,k)+=a->vof(i,j,k);
                        a->vof(i,j,k)=0.0;
                        if(a->vof(i,j-1,k)>1.0)
                        {
                            Watersafe(i,j-1,k)=(a->vof(i,j-1,k)-1.0)*p->DXN[JM1]*p->DYN[JM1]*p->DZN[JM1];
                            a->vof(i,j-1,k)=1.0;
                        }
                    }
                    else if(a->vof(i,j+2,k)>0.999)
                    {
                        a->vof(i,j+1,k)+=a->vof(i,j,k);
                        a->vof(i,j,k)=0.0;
                        if(a->vof(i,j+1,k)>1.0)
                        {
                            Watersafe(i,j+1,k)=(a->vof(i,j+1,k)-1.0)*p->DXN[JP1]*p->DYN[JP1]*p->DZN[JP1];
                            a->vof(i,j+1,k)=1.0;
                        }
                    }
                    else
                        cout<<"water lost along sweep:"<<sweep<<endl;
                }
            }
            else
                {
                    if(a->vof(i,j,k-2)>0.999)
                    {
                        a->vof(i,j,k-1)+=a->vof(i,j,k);
                        a->vof(i,j,k)=0.0;
                        if(a->vof(i,j,k-1)>1.0)
                        {
                            Watersafe(i,j,k-1)=(a->vof(i,j,k-1)-1.0)*p->DXN[KM1]*p->DYN[KM1]*p->DZN[KM1];
                            a->vof(i,j-1,k)=1.0;
                        }
                    }
                    else if(a->vof(i,j,k+2)>0.999)
                    {
                        a->vof(i,j,k+1)+=a->vof(i,j,k);
                        a->vof(i,j,k)=0.0;
                        if(a->vof(i,j,k+1)>1.0)
                        {
                            Watersafe(i,j,k+1)=(a->vof(i,j,k+1)-1.0)*p->DXN[KP1]*p->DYN[KP1]*p->DZN[KP1];
                            a->vof(i,j,k+1)=1.0;
                        }
                    }
                    else
                        cout<<"water lost along sweep:"<<sweep<<endl;
                }
        }
    }*/
}

void VOF_PLIC::updateVOF_sweepless(fdm* a, lexer* p)
{
    double C_c, dudx,dwdz;
	LOOP
    {
        dudx=(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
        dwdz=(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
        if(vof_old(i,j,k)>0.5)
            C_c=0.0;
        else
            C_c=0.0;
            
        a->vof(i,j,k)=vof_old(i,j,k)+V_w_update(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP])+C_c*dudx*p->dt+C_c*dwdz*p->dt+Watersafe(i,j,k);
            if(a->vof(i,j,k)>1.0)
            {
                cout<<"overflow"<<endl;
                //Watersafe(i,j,k)=1.0-a->vof(i,j,k);
                a->vof(i,j,k)=1.0;
            }
            else if(a->vof(i,j,k)<0.0)
            {
                cout<<"underflow"<<endl;
               // Watersafe(i,j,k)=a->vof(i,j,k);
                a->vof(i,j,k)=0.0;
            }
            
        
        
    }
    /*
    LOOP
    {
        if((a->vof(i,j,k)>=0.001 && a->vof(i,j,k)<=0.999) && (!(a->vof(i+1,j,k)>0.999 || a->vof(i-1,j,k)>0.999 || a->vof(i,j,k+1)>0.999 || a->vof(i,j,k-1)>0.999 || a->vof(i+1,j,k+1)>0.999 || a->vof(i+1,j,k-1)>0.999 || a->vof(i-1,j,k+1)>0.999 || a->vof(i-1,j,k-1)>0.999)))
        {
            if(fabs(a->u(i,j,k))>fabs(a->w(i,j,k)))
            {
                if(a->u(i,j,k)>0.0)
                {
                    if(a->vof(i-1,j,k)>=0.001)
                    {
                        a->vof(i-1,j,k)+=a->vof(i,j,k)*p->DXN[IP]/p->DXN[IM1];
                        if(a->vof(i-1,j,k)>1.0)
                            a->vof(i-1,j,k)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else if(a->vof(i+1,j,k)>=0.001)
                    {
                        a->vof(i+1,j,k)+=a->vof(i,j,k)*p->DXN[IP]/p->DXN[IM1];
                        if(a->vof(i+1,j,k)>1.0)
                            a->vof(i+1,j,k)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else if(a->vof(i,j,k-1)>=0.001)
                    {
                        a->vof(i,j,k-1)+=a->vof(i,j,k)*p->DZN[KP]/p->DZN[KM1];
                        if(a->vof(i,j,k-1)>1.0)
                            a->vof(i,j,k-1)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else if(a->vof(i,j,k+1)>=0.001)
                    {
                        a->vof(i,j,k+1)+=a->vof(i,j,k)*p->DZN[KP]/p->DZN[KM1];
                        if(a->vof(i,j,k+1)>1.0)
                            a->vof(i,j,k+1)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else
                        cout<<"problem finding water"<<endl;
                }
                else
                {
                    if(a->vof(i+1,j,k)>=0.001)
                    {
                        a->vof(i+1,j,k)+=a->vof(i,j,k)*p->DXN[IP]/p->DXN[IM1];
                        if(a->vof(i+1,j,k)>1.0)
                            a->vof(i+1,j,k)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else if(a->vof(i-1,j,k)>=0.001)
                    {
                        a->vof(i-1,j,k)+=a->vof(i,j,k)*p->DXN[IP]/p->DXN[IM1];
                        if(a->vof(i-1,j,k)>1.0)
                            a->vof(i-1,j,k)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else if(a->vof(i,j,k-1)>=0.001)
                    {
                        a->vof(i,j,k-1)+=a->vof(i,j,k)*p->DZN[KP]/p->DZN[KM1];
                        if(a->vof(i,j,k-1)>1.0)
                            a->vof(i,j,k-1)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else if(a->vof(i,j,k+1)>=0.001)
                    {
                        a->vof(i,j,k+1)+=a->vof(i,j,k)*p->DZN[KP]/p->DZN[KM1];
                        if(a->vof(i,j,k+1)>1.0)
                            a->vof(i,j,k+1)=1.0;
                        a->vof(i,j,k)=0.0;
                    }
                    else
                        cout<<"problem finding water"<<endl;
                }
            }
            else
            {   if(a->vof(i,j,k-1)>=0.001)
                {
                    a->vof(i,j,k-1)+=a->vof(i,j,k)*p->DZN[KP]/p->DZN[KM1];
                    if(a->vof(i,j,k-1)>1.0)
                        a->vof(i,j,k-1)=1.0;
                    a->vof(i,j,k)=0.0;
                }
                else if(a->vof(i,j,k+1)>=0.001)
                {
                    a->vof(i,j,k+1)+=a->vof(i,j,k)*p->DZN[KP]/p->DZN[KM1];
                    if(a->vof(i,j,k+1)>1.0)
                    a->vof(i,j,k+1)=1.0;
                    a->vof(i,j,k)=0.0;
                }
                else if(a->vof(i-1,j,k)>=0.001)
                {
                    a->vof(i-1,j,k)+=a->vof(i,j,k)*p->DXN[IP]/p->DXN[IM1];
                    if(a->vof(i-1,j,k)>1.0)
                        a->vof(i-1,j,k)=1.0;
                    a->vof(i,j,k)=0.0;
                }
                else if(a->vof(i+1,j,k)>=0.001)
                {
                    a->vof(i+1,j,k)+=a->vof(i,j,k)*p->DXN[IP]/p->DXN[IM1];
                    if(a->vof(i+1,j,k)>1.0)
                        a->vof(i+1,j,k)=1.0;
                    a->vof(i,j,k)=0.0;
                }
                else
                        cout<<"problem finding water"<<endl;
                
            }
        }
    }*/
}

void VOF_PLIC::updateVOF_Weymouth
(
    fdm* a,
    lexer* p,
    int sweep
)
{ //2D 
    double C_c, dudi;
    if(sweep==0)
    {
        LOOP
        {
            if(vof_old(i,j,k)>0.5)
                C_c=1.0;
            else
                C_c=1.0;
            
            dudi=(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
            a->vof(i,j,k)=vof_old(i,j,k)+V_w_update(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP])+C_c*dudi*p->dt+Watersafe(i,j,k);
            if(a->vof(i,j,k)>1.0)
            {
                cout<<"overflow"<<endl;
                Watersafe(i,j,k)=1.0-a->vof(i,j,k);
                a->vof(i,j,k)=1.0;
            }
            else if(a->vof(i,j,k)<0.0)
            {
                cout<<"underflow"<<endl;
                Watersafe(i,j,k)=a->vof(i,j,k);
                a->vof(i,j,k)=0.0;
            }
        }
    }
    else if(sweep==2)
    {
        LOOP
        {
            if(vof_old(i,j,k)>0.5)
                C_c=1.0;
            else
                C_c=1.0;
            dudi=(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
            a->vof(i,j,k)=vof_old(i,j,k)+V_w_update(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP])+C_c*dudi*p->dt+Watersafe(i,j,k);
            if(a->vof(i,j,k)>1.0)
            {
                cout<<"overflow"<<endl;
                Watersafe(i,j,k)=1.0-a->vof(i,j,k);
                a->vof(i,j,k)=1.0;
            }
            else if(a->vof(i,j,k)<0.0)
            {
                cout<<"underflow"<<endl;
                Watersafe(i,j,k)=a->vof(i,j,k);
                a->vof(i,j,k)=0.0;
            }
        }
    }
}
                

void VOF_PLIC::updateVOF_MACHO2D
(
    fdm* a,
    lexer* p,
    int sweep,
    int nSweep
)
{   
    if(nSweep<2)
    {
        if(sweep==0)
        {
            LOOP
            {
                vof_old(i,j,k)=vof_old(i,j,k)
                                -(V_w_p(i,j,k)-V_w_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP])
                                +p->dt*vof_old(i,j,k)*(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
                V_w_p_star(i,j,k)=V_w_p(i,j,k);
                V_w_m_star(i,j,k)=V_w_m(i,j,k);
            }
        }
        else if(sweep==2)
        {
            LOOP
            {
                vof_old(i,j,k)=vof_old(i,j,k)
                                -(V_w_p(i,j,k)-V_w_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP])
                                +p->dt*vof_old(i,j,k)*(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
                V_w_p_star(i,j,k)=V_w_p(i,j,k);
                V_w_p_star(i,j,k)=V_w_p(i,j,k);
            }
        }
    }
    else
    {
        LOOP
        {
        vof_old(i,j,k)=vof_prevstep(i,j,k)
                        -(V_w_p(i,j,k)-V_w_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP])
                        -(V_w_p_star(i,j,k)-V_w_p_star(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
        a->vof(i,j,k)=vof_old(i,j,k);
        vof_prevstep(i,j,k)=vof_old(i,j,k);
        V_w_p_star(i,j,k)=0.0;
        V_w_p_star(i,j,k)=0.0;
                        
        }
    }
}

