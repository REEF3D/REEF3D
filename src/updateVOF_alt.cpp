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
    double V_w_celltotal, V_a_celltotal;
	LOOP
    {
        V_w_celltotal=V_w_old(i,j,k)+V_w_update(i,j,k)+Watersafe(i,j,k);
        V_a_celltotal=V_a_old(i,j,k)+V_a_update(i,j,k);
        if(V_w_celltotal<0.0)
        {
          //  cout<<"water vol loss"<<endl;
            a->vof(i,j,k)=0.0;
        }
        else if(V_a_celltotal<0.0)
        {
           // cout<<"air vol loss"<<endl;
            a->vof(i,j,k)=1.0;
        } 
        else
        {
            if(V_w_celltotal+V_a_celltotal > p->DXN[IP]*p->DYN[JP]*p->DZN[KP])
         //       cout<<"overflow"<<endl;
         
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
    }
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