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
Author Fabian Knoblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void VOF_PLIC::surface_tension2D(lexer* p,fdm* a,ghostcell*pgc,int gcval)
{
    double calccurv, dHp, dHm, nx_loc, nz_loc;
    
	n=0;
	
	if(gcval==10 && p->W5>1.0e-10)
    {
        ULOOP
        {
            if(a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres)
            {
                calcNormalMYC3D(a,p,a->vof);
                if(nx(i,j,k)!=nx(i,j,k))
                    cout<<"nxNAN"<<endl;
                if(nz(i,j,k)!=nz(i,j,k))
                    cout<<"nzNAN"<<endl;    
                if(fabs(nz(i,j,k))>=fabs(nx(i,j,k)))
                {
                    dHp=((a->vof(i+1,j,k-1)*p->DZN[KM1]+a->vof(i+1,j,k)*p->DZN[KP]+a->vof(i+1,j,k+1)*p->DZN[KP1])
                        - (a->vof(i,j,k-1)*p->DZN[KM1]+a->vof(i,j,k)*p->DZN[KP]+a->vof(i,j,k+1)*p->DZN[KP1]))
                        /p->DXP[IP];
                        
                    dHm=((a->vof(i,j,k-1)*p->DZN[KM1]+a->vof(i,j,k)*p->DZN[KP]+a->vof(i,j,k+1)*p->DZN[KP1])
                        - (a->vof(i-1,j,k-1)*p->DZN[KM1]+a->vof(i-1,j,k)*p->DZN[KP]+a->vof(i-1,j,k+1)*p->DZN[KP1]))
                        /p->DXP[IM1];
                    
                    curv(i,j,k)=1.0/p->DXN[IP]*( dHp/pow((1.0+dHp*dHp),(3.0/2.0)) - dHm/pow((1.0+dHm*dHm),(3.0/2.0) ));
                }
                else
                {
                    dHp=((a->vof(i-1,j,k+1)*p->DXN[IM1]+a->vof(i,j,k+1)*p->DXN[IP]+a->vof(i+1,j,k+1)*p->DXN[IP1])
                       - (a->vof(i-1,j,k)*p->DXN[IM1]+a->vof(i,j,k)*p->DXN[IP]+a->vof(i+1,j,k)*p->DXN[IP1]))
                       /p->DZP[KP];
                       
                    dHm=((a->vof(i-1,j,k)*p->DXN[IM1]+a->vof(i,j,k)*p->DXN[IP]+a->vof(i+1,j,k)*p->DXN[IP1])
                        - (a->vof(i-1,j,k-1)*p->DXN[IM1]+a->vof(i,j,k-1)*p->DXN[IP]+a->vof(i+1,j,k-1)*p->DXN[IP1]))
                        /p->DZP[KM1];
                    
                    curv(i,j,k)=1.0/p->DZN[KP]*( dHp/pow((1.0+dHp*dHp),(3.0/2.0)) - dHm/pow((1.0+dHm*dHm),(3.0/2.0) ));
                }
            }
            else
                curv(i,j,k)=0.0;
        }
        pgc->start4(p,curv,1);
        pgc->start4(p,nx,1);
        pgc->start4(p,nz,1);
        ULOOP
        {   if((a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres) || (a->vof(i+1,j,k)<=w_thres && a->vof(i+1,j,k)>=a_thres))
            {
                if((a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres) && (a->vof(i+1,j,k)<=w_thres && a->vof(i+1,j,k)>=a_thres))
                {
                    calccurv=0.5*(curv(i,j,k)+curv(i+1,j,k));
                    nx_loc=0.5*(nx(i,j,k)+nx(i+1,j,k));
                }
                else if ((a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres) && ( !(a->vof(i+1,j,k)<=w_thres && a->vof(i+1,j,k)>=a_thres)))
                {
                    calccurv=curv(i,j,k);
                    nx_loc=nx(i,j,k);
                }
                else if ((!(a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres)) && (a->vof(i+1,j,k)<=w_thres && a->vof(i+1,j,k)>=a_thres))
                {
                    calccurv=curv(i+1,j,k);
                    nx_loc=nx(i+1,j,k);
                }
                else
                {
                    calccurv=0.0;
                    cout<<"bug_u_tens"<<endl;
                }
                if(nx_loc != nx_loc)
                    cout<<"nxLoc_NAN"<<endl;
                a->rhsvec.V[n]-=(calccurv*p->W5*(-nx_loc))/(0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)));
            }
            ++n;
        }
    }

	if(gcval==11 && p->W5>1.0e-10)
    {
       
    }

	if(gcval==12 && p->W5>1.0e-10)
    {
        WLOOP
        {
            if(a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres)
            {
                calcNormalMYC3D(a,p,a->vof);
                if(fabs(nz(i,j,k))>=fabs(nx(i,j,k)))
                {
                    dHp=((a->vof(i+1,j,k-1)*p->DZN[KM1]+a->vof(i+1,j,k)*p->DZN[KP]+a->vof(i+1,j,k+1)*p->DZN[KP1])
                        - (a->vof(i,j,k-1)*p->DZN[KM1]+a->vof(i,j,k)*p->DZN[KP]+a->vof(i,j,k+1)*p->DZN[KP1]))
                        /p->DXP[IP];
                        
                    dHm=((a->vof(i,j,k-1)*p->DZN[KM1]+a->vof(i,j,k)*p->DZN[KP]+a->vof(i,j,k+1)*p->DZN[KP1])
                        - (a->vof(i-1,j,k-1)*p->DZN[KM1]+a->vof(i-1,j,k)*p->DZN[KP]+a->vof(i-1,j,k+1)*p->DZN[KP1]))
                        /p->DXP[IM1];
                    
                    curv(i,j,k)=1.0/p->DXN[IP]*( dHp/pow((1.0+dHp*dHp),(3.0/2.0)) - dHm/pow((1.0+dHm*dHm),(3.0/2.0) ));
                }
                else
                {
                    dHp=((a->vof(i-1,j,k+1)*p->DXN[IM1]+a->vof(i,j,k+1)*p->DXN[IP]+a->vof(i+1,j,k+1)*p->DXN[IP1])
                       - (a->vof(i-1,j,k)*p->DXN[IM1]+a->vof(i,j,k)*p->DXN[IP]+a->vof(i+1,j,k)*p->DXN[IP1]))
                       /p->DZP[KP];
                       
                    dHm=((a->vof(i-1,j,k)*p->DXN[IM1]+a->vof(i,j,k)*p->DXN[IP]+a->vof(i+1,j,k)*p->DXN[IP1])
                        - (a->vof(i-1,j,k-1)*p->DXN[IM1]+a->vof(i,j,k-1)*p->DXN[IP]+a->vof(i+1,j,k-1)*p->DXN[IP1]))
                        /p->DZP[KM1];
                        
                    curv(i,j,k)=1.0/p->DZN[KP]*( dHp/pow((1.0+dHp*dHp),(3.0/2.0)) - dHm/pow((1.0+dHm*dHm),(3.0/2.0) ));
                }
            }
            else
                curv(i,j,k)=0.0;
        }
        pgc->start4(p,curv,1);
        WLOOP
        {   if((a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres) || (a->vof(i,j,k+1)<=w_thres && a->vof(i,j,k+1)>=a_thres))
            {
                if((a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres) && (a->vof(i,j,k+1)<=w_thres && a->vof(i,j,k+1)>=a_thres))
                {
                    calccurv=0.5*(curv(i,j,k)+curv(i,j,k+1));
                    nz_loc=0.5*(nz(i,j,k)+nz(i,j,k+1));
                }
                else if((a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres) && (!(a->vof(i,j,k+1)<=w_thres && a->vof(i,j,k+1)>=a_thres)))
                {
                    calccurv=curv(i,j,k);
                    nz_loc=nz(i,j,k);
                }
                else if((!(a->vof(i,j,k)<=w_thres && a->vof(i,j,k)>=a_thres)) && (a->vof(i,j,k+1)<=w_thres && a->vof(i,j,k+1)>=a_thres))
                {
                    calccurv=curv(i,j,k+1);
                    nz_loc=nz(i,j,k+1);
                }
                else
                {
                    calccurv=0.0;
                    cout<<"bug_w_tens"<<endl;
                }
                if(nz_loc != nz_loc)
                    cout<<"nzLoc_NAN"<<endl;
                a->rhsvec.V[n]-=(calccurv*p->W5*(-nz_loc))/(0.5*(a->ro(i,j,k)+a->ro(k+1,j,k)));
            }
            ++n;
        }
    }

}

