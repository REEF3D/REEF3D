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
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"lexer.h"
#include"gradient.h"

void VOF_PLIC::calculateSubFractions(lexer* p, fdm* a, ghostcell* pgc, field& voffield)
{
    LOOP
    {
        if(a->phasemarker(i,j,k)>-1E-06 && a->phasemarker(i,j,k)<1E-06)
        {
            a->vof_nt(i,j,k)=0.0;
            a->vof_nb(i,j,k)=0.0;
            a->vof_st(i,j,k)=0.0;
            a->vof_sb(i,j,k)=0.0;
        }
        else if(a->phasemarker(i,j,k)>10.0 -1E-06 && a->phasemarker(i,j,k)<10.0+ 1E-06)
        {
            a->vof_nt(i,j,k)=1.0;
            a->vof_nb(i,j,k)=1.0;
            a->vof_st(i,j,k)=1.0;
            a->vof_sb(i,j,k)=1.0;
        }
        else
        {   
            double r0xz,scaledVol,recheck;
            reconstructPlane_alt(a,p,voffield);
            
            //NT
            r0xz=-(nx(i,j,k)*(0.25*p->DXN[IP])+nz(i,j,k)*(0.25*p->DZN[KP])-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                a->vof_nt(i,j,k)=scaledVol;
            }
            else
            {
                
                    
                if((nx(i,j,k)>0.0 && nz(i,j,k)<0.0) || (nx(i,j,k)<0.0 && nz(i,j,k)>0.0))
                {
                    if(alpha(i,j,k)>0.0)
                        a->vof_nt(i,j,k)=1.0;
                    else
                        a->vof_nt(i,j,k)=0.0;
                }
                else if(nx(i,j,k)>=0.0 && nz(i,j,k)>=0.0)
                    a->vof_nt(i,j,k)=0.0;
                else
                    a->vof_nt(i,j,k)=1.0;
            }
            
            //NB
            r0xz=-(nx(i,j,k)*(0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                a->vof_nb(i,j,k)=scaledVol;
            }
            else
            {
                if((nx(i,j,k)>0.0 && nz(i,j,k)>0.0) || (nx(i,j,k)<0.0 && nx(i,j,k)<0.0))
                {
                    if(alpha(i,j,k)>0.0)
                        a->vof_nb(i,j,k)=1.0;
                    else
                        a->vof_nb(i,j,k)=0.0;
                }
                else if(nx(i,j,k)>=0.0 && nz(i,j,k)<=0.0)
                    a->vof_nb(i,j,k)=0.0;
                else
                    a->vof_nb(i,j,k)=1.0;
            }
            
            //ST
            r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(0.25*p->DZN[KP])-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                a->vof_st(i,j,k)=scaledVol;
            }
            else
            {
                if((nx(i,j,k)>0.0 && nz(i,j,k)>0.0) || (nx(i,j,k)<0.0 && nx(i,j,k)<0.0))
                {
                    if(alpha(i,j,k)>0.0)
                        a->vof_st(i,j,k)=1.0;
                    else
                        a->vof_st(i,j,k)=0.0;
                }
                else if(nx(i,j,k)<=0.0 && nz(i,j,k)>=0.0)
                    a->vof_st(i,j,k)=0.0;
                else
                    a->vof_st(i,j,k)=1.0;
            }
            
            //SB
            r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                a->vof_sb(i,j,k)=scaledVol;
            }
            else
            {
                if((nx(i,j,k)>0.0 && nz(i,j,k)<0.0) || (nx(i,j,k)<0.0 && nz(i,j,k)>0.0))
                {
                    if(alpha(i,j,k)>0.0)
                        a->vof_sb(i,j,k)=1.0;
                    else
                        a->vof_sb(i,j,k)=0.0;
                }
                else if(nx(i,j,k)<=0.0 && nz(i,j,k)<=0.0)
                    a->vof_sb(i,j,k)=0.0;
                else
                    a->vof_sb(i,j,k)=1.0;
            }
        }
    }
    pgc->start4(p,a->vof_nt,1);
    pgc->start4(p,a->vof_nb,1);
    pgc->start4(p,a->vof_st,1);
    pgc->start4(p,a->vof_sb,1);
    
    return;
}