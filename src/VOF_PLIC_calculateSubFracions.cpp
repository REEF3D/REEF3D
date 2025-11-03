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
    if(p->j_dir<1)
    {
        LOOP
        {
            if(a->vof(i,j,k)<a_thres)
            {
                a->vof_nt(i,j,k)=0.0;
                a->vof_nb(i,j,k)=0.0;
                a->vof_st(i,j,k)=0.0;
                a->vof_sb(i,j,k)=0.0;
            }
            else if(a->vof(i,j,k)>w_thres)
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
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_nt(i,j,k)=scaledVol;
                }
                else
                {
                
                    
                    //if((nx(i,j,k)>0.0 && nz(i,j,k)<0.0) || (nx(i,j,k)<0.0 && nz(i,j,k)>0.0))
                   // {
                        if(alpha(i,j,k)>0.0)
                            a->vof_nt(i,j,k)=1.0;
                        else
                            a->vof_nt(i,j,k)=0.0;
                   /* }
                    else if(nx(i,j,k)>=0.0 && nz(i,j,k)>=0.0)
                        a->vof_nt(i,j,k)=0.0;
                    else
                        a->vof_nt(i,j,k)=1.0;*/
                }
            
                //NB
                r0xz=-(nx(i,j,k)*(0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_nb(i,j,k)=scaledVol;
                }
                else
                {
                   // if((nx(i,j,k)>0.0 && nz(i,j,k)>0.0) || (nx(i,j,k)<0.0 && nx(i,j,k)<0.0))
                   // {
                        if(alpha(i,j,k)>0.0)
                            a->vof_nb(i,j,k)=1.0;
                        else
                            a->vof_nb(i,j,k)=0.0;
                 /*   }
                    else if(nx(i,j,k)>=0.0 && nz(i,j,k)<=0.0)
                        a->vof_nb(i,j,k)=0.0;
                    else
                        a->vof_nb(i,j,k)=1.0; */
                }
            
                //ST
                r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(0.25*p->DZN[KP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_st(i,j,k)=scaledVol;
                }
                else
                {
                   // if((nx(i,j,k)>0.0 && nz(i,j,k)>0.0) || (nx(i,j,k)<0.0 && nx(i,j,k)<0.0))
                   // {
                        if(alpha(i,j,k)>0.0)
                            a->vof_st(i,j,k)=1.0;
                        else
                            a->vof_st(i,j,k)=0.0;
                  /*  }
                    else if(nx(i,j,k)<=0.0 && nz(i,j,k)>=0.0)
                        a->vof_st(i,j,k)=0.0;
                    else
                        a->vof_st(i,j,k)=1.0;*/
                }
            
                //SB
                r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_sb(i,j,k)=scaledVol;
                }
                else
                {
                   // if((nx(i,j,k)>0.0 && nz(i,j,k)<0.0) || (nx(i,j,k)<0.0 && nz(i,j,k)>0.0))
                   // {
                        if(alpha(i,j,k)>0.0)
                            a->vof_sb(i,j,k)=1.0;
                        else
                            a->vof_sb(i,j,k)=0.0;
                   /* }
                    else if(nx(i,j,k)<=0.0 && nz(i,j,k)<=0.0)
                        a->vof_sb(i,j,k)=0.0;
                    else
                        a->vof_sb(i,j,k)=1.0;*/
                }
            }
        }
        pgc->start4(p,a->vof_nt,1);
        pgc->start4(p,a->vof_nb,1);
        pgc->start4(p,a->vof_st,1);
        pgc->start4(p,a->vof_sb,1);
    }
    else
    {
        LOOP
        {
            if(a->vof(i,j,k)<a_thres)
            {
                a->vof_nte(i,j,k)=0.0;
                a->vof_ntw(i,j,k)=0.0;
                a->vof_nbe(i,j,k)=0.0;
                a->vof_nbw(i,j,k)=0.0;
                a->vof_ste(i,j,k)=0.0;
                a->vof_stw(i,j,k)=0.0;
                a->vof_sbe(i,j,k)=0.0;
                a->vof_sbw(i,j,k)=0.0;
            }
            else if(a->vof(i,j,k)>w_thres)
            {
                a->vof_nte(i,j,k)=1.0;
                a->vof_ntw(i,j,k)=1.0;
                a->vof_nbe(i,j,k)=1.0;
                a->vof_nbw(i,j,k)=1.0;
                a->vof_ste(i,j,k)=1.0;
                a->vof_stw(i,j,k)=1.0;
                a->vof_sbe(i,j,k)=1.0;
                a->vof_sbw(i,j,k)=1.0;
            }
            else
            {   
                double r0xz,scaledVol,recheck;
                reconstructPlane_alt(a,p,voffield);
            
                //NTE
                r0xz=-(nx(i,j,k)*(0.25*p->DXN[IP])+nz(i,j,k)*(0.25*p->DZN[KP])+ny(i,j,k)*(-0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_nte(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_nte(i,j,k)=1.0;
                        else
                            a->vof_nte(i,j,k)=0.0;
                }
                
                //NTW
                r0xz=-(nx(i,j,k)*(0.25*p->DXN[IP])+nz(i,j,k)*(0.25*p->DZN[KP])+ny(i,j,k)*(0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_ntw(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_ntw(i,j,k)=1.0;
                        else
                            a->vof_ntw(i,j,k)=0.0;
                }
                
                //NBE
                r0xz=-(nx(i,j,k)*(0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])+ny(i,j,k)*(-0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_nbe(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_nbe(i,j,k)=1.0;
                        else
                            a->vof_nbe(i,j,k)=0.0;
                }
                
                //NBW
                r0xz=-(nx(i,j,k)*(0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])+ny(i,j,k)*(0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_nbw(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_nbw(i,j,k)=1.0;
                        else
                            a->vof_nbw(i,j,k)=0.0;
                }
                
                //STE
                r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(0.25*p->DZN[KP])+ny(i,j,k)*(-0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_ste(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_ste(i,j,k)=1.0;
                        else
                            a->vof_ste(i,j,k)=0.0;
                }
                
                //STW
                r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(0.25*p->DZN[KP])+ny(i,j,k)*(0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_stw(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_stw(i,j,k)=1.0;
                        else
                            a->vof_stw(i,j,k)=0.0;
                }
                
                //SBE
                r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])+ny(i,j,k)*(-0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_sbe(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_sbe(i,j,k)=1.0;
                        else
                            a->vof_sbe(i,j,k)=0.0;
                }
                
                //SBW
                r0xz=-(nx(i,j,k)*(-0.25*p->DXN[IP])+nz(i,j,k)*(-0.25*p->DZN[KP])+ny(i,j,k)*(0.25*p->DYN[JP])-alpha(i,j,k));
                recheck=0.5*(fabs(nx(i,j,k))*0.5*p->DXN[IP]+fabs(ny(i,j,k))*0.5*p->DYN[JP]+fabs(nz(i,j,k))*0.5*p->DZN[KP])-fabs(r0xz);
                if(recheck>1E-20)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),0.5*p->DXN[IP],0.5*p->DYN[JP],0.5*p->DZN[KP],r0xz);
                    a->vof_sbw(i,j,k)=scaledVol;
                }
                else
                {
                        if(alpha(i,j,k)>0.0)
                            a->vof_sbw(i,j,k)=1.0;
                        else
                            a->vof_sbw(i,j,k)=0.0;
                }
                    
                
                
                
            }   
        }
        pgc->start4(p,a->vof_nte,1);
        pgc->start4(p,a->vof_ntw,1);
        pgc->start4(p,a->vof_nbe,1);
        pgc->start4(p,a->vof_nbw,1);
        pgc->start4(p,a->vof_ste,1);
        pgc->start4(p,a->vof_stw,1);
        pgc->start4(p,a->vof_sbe,1);
        pgc->start4(p,a->vof_sbw,1);

    }
}