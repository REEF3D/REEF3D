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

void VOF_PLIC::updatePhasemarkers( lexer* p, fdm* a, ghostcell* pgc,field& voffield)
{
    LOOP
        a->phasemarker(i,j,k)=0.0;
        
    LOOP
    {
        if(voffield(i,j,k)>w_thres)
        {
            a->phasemarker(i,j,k)=10.0;
        }
    }
    
    pgc->start4(p,a->phasemarker,1);
    
    LOOP
    {
        if(a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<1.0)
        {   
            if(searchMarkerInVicinity(p,a,1,10.0,i,j,k)>=1)
            {
                a->phasemarker(i,j,k)=6.0;
            }
        }
    }
    
    pgc->start4(p,a->phasemarker,1);
    /*
    LOOP
    {
        if(a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<0.1)
        {
            if((searchMarkerAlongDims(p,a,1,6.0,i,j,k)>=2) && (searchMarkerInVicinity(p,a,1,10.0,i,j,k)>=1))
            {
                a->phasemarker(i,j,k)=4.0;
            }
        }
    }
    
    pgc->start4(p,a->phasemarker,1);*/

} 

void VOF_PLIC::updatePhasemarkersCompression( lexer* p, fdm* a, ghostcell* pgc,field& voffield)
{
    LOOP
    {
        a->phasemarker(i,j,k)=0.0;
        compressvol(i,j,k)=0.0;
    }
    LOOP
    {
        if(voffield(i,j,k)>corr_thres)
        {
            a->phasemarker(i,j,k)=10.0;
        }
    }
    
    pgc->start4(p,a->phasemarker,1);
    
    LOOP
    {
        if(a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<1.0)
        {   
            if(searchMarkerInVicinity(p,a,2,10.0,i,j,k)>=1)
            {
                a->phasemarker(i,j,k)=6.0;
            }
        }
    }
    
    pgc->start4(p,a->phasemarker,1);
    pgc->start4(p,compressvol,1);
    
    LOOP
    {
        if( (voffield(i,j,k)>=a_thres && voffield(i,j,k) <= w_thres) && (a->phasemarker(i,j,k) > -0.1 && a->phasemarker(i,j,k) < 0.1) )
        {
            switch(p->F88)
            {
                case 10:
                    calcNormalMYC2D(a,p,voffield);
                    break;
                case 11:
                    calcNormalMYC2D_V2(a,p,voffield);
                    break;
                case 12:
                    calcNormalMYC2D_V3(a,p,voffield);
                    break;
                case 13:
                    calcNormalMYC2D_V4(a,p,voffield);
                    break;
            }
            /*if(fabs(nz(i,j,k))>=fabs(nx(i,j,k)))
            {
                if(nz(i,j,k)>=0.0)
                    compressvol(i,j,k-1)+=voffield(i,j,k);
                else
                    compressvol(i,j,k+1)+=voffield(i,j,k);
            }
            else
            {
                if(nx(i,j,k)>=0.0)
                    compressvol(i-1,j,k)+=voffield(i,j,k);
                else
                    compressvol(i+1,j,k)+=voffield(i,j,k);
            }*/
            
            if(nz(i,j,k)>=0.0)
            {
                if(nx(i,j,k)>=0.0)
                {
                    compressvol(i,j,k-1)+=fabs(nz(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                    compressvol(i-1,j,k)+=fabs(nx(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                }
                else
                {
                    compressvol(i,j,k-1)+=fabs(nz(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                    compressvol(i+1,j,k)+=fabs(nx(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                }
            }
            else
            {
                 if(nx(i,j,k)>=0.0)
                {
                    compressvol(i,j,k+1)+=fabs(nz(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                    compressvol(i-1,j,k)+=fabs(nx(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                }
                else
                {
                    compressvol(i,j,k+1)+=fabs(nz(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                    compressvol(i+1,j,k)+=fabs(nx(i,j,k))/(fabs(nx(i,j,k))+fabs(nz(i,j,k)))*voffield(i,j,k);
                }
            }
            
            a->phasemarker(i,j,k)=5.0;
        }
    }
    pgc->start4(p,a->phasemarker,1);
    pgc->start4(p,compressvol,1);
    
    LOOP
    {
        voffield(i,j,k)=voffield(i,j,k)+compressvol(i,j,k);
        if(a->phasemarker(i,j,k)>4.9 && a->phasemarker(i,j,k)<5.1)
            voffield(i,j,k)=0.0;
    }
    pgc->start4(p,a->phasemarker,1);
    pgc->start4(p,voffield,1);
    
    LOOP
    {
        a->phasemarker(i,j,k)=0.0;
        compressvol(i,j,k)=0.0;
    }
    LOOP
    {
        if(voffield(i,j,k)>corr_thres)
        {
            a->phasemarker(i,j,k)=10.0;
        }
    }
    
    pgc->start4(p,a->phasemarker,1);
    
    LOOP
    {
        if(a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<1.0)
        {   
            if(searchMarkerInVicinity(p,a,2,10.0,i,j,k)>=1)
            {
                a->phasemarker(i,j,k)=6.0;
            }
        }
    }
    
    LOOP
    {
        if(a->phasemarker(i,j,k)>9.9 && a->phasemarker(i,j,k)<10.1)
            voffield(i,j,k)=1.0;
            
        if((a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<0.1))
            voffield(i,j,k)=0.0;
    }
    
    pgc->start4(p,voffield,1); 
    
    return;
    
}

void VOF_PLIC::updatePhasemarkersCorrection( lexer* p, fdm* a, ghostcell* pgc,field& voffield)
{
    LOOP
        a->phasemarker(i,j,k)=0.0;
        
    LOOP
    {
        if(voffield(i,j,k)>corr_thres)
        {
            a->phasemarker(i,j,k)=10.0;
        }
    }
    
    pgc->start4(p,a->phasemarker,1);
    
    LOOP
    {
        if(a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<0.1)
        {   
            if(searchMarkerInVicinity(p,a,12,10.0,i,j,k)>=1)
            {
                a->phasemarker(i,j,k)=6.0;
            }
        }
    }
    
    pgc->start4(p,a->phasemarker,1);
    LOOP
    {
        if(a->phasemarker(i,j,k)>9.9 && a->phasemarker(i,j,k)<10.1)
        if(searchMarkerInVicinity(p,a,1,6.0,i,j,k)>=1)
        {
            a->phasemarker(i,j,k)=4.0;
        }
        
    }
    pgc->start4(p,a->phasemarker,1);
    
    /*LOOP
    {
        if(a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<0.1)
        {
            if((searchMarkerAlongDims(p,a,1,6.0,i,j,k)>=2) && (searchMarkerInVicinity(p,a,1,10.0,i,j,k)>=1))
            {
                a->phasemarker(i,j,k)=4.0;
            }
        }
    }
    
    pgc->start4(p,a->phasemarker,1);*/
    
    LOOP
    {
        if(a->phasemarker(i,j,k)>9.9 && a->phasemarker(i,j,k)<10.1)
            voffield(i,j,k)=1.0;
            
        if((a->phasemarker(i,j,k)>-0.1 && a->phasemarker(i,j,k)<0.1))
            voffield(i,j,k)=0.0;
    }
    
    pgc->start4(p,voffield,1); 
    
    return;
}

int VOF_PLIC::searchMarkerInVicinity(lexer*p, fdm* a, int dist, double markernum,int ii, int jj, int kk)
{
    int returnnum = 0;
    
    if(p->j_dir>0)
    {
        for (int i_s = -dist; i_s < dist+1 ; i_s++)
        {
            for(int j_s = -dist; j_s < dist+1 ; j_s++)
            {
                for(int k_s= -dist; k_s < dist+1 ; k_s++)
                {
                    if(a->phasemarker(ii+i_s,jj+j_s,kk+k_s) > markernum -0.1 && a->phasemarker(ii+i_s,jj+j_s,kk+k_s) < markernum +0.1)
                        returnnum++;
                }
            }
        }
    }
    else
    {
        for (int i_s = -dist; i_s < dist+1 ; i_s++)
        {
            
            for(int k_s= -dist; k_s < dist+1 ; k_s++)
            {
                if(a->phasemarker(ii+i_s,jj,kk+k_s) > markernum -0.1 && a->phasemarker(ii+i_s,jj,kk+k_s) < markernum +0.1)
                    returnnum++;
            }
        }
    }
    
    return returnnum;
    
}

int VOF_PLIC::searchMarkerAlongDims(lexer* p, fdm* a, int dist, double markernum, int ii, int jj, int kk)
{
    int returnnum = 0;
    
    for(int i_s = -dist; i_s < dist+1; i_s++)
    {
        if(a->phasemarker(ii+i_s,jj,kk) > markernum -0.1 && a->phasemarker(ii+i_s,jj,kk) < markernum +0.1)
            returnnum++;
    }
    
    if(p->j_dir>0)
    {
        for(int j_s = -dist; j_s < dist+1; j_s++)
        {
            if(a->phasemarker(ii,jj+j_s,kk) > markernum -0.1 && a->phasemarker(ii,jj+j_s,kk) < markernum +0.1)
                returnnum++;
        }
    }
    
    for(int k_s = -dist; k_s < dist+1; k_s++)
    {
        if(a->phasemarker(ii,jj,kk+k_s) > markernum -0.1 && a->phasemarker(ii,jj,kk+k_s) < markernum +0.1)
            returnnum++;
    }
    
    return returnnum;
}

