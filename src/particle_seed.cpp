/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>

void particle_f::seed_ini(lexer* p, fdm* a, ghostcell* pgc)
{
    // ini
    LOOP
    active(i,j,k) = 0.0;
    
    ppcell = 0;
    
    // Box
    cellcount=0;
    for(qn=0;qn<p->Q110;++qn)
    LOOP
	if(p->XN[IP]>=p->Q110_xs[qn] && p->XN[IP]<p->Q110_xe[qn]
	&& p->YN[JP]>=p->Q110_ys[qn] && p->YN[JP]<p->Q110_ye[qn]
	&& p->ZN[KP]>=p->Q110_zs[qn] && p->ZN[KP]<p->Q110_ze[qn])
	{
	active(i,j,k) = 1.0;
    ++cellcount;
	}
    
    // guess particle demand
    if(p->Q24>0)
    ppcell = p->Q24;
    
    partnum = cellcount * ppcell;
    
}

void particle_f::seed(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q110>0)
    posseed(p,a,pgc);
		

}


void particle_f::posseed(lexer* p, fdm* a, ghostcell* pgc)
{
	if(p->Q29>0)
    srand(p->Q29);

    if(p->Q29==0)
    srand((unsigned)time(0));
	
    LOOP
    if(active(i,j,k)>0.0)
    {
        
            for(qn=0;qn<ppcell;++qn)
            if(pactive<maxparticle)
            {
                pos[pactive][0] = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                pos[pactive][1] = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                pos[pactive][2] = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
                pos[pactive][3] = 0.0;
                posflag[pactive]=1;
                ++pactive;
            }
    }
    posactive=pactive;
    pcount=pactive;
}


void particle_f::posseed_topo(lexer* p, fdm* a, ghostcell* pgc)
{

        // // POS
        //     if(pcount>0)
        //     {
        //         reseeded++;
        //         pcount--;

        //         pos[PC][0] = (double(i) + (rand()%(irand))/drand)*dx;
        //         pos[PC][1] = (double(j) + (rand()%(irand))/drand)*dx;
        //         pos[PC][2] = (double(k) + (rand()%(irand))/drand)*dx;
        //         pos[PC][3] = phipol(p,a,pos[PC][0],pos[PC][1],pos[PC][2]);
        //         posflag[PC]=3;

        //         phival=MAX(((rand()%(irand))/drand)*epsi,rmin);

        //         lambda=1.0;
        //         qq=0;

        //         do
        //         {
        //         normal(a,pos[PC][0],pos[PC][1],pos[PC][2],pos[PC][3]);
        //         pos[PC][0] += lambda*(phival - pos[PC][3])*nvec[0];
        //         pos[PC][1] += lambda*(phival - pos[PC][3])*nvec[1];
        //         pos[PC][2] += lambda*(phival - pos[PC][3])*nvec[2];

        //         ii=int((pos[PC][0])/dx);
        //         jj=int((pos[PC][1])/dx);
        //         kk=int((pos[PC][2])/dx);
        //         check=boundcheck(p,a,ii,jj,kk,0);
        //         if(check==0)
        //         break;

        //         pos[PC][3] = phipol(p,a,pos[PC][0],pos[PC][1],pos[PC][2]);

        //         lambda/=2.0;
        //         ++qq;
        //         }while((pos[PC][3]>epsi || pos[PC][3]<rmin)&& qq<15);
				
		// 		//posradius(p,a,PC);

        //         if((pos[PC][3]>epsi || pos[PC][3]<rmin) || check==0)
        //         {
        //         posflag[PC]=0;
        //         pcount++;
        //         reseeded--;
        //         }
        //     }
			
        //     if(pcount==0 && posactive<maxparticle)
        //     {	
        //         pos[posactive][0] = (double(i)  + (rand()%(irand))/drand)*dx;
        //         pos[posactive][1] = (double(j)  + (rand()%(irand))/drand)*dx;
        //         pos[posactive][2] = (double(k)  + (rand()%(irand))/drand)*dx;
        //         pos[posactive][3] = phipol(p,a,pos[posactive][0],pos[posactive][1],pos[posactive][2]);
        //         posflag[posactive]=3;

        //         phival=MAX(((rand()%(irand))/drand)*epsi,rmin);

        //         lambda=1.0;
        //         qq=0;

        //         do
        //         {
        //         normal(a,pos[posactive][0],pos[posactive][1],pos[posactive][2],pos[posactive][3]);
        //         pos[posactive][0] += lambda*(phival - pos[posactive][3])*nvec[0];
        //         pos[posactive][1] += lambda*(phival - pos[posactive][3])*nvec[1];
        //         pos[posactive][2] += lambda*(phival - pos[posactive][3])*nvec[2];

        //         ii=int((pos[posactive][0])/dx);
        //         jj=int((pos[posactive][1])/dx);
        //         kk=int((pos[posactive][2])/dx);
        //         check=boundcheck(p,a,ii,jj,kk,0);
        //         if(check==0)
        //         break;

        //         pos[posactive][3] = phipol(p,a,pos[posactive][0],pos[posactive][1],pos[posactive][2]);
        //         lambda/=2.0;
        //         ++qq;
        //         }while((pos[posactive][3]>epsi || pos[posactive][3]<rmin) && qq<15);


        //         if(pos[posactive][3]<=epsi && pos[posactive][3]>=rmin && check==1)
        //         {
		// 		//posradius(p,a,posactive);
        //         posactive++;
        //         reseeded++;
        //         }
        //     }
			


}


