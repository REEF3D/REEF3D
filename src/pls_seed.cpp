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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>

void particle_pls::seed(lexer* p, fdm* a, ghostcell* pgc, double fraction,double fac)
{
    field4 partcount(p);

    LOOP
    partcount(i,j,k)=0.0;

    reseeded=0;

    // POS
    for(n=0;n<posactive;++n)
    if(posflag[n]>0)
    {
        i=int((pos[n][0])/dx);
        j=int((pos[n][1])/dx);
        k=int((pos[n][2])/dx);

    partcount(i,j,k)+=1.0;
    }

    // NEG
    for(n=0;n<negactive;++n)
    if(negflag[n]>0)
    {
        i=int((neg[n][0])/dx);
        j=int((neg[n][1])/dx);
        k=int((neg[n][2])/dx);

    partcount(i,j,k)+=1.0;
    }

    LOOP
    if(fabs(a->phi(i,j,k))<dx*fac)
    {

        // POS
        if(partcount(i,j,k)<double(pnum)*fraction && a->phi(i,j,k)>0.0)
        {
			while(partcount(i,j,k)<double(pnum)*fraction)
			{
			check=posseed(p,a,pgc,1.0);
			
			if(check==1)
			partcount(i,j,k)+=1.0;
			}
			
		}	

        //NEG
        if(partcount(i,j,k)<double(pnum)*fraction && a->phi(i,j,k)<0.0)
        {

            while(partcount(i,j,k)<double(pnum)*fraction)
            {
			check=negseed(p,a,pgc,1.0);
			
			if(check==1)
			partcount(i,j,k)+=1.0;
			}
			
		}	
	}

}


int particle_pls::posseed(lexer* p, fdm* a, ghostcell* pgc, double factor)
{
	int success=1;
	
        // POS
            if(pcount>0)
            {
                reseeded++;
                pcount--;

                pos[posmem[pcount]][0] = (double(i) + (rand()%(irand))/drand)*dx;
                pos[posmem[pcount]][1] = (double(j) + (rand()%(irand))/drand)*dx;
                pos[posmem[pcount]][2] = (double(k) + (rand()%(irand))/drand)*dx;
                pos[posmem[pcount]][3] = phipol(p,a,pos[posmem[pcount]][0],pos[posmem[pcount]][1],pos[posmem[pcount]][2]);
                posflag[posmem[pcount]]=3;

                phival=MAX(((rand()%(irand))/drand)*epsi,rmin*factor);

                lambda=1.0;
                qq=0;

                do
                {
                normal(a,pos[posmem[pcount]][0],pos[posmem[pcount]][1],pos[posmem[pcount]][2],pos[posmem[pcount]][3]);
                pos[posmem[pcount]][0] += lambda*(phival - pos[posmem[pcount]][3])*nvec[0];
                pos[posmem[pcount]][1] += lambda*(phival - pos[posmem[pcount]][3])*nvec[1];
                pos[posmem[pcount]][2] += lambda*(phival - pos[posmem[pcount]][3])*nvec[2];

                ii=int((pos[posmem[pcount]][0])/dx);
                jj=int((pos[posmem[pcount]][1])/dx);
                kk=int((pos[posmem[pcount]][2])/dx);
                check=boundcheck(p,a,ii,jj,kk,0);
                if(check==0)
                break;

                pos[posmem[pcount]][3] = phipol(p,a,pos[posmem[pcount]][0],pos[posmem[pcount]][1],pos[posmem[pcount]][2]);

                lambda/=2.0;
                ++qq;
                }while((pos[posmem[pcount]][3]>epsi || pos[posmem[pcount]][3]<rmin)&& qq<15);
				
				//posradius(p,a,posmem[pcount]);

                if((pos[posmem[pcount]][3]>epsi || pos[posmem[pcount]][3]<rmin) || check==0)
                {
                posflag[posmem[pcount]]=0;
				 pcount++;
                reseeded--;
				success=0;
                }
            }
			
            if(pcount==0 && posactive<maxparticle)
            {	
                pos[posactive][0] = (double(i)  + (rand()%(irand))/drand)*dx;
                pos[posactive][1] = (double(j)  + (rand()%(irand))/drand)*dx;
                pos[posactive][2] = (double(k)  + (rand()%(irand))/drand)*dx;
                pos[posactive][3] = phipol(p,a,pos[posactive][0],pos[posactive][1],pos[posactive][2]);
                posflag[posactive]=3;

                phival=MAX(((rand()%(irand))/drand)*epsi,rmin*factor);

                lambda=1.0;
                qq=0;

                do
                {
                normal(a,pos[posactive][0],pos[posactive][1],pos[posactive][2],pos[posactive][3]);
                pos[posactive][0] += lambda*(phival - pos[posactive][3])*nvec[0];
                pos[posactive][1] += lambda*(phival - pos[posactive][3])*nvec[1];
                pos[posactive][2] += lambda*(phival - pos[posactive][3])*nvec[2];

                ii=int((pos[posactive][0])/dx);
                jj=int((pos[posactive][1])/dx);
                kk=int((pos[posactive][2])/dx);
                check=boundcheck(p,a,ii,jj,kk,0);
                if(check==0)
                break;

                pos[posactive][3] = phipol(p,a,pos[posactive][0],pos[posactive][1],pos[posactive][2]);
                lambda/=2.0;
                ++qq;
                }while((pos[posactive][3]>epsi || pos[posactive][3]<rmin) && qq<15);


                if(pos[posactive][3]<=epsi && pos[posactive][3]>=rmin && check==1)
                {
				//posradius(p,a,posactive);
                posactive++;
                reseeded++;
				success=1;
                }
            }
			
			
	return success;

}


int particle_pls::negseed(lexer* p, fdm* a, ghostcell* pgc, double factor)
{	
	int success=1;
	
	
            if(ncount>0)
            {
                reseeded++;
                ncount--;

                neg[negmem[ncount]][0] = (double(i)  + (rand()%(irand))/drand)*dx;
                neg[negmem[ncount]][1] = (double(j)  + (rand()%(irand))/drand)*dx;
                neg[negmem[ncount]][2] = (double(k)  + (rand()%(irand))/drand)*dx;
                neg[negmem[ncount]][3] = phipol(p,a,neg[negmem[ncount]][0],neg[negmem[ncount]][1],neg[negmem[ncount]][2]);
                negflag[negmem[ncount]]=3;

                phival=MIN(-((rand()%(irand))/drand)*epsi,-rmin*factor);

                lambda=1.0;
                qq=0;

                do
                {
                normal(a,neg[negmem[ncount]][0],neg[negmem[ncount]][1],neg[negmem[ncount]][2],neg[negmem[ncount]][3]);
                neg[negmem[ncount]][0] += lambda*(phival - neg[negmem[ncount]][3])*nvec[0];
                neg[negmem[ncount]][1] += lambda*(phival - neg[negmem[ncount]][3])*nvec[1];
                neg[negmem[ncount]][2] += lambda*(phival - neg[negmem[ncount]][3])*nvec[2];

                ii=int((neg[negmem[ncount]][0])/dx);
                jj=int((neg[negmem[ncount]][1])/dx);
                kk=int((neg[negmem[ncount]][2])/dx);
                check=boundcheck(p,a,ii,jj,kk,0);
                if(check==0)
                break;

                neg[negmem[ncount]][3] = phipol(p,a,neg[negmem[ncount]][0],neg[negmem[ncount]][1],neg[negmem[ncount]][2]);
                lambda/=2.0;
                ++qq;
                }while((neg[negmem[ncount]][3]<-epsi || neg[negmem[ncount]][3]>-rmin)&& qq<15);
				
				//negradius(p,a,negmem[ncount]);
				
                if((neg[negmem[ncount]][3]<-epsi || neg[negmem[ncount]][3]>-rmin)||check==0)
                {
                
                negflag[negmem[ncount]]=0;
				++ncount;
                --reseeded;
				success=0;
                }
            }

            if(ncount==0 && negactive<maxparticle)
            {	
                neg[negactive][0] = (double(i) + (rand()%(irand))/drand)*dx;
                neg[negactive][1] = (double(j) + (rand()%(irand))/drand)*dx;
                neg[negactive][2] = (double(k) + (rand()%(irand))/drand)*dx;
                neg[negactive][3] = phipol(p,a,neg[negactive][0],neg[negactive][1],neg[negactive][2]);
                negflag[negactive]=3;

                lambda=1.0;

                phival=MIN(-((rand()%(irand))/drand)*epsi,-rmin*factor);
                qq=0;

                do
                {

                normal(a,neg[negactive][0],neg[negactive][1],neg[negactive][2],neg[negactive][3]);
                neg[negactive][0] += lambda*(phival - neg[negactive][3])*nvec[0];
                neg[negactive][1] += lambda*(phival - neg[negactive][3])*nvec[1];
                neg[negactive][2] += lambda*(phival - neg[negactive][3])*nvec[2];

                ii=int((neg[negactive][0])/dx);
                jj=int((neg[negactive][1])/dx);
                kk=int((neg[negactive][2])/dx);
                check=boundcheck(p,a,ii,jj,kk,0);
                if(check==0)
                break;

                neg[negactive][3] = phipol(p,a,neg[negactive][0],neg[negactive][1],neg[negactive][2]);
                lambda/=2.0;
                ++qq;
                }while((neg[negactive][3]<-epsi || neg[negactive][3]>-rmin) && qq<15);


                if(neg[negactive][3]>=-epsi && neg[negactive][3]<=-rmin && check==1)
                {
				//negradius(p,a,negactive);
                ++negactive;
                ++reseeded;
				success=1;
                }
            }
	
	return success;
}
	
