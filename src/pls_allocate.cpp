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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_pls::allocate(lexer* p,fdm* a,ghostcell* pgc)
{
     maxparticle = int(p->F33*double(p->cellnum*pnum/2));
	 
	 p->Darray(neg,maxparticle,5);
	 p->Darray(pos,maxparticle,5);
	 
	 p->Iarray(negflag,maxparticle);
	 p->Iarray(posflag,maxparticle);
	 
	 p->Iarray(negmem,maxparticle);
	 p->Iarray(posmem,maxparticle);

     for(n=0;n<maxparticle;++n)
     {

         negflag[n] = 0;
         posflag[n] = 0;

         negmem[n] = 0;
         posmem[n] = 0;
     }

    // parallel
	p->Iarray(pxs,6);
	p->Iarray(nxs,6);
	
	p->Iarray(pxr,6);
	p->Iarray(nxr,6);


    posxs= new double*[6];
    negxs= new double*[6];

    for(n=0;n<6;++n)
	{	
	if(p->gcpara1_count>0)
	posxs[0] = new double[5*p->gcpara1_count*pnum];
	if(p->gcpara2_count>0)
    posxs[1] = new double[5*p->gcpara2_count*pnum];
	if(p->gcpara3_count>0)
	posxs[2] = new double[5*p->gcpara3_count*pnum];
	if(p->gcpara4_count>0)
    posxs[3] = new double[5*p->gcpara4_count*pnum];
	if(p->gcpara5_count>0)
    posxs[4] = new double[5*p->gcpara5_count*pnum];
	if(p->gcpara6_count>0)
    posxs[5] = new double[5*p->gcpara6_count*pnum];
	
	if(p->gcpara1_count>0)
    negxs[0] = new double[5*p->gcpara1_count*pnum];
	if(p->gcpara2_count>0)
    negxs[1] = new double[5*p->gcpara2_count*pnum];
	if(p->gcpara3_count>0)
    negxs[2] = new double[5*p->gcpara3_count*pnum];
	if(p->gcpara4_count>0)
    negxs[3] = new double[5*p->gcpara4_count*pnum];
	if(p->gcpara5_count>0)
    negxs[4] = new double[5*p->gcpara5_count*pnum];
	if(p->gcpara6_count>0)
    negxs[5] = new double[5*p->gcpara6_count*pnum];
    }

    posxr= new double*[6];
    negxr= new double*[6];

    for(n=0;n<6;++n)
    {
	if(p->gcpara1_count>0)
	posxr[0] = new double[5*p->gcpara1_count*pnum];
	if(p->gcpara2_count>0)
    posxr[1] = new double[5*p->gcpara2_count*pnum];
	if(p->gcpara3_count>0)
    posxr[2] = new double[5*p->gcpara3_count*pnum];
	if(p->gcpara4_count>0)
	posxr[3] = new double[5*p->gcpara4_count*pnum];
	if(p->gcpara5_count>0)
	posxr[4] = new double[5*p->gcpara5_count*pnum];
	if(p->gcpara6_count>0)
	posxr[5] = new double[5*p->gcpara6_count*pnum];
	
	if(p->gcpara1_count>0)
    negxr[0] = new double[5*p->gcpara1_count*pnum];
	if(p->gcpara2_count>0)
    negxr[1] = new double[5*p->gcpara2_count*pnum];
	if(p->gcpara3_count>0)
    negxr[2] = new double[5*p->gcpara3_count*pnum];
	if(p->gcpara4_count>0)
    negxr[3] = new double[5*p->gcpara4_count*pnum];
	if(p->gcpara5_count>0)
    negxr[4] = new double[5*p->gcpara5_count*pnum];
	if(p->gcpara6_count>0)
    negxr[5] = new double[5*p->gcpara6_count*pnum];
    }

}
