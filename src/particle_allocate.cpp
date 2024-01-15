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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_f::allocate(lexer* p,fdm* a,ghostcell* pgc)
{
	size_t maxparticle = ceil(p->Q25*double(gpartnum));
	
	PP.reserve(maxparticle);

    // parallel
	p->Iarray(pxs,6);
	p->Iarray(pxr,6);


    posxs= new double*[6];
    int fac=2;

    for(n=0;n<6;++n)
	{	
		if(p->gcpara1_count>0)
			posxs[0] = new double[PP.entries*p->gcpara1_count*ppcell*fac];
		if(p->gcpara2_count>0)
			posxs[1] = new double[PP.entries*p->gcpara2_count*ppcell*fac];
		if(p->gcpara3_count>0)
			posxs[2] = new double[PP.entries*p->gcpara3_count*ppcell*fac];
		if(p->gcpara4_count>0)
			posxs[3] = new double[PP.entries*p->gcpara4_count*ppcell*fac];
		if(p->gcpara5_count>0)
			posxs[4] = new double[PP.entries*p->gcpara5_count*ppcell*fac];
		if(p->gcpara6_count>0)
			posxs[5] = new double[PP.entries*p->gcpara6_count*ppcell*fac];
    }

    posxr= new double*[6];

    for(n=0;n<6;++n)
    {
		if(p->gcpara1_count>0)
			posxr[0] = new double[PP.entries*p->gcpara1_count*ppcell*fac];
		if(p->gcpara2_count>0)
			posxr[1] = new double[PP.entries*p->gcpara2_count*ppcell*fac];
		if(p->gcpara3_count>0)
			posxr[2] = new double[PP.entries*p->gcpara3_count*ppcell*fac];
		if(p->gcpara4_count>0)
			posxr[3] = new double[PP.entries*p->gcpara4_count*ppcell*fac];
		if(p->gcpara5_count>0)
			posxr[4] = new double[PP.entries*p->gcpara5_count*ppcell*fac];
		if(p->gcpara6_count>0)
			posxr[5] = new double[PP.entries*p->gcpara6_count*ppcell*fac];
    }

}
