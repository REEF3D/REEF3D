/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"print_runup.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void print_runup::cone(lexer *p, fdm* a, ghostcell *pgc)
{	
	int ds_num,q,check;
	int flagval1, flagval2;
	double xs,ys,zs,xe,ye,ze;
	double pos1_x,pos1_y,pos1_z,pos2_x,pos2_y,pos2_z;
	double length,ds,lsv1,lsv2;
	double step,fac,dnom;

	for(n=0;n<line_num;++n)
	{
	cut_x[n] = 0.0;
	cut_y[n] = 0.0;
	cut_z[n] = -100000.0;
	cut_active[n] = -1;
	}
	
	for(n=0;n<line_num*p->M10;++n)
	{
	cutall_x[n] = 0.0;
	cutall_y[n] = 0.0;
	cutall_z[n] = -100000.0;
	}
	
	
		cut_count=0;
		for(n=0;n<line_num;++n)
		{
			xs = line[n][0];
			ys = line[n][1];
			zs = line[n][2];
			
			xe = line[n][3];
			ye = line[n][4];
			ze = line[n][5];
			
				length = sqrt(pow(xe-xs,2.0) + pow(ye-ys,2.0) + pow(ze-zs,2.0));
				ds_num = int(length/(0.025*p->dx)) + 1;
				ds = length/double(ds_num);
				
			step=0.0;
			check=0;
			for(q=0;q<ds_num-1;++q)
			{
			pos1_x = xs + step*ds*(xe-xs);
			pos1_y = ys + step*ds*(ye-ys);
			pos1_z = zs + step*ds*(ze-zs);
			
			pos2_x = xs + (step+1.0)*ds*(xe-xs);
			pos2_y = ys + (step+1.0)*ds*(ye-ys);
			pos2_z = zs + (step+1.0)*ds*(ze-zs);
			
			
			
			

			
				if(positioncheck(p,a,pos1_x,pos1_y,pos1_z,0)>0 && positioncheck(p,a,pos2_x,pos2_y,pos2_z,0)>0)
				{
					i = int((pos1_x-p->originx)/p->dx);
					j = int((pos1_y-p->originy)/p->dx);
					k = int((pos1_z-p->originz)/p->dx);
					
					flagval1 = p->flag5[IJK];

					i = int((pos2_x-p->originx)/p->dx);
					j = int((pos2_y-p->originy)/p->dx);
					k = int((pos2_z-p->originz)/p->dx);
					
					flagval2 = p->flag5[IJK];
					
					//cout<<p->mpirank<<" ijk: "<<i<<" "<<j<<" "<<k<<"     "<<flagval2<<endl;
					
					//cout<<"flagval1: "<<flagval1<<"  "<<flagval2<<endl;
					
					if((flagval1>=0 || flagval1<=-10) && (flagval2>=0 || flagval2<=-10))
					{
					lsv1 = p->ccipol4phi(a,a->phi,pos1_x-p->originx,pos1_y-p->originy,pos1_z-p->originz);
					lsv2 = p->ccipol4phi(a,a->phi,pos2_x-p->originx,pos2_y-p->originy,pos2_z-p->originz);
					
					if((lsv1>=0.0 && lsv2<0.0) || (lsv1<0.0 && lsv2>=0.0))
					{
					dnom=lsv2-lsv1;
					dnom=fabs(dnom)>1.0e-20?dnom:1.0e-20;
					fac = lsv1/dnom;
					
						if(pos1_z + fac*ds*(ze-zs)>=zs && pos1_z + fac*ds*(ze-zs)<=ze && check==0)
						{
						cut_x[cut_count] = pos1_x + fac*ds*(xe-xs);
						cut_y[cut_count] = pos1_y + fac*ds*(ye-ys);
						cut_z[cut_count] = pos1_z + fac*ds*(ze-zs);
						cut_ID[cut_count] = n;
						
						++cut_count;
						check=1;
						}
					}
					}
				}
			step+=1.0;
			}
		}
		
	pgc->gather_int(&cut_count,1,cutall_count,1);
	
	for(q=0;q<p->M10;++q)
	displ[q]=0;
	
	if(p->mpirank==0)
	for(q=1;q<p->M10;++q)
	displ[q]=displ[q-1]+cutall_count[q-1];
	
	pgc->gatherv_int(cut_ID,cut_count,cutall_ID,cutall_count,displ);
	pgc->gatherv_double(cut_x,cut_count,cutall_x,cutall_count,displ);
	pgc->gatherv_double(cut_y,cut_count,cutall_y,cutall_count,displ);
	pgc->gatherv_double(cut_z,cut_count,cutall_z,cutall_count,displ);
	
	for(n=0;n<p->dgc4_count;++n)
    {
    i=p->dgc4[n][0];
    j=p->dgc4[n][1];
    k=p->dgc4[n][2];
	}
	
	
	if(p->mpirank==0)
	{
		cuttotal_count=0;
		for(q=0;q<p->M10;++q)
		cuttotal_count+=cutall_count[q];
		
		for(n=0;n<cuttotal_count;++n)
		{	
			
			if(cutall_z[n]>runup_z[cutall_ID[n]])
			{
			runup_x[cutall_ID[n]] = cutall_x[n];
			runup_y[cutall_ID[n]] = cutall_y[n];
			runup_z[cutall_ID[n]] = cutall_z[n];
			runup_active[cutall_ID[n]] = 1;
			}
		}

		for(n=0;n<cuttotal_count;++n)
		{
		cut_x[cutall_ID[n]] = cutall_x[n];
		cut_y[cutall_ID[n]] = cutall_y[n];
		cut_z[cutall_ID[n]] = cutall_z[n];
		cut_active[cutall_ID[n]] = 1;
		}
		
	}
}