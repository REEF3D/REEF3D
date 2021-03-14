/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"

void sflow_sediment_f::sandslide(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    fac1 = (1.0/6.0);
	fac2 = (1.0/12.0);

    if(p->S90==1)
    for(int qn=0; qn<p->S91; ++qn)
    {
        slidecount=0;
        
        SLICELOOP4
        fh(i,j)=b->bed(i,j);

        // slide loop
        SLICELOOP4
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
        {
            slide(p,b,pgc);
        }
        
        SLICELOOP4
        b->bed(i,j)=fh(i,j);

        pgc->gcsl_start4(p,b->bed,1);

        slidecount=pgc->globalimax(slidecount);

        p->slidecells=slidecount;
        
        if(p->slidecells==0)
        break;

        if(p->mpirank==0)
        cout<<"sandslide_f corrections: "<<p->slidecells<<endl;
    }
}

void sflow_sediment_f::slide(lexer *p, fdm2D *b, ghostcell *pgc)
{	
    double dh,dh_corr,maxdh,maxdhs;
    
        // 1
        dh = b->bed(i,j) - b->bed(i-1,j);
        dh_corr = dh + tan(p->S93*(PI/180.0))*p->DXP[IM1];
        
        maxdh = tan(phi(i,j))*p->DXP[IM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{   
            //dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*p->DXP[IM1];
			fh(i,j)-= fac1*dh_corr;
            fh(i-1,j)+= fac1*dh_corr;
            
		++slidecount;
		}

        // 2
        dh = b->bed(i,j) - b->bed(i+1,j);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DXP[IP];
        
        maxdh = tan(phi(i,j))*p->DXP[IP];
		
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
            //dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*p->DXP[IP];
            
			fh(i,j)-= fac1*dh_corr;
            fh(i+1,j)+= fac1*dh_corr;
			
        ++slidecount;
        }

        // 3
        dh = b->bed(i,j) - b->bed(i,j-1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DYP[JM1];
        
        maxdh = tan(phi(i,j))*p->DYP[JM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
            //dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*p->DYP[JM1];
            
			 fh(i,j)-= fac1*dh_corr;
            fh(i,j-1)+= fac1*dh_corr;
			
        ++slidecount;
        }

        // 4
        dh = b->bed(i,j) - b->bed(i,j+1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DXP[JP];
        
        maxdh = tan(phi(i,j))*p->DYP[JP];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
			//dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*p->DXP[JP];
            
            fh(i,j)-= fac1*dh_corr;
            fh(i,j+1)+= fac1*dh_corr;

        ++slidecount;
        }
		
		
        // 5
        dh = b->bed(i,j) - b->bed(i-1,j-1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);
        
        maxdhs = tan(phi(i,j))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
        {
			//dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);
            fh(i,j)-= fac2*dh_corr;
            fh(i-1,j-1)+= fac2*dh_corr;
        
        ++slidecount;
        }
    

        // 6
        dh = b->bed(i,j) - b->bed(i-1,j+1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);
        
        maxdhs = tan(phi(i,j))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{            
			//dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);
            fh(i,j)-= fac2*dh_corr;
            fh(i-1,j+1)+= fac2*dh_corr;
			
        ++slidecount;
        }

        // 7
        dh = b->bed(i,j) - b->bed(i+1,j-1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
        
        maxdhs = tan(phi(i,j))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{
			//dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
            fh(i,j)-= fac2*dh_corr;
            fh(i+1,j-1)+= fac2*dh_corr;

        ++slidecount;
        }
    
        // 8
        dh = b->bed(i,j) - b->bed(i+1,j+1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);
        
        maxdhs = tan(phi(i,j))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
		{            
			//dh_corr = dh-maxdh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);
            fh(i,j)-= fac2*dh_corr;
            fh(i+1,j+1)+= fac2*dh_corr;
            
        ++slidecount;
        }
    
}
