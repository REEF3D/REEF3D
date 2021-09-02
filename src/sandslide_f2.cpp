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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sandslide_f2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sliceint.h"

sandslide_f2::sandslide_f2(lexer *p) : norm_vec(p), bedslope(p), fh(p)
{
    if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;

	fac1 = (1.0/6.0);
	fac2 = (1.0/12.0);
}

sandslide_f2::~sandslide_f2()
{
}

void sandslide_f2::start(lexer *p, fdm * a, ghostcell *pgc, sediment_fdm *s)
{
    // mainloop
    for(int qn=0; qn<p->S91; ++qn)
    {
        count=0;
        
        // fill
        SLICELOOP4
        fh(i,j)=0.0;
        
        pgc->gcsl_start4(p,fh,1);
        
        // slide loop
        SLICELOOP4
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
        {
            slide(p,a,pgc);
        }
        pgc->gcslparax_fh(p,fh,4);
        
        // fill back
        SLICELOOP4
        a->bedzh(i,j)+=fh(i,j);

        pgc->gcsl_start4(p,a->bedzh,1);

        count=pgc->globalimax(count);

        p->slidecells=count;
        
        if(p->slidecells==0)
        break;

        if(p->mpirank==0)
        cout<<"sandslide_f2 corrections: "<<p->slidecells<<endl;
    }
}

void sandslide_f2::slide(lexer *p, fdm * a, ghostcell *pgc)
{
		k = a->bedk(i,j);
		
		slope(p,a,pgc,teta,alpha,gamma,phi);

			
        // 1
        dh = a->bedzh(i,j) - a->bedzh(i-1,j);
        
        maxdh = tan(phi)*p->DXP[IM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
            dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DXP[IM1];
            
            fh(i,j)-= fac1*dh_corr;
            fh(i-1,j)+= fac1*dh_corr;
            
		++count;
		}

        // 2
        dh = a->bedzh(i,j) - a->bedzh(i+1,j);
        
        maxdh = tan(phi)*p->DXP[IP];
		
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
			dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DXP[IP];
            
            fh(i,j)-= fac1*dh_corr;
            fh(i+1,j)+= fac1*dh_corr;
			
        ++count;
        }

        // 3
        dh = a->bedzh(i,j) - a->bedzh(i,j-1);
        
        maxdh = tan(phi)*p->DYP[JM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{          
            dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DYP[JM1];
            
            fh(i,j)-= fac1*dh_corr;
            fh(i,j-1)+= fac1*dh_corr;
			
        ++count;
        }

        // 4
        dh = a->bedzh(i,j) - a->bedzh(i,j+1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DXP[IP];
        
        maxdh = tan(phi)*p->DYP[JP];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
            dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DXP[IP];
            
            fh(i,j)-= fac1*dh_corr;
            fh(i,j+1)+= fac1*dh_corr;

        ++count;
        }
		
		
        // 5
        dh = a->bedzh(i,j) - a->bedzh(i-1,j-1);
        
        maxdhs = tan(phi)*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
        {
            dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);
            
            fh(i,j)-= fac2*dh_corr;
            fh(i-1,j-1)+= fac2*dh_corr;
        
        ++count;
        }
    

        // 6
        dh = a->bedzh(i,j) - a->bedzh(i-1,j+1);

        maxdhs = tan(phi)*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{   
            dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);         
            
            fh(i,j)-= fac2*dh_corr;
            fh(i-1,j+1)+= fac2*dh_corr;
			
        ++count;
        }

        // 7
        dh = a->bedzh(i,j) - a->bedzh(i+1,j-1);
 
        maxdhs = tan(phi)*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{
            dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
            
			fh(i,j)-= fac2*dh_corr;
            fh(i+1,j-1)+= fac2*dh_corr;

        ++count;
        }
    
        // 8
        dh = a->bedzh(i,j) - a->bedzh(i+1,j+1);
  
        maxdhs = tan(phi)*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
		{   
            dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);   
      
            fh(i,j)-= fac2*dh_corr;
            fh(i+1,j+1)+= fac2*dh_corr;
            
        ++count;
        }
        
}

void sandslide_f2::topo_zh_update(lexer *p, fdm *a,ghostcell *pgc)
{
	pgc->gcsl_start4(p,a->bedzh,1);
	
    ALOOP
    {
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    a->topo(i,j,k)=-a->bedzh(i,j)+p->pos_z();
    }
	
	pgc->start4a(p,a->topo,150);
}

