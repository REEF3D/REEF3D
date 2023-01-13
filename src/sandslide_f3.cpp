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

#include"sandslide_f3.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sliceint.h"

sandslide_f3::sandslide_f3(lexer *p) : norm_vec(p), bedslope(p), fh(p)
{
    if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;

	fac1 = p->S92*(1.0/6.0);
	fac2 = p->S92*(1.0/12.0);
}

sandslide_f3::~sandslide_f3()
{
}

void sandslide_f3::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    SLICELOOP4
    s->slideflag(i,j)=0.0;
    
    if(p->mpirank==0)
    cout<<"sandslide_f3"<<endl;
    
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
            slide(p,pgc,s);
        }
        
        pgc->gcslparax_fh(p,fh,4);

        pgc->gcsl_start4(p,s->bedzh,1);

        count=pgc->globalimax(count);

        p->slidecells=count;
        
        if(p->slidecells==0)
        break;

        if(p->mpirank==0)
        cout<<"sandslide_f3 corrections: "<<p->slidecells<<endl;
    }
}

void sandslide_f3::slide(lexer *p, ghostcell *pgc, sediment_fdm *s)
{   
        
        double dzp=0.0;
        int Iup=0;
        int id[8];
		k = s->bedk(i,j);
        
        
        for(int qn=0; qn<8; ++qn)
        id[qn]=0;
        
			
        // 1
        dh = s->bedzh(i-1,j) - s->bedzh(i,j);
        
        maxdh = tan(s->phi(i,j))*p->DXP[IM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
            
            dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DXP[IM1];

            dzp += dh_corr;
            ++Iup;
            
            id[0]=1;
            
		++count;
		}

        // 2
        dh = s->bedzh(i+1,j) - s->bedzh(i,j);
        
        maxdh = tan(s->phi(i,j))*p->DXP[IP];
		
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
			dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DXP[IP];

            dzp += dh_corr;
            ++Iup;
            
            id[1]=1;
            
        ++count;
        }

        // 3
        dh = s->bedzh(i,j-1) - s->bedzh(i,j);
        
        maxdh = tan(s->phi(i,j))*p->DYP[JM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{          
            dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DYP[JM1];

            dzp += dh_corr;
            ++Iup;
            
            id[2]=1;
            
        ++count;
        }

        // 4
        dh = s->bedzh(i,j+1) - s->bedzh(i,j);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DYP[JP];
        
        maxdh = tan(s->phi(i,j))*p->DYP[JP];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
            dh_corr = (dh-maxdh) + tan(p->S93*(PI/180.0))*p->DYP[JP];

            dzp += dh_corr;
            ++Iup;
            
            id[3]=1;

        ++count;
        }
		
		
        // 5
        dh = s->bedzh(i-1,j-1) - s->bedzh(i,j);
        
        maxdhs = tan(s->phi(i,j))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
        {
            dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);

            dzp += dh_corr;
            ++Iup;
            
            id[4]=1;
              
        ++count;
        }
    

        // 6
        dh = s->bedzh(i-1,j+1) - s->bedzh(i,j);

        maxdhs = tan(s->phi(i,j))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{   
             dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);  

            dzp += dh_corr;
            ++Iup;
            
            id[5]=1;
            	
        ++count;
        }

        // 7
        dh = s->bedzh(i+1,j-1) - s->bedzh(i,j);
 
        maxdhs = tan(s->phi(i,j))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{
            dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);

            dzp += dh_corr;
            ++Iup;
            
            id[6]=1;
            
        ++count;
        }
    
    
        // 8
        dh =  s->bedzh(i+1,j+1) - s->bedzh(i,j);
  
        maxdhs = tan(s->phi(i,j))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
		{   
            dh_corr = (dh-maxdhs) + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);   

            dzp += dh_corr;
            ++Iup;
            
            id[7]=1;
                      
        ++count;
        }
        
        dzp = dzp/double(Iup + 1);
        
    
        
        if(id[0]==1)
        s->bedzh(i-1,j) += dzp + s->bedzh(i,j) - s->bedzh(i-1,j) + tan(s->phi(i,j))*p->DXP[IM1] - tan(p->S93*(PI/180.0))*p->DXP[IM1]; 
        
        if(id[1]==1)
        s->bedzh(i+1,j) += dzp + s->bedzh(i,j) - s->bedzh(i+1,j) + tan(s->phi(i,j))*p->DXP[IP] - tan(p->S93*(PI/180.0))*p->DXP[IP];
        
        if(id[2]==1)
        s->bedzh(i,j-1) += dzp + s->bedzh(i,j) - s->bedzh(i,j-1) + tan(s->phi(i,j))*p->DYP[JM1] - tan(p->S93*(PI/180.0))*p->DYP[JM1];
        
        if(id[3]==1)
        s->bedzh(i,j+1) += dzp + s->bedzh(i,j) - s->bedzh(i,j+1) + tan(s->phi(i,j))*p->DYP[JP] - tan(p->S93*(PI/180.0))*p->DYP[JP];
        
        if(id[4]==1)
        s->bedzh(i-1,j-1) += dzp + s->bedzh(i,j) - s->bedzh(i-1,j-1) + tan(s->phi(i,j))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]) - tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);
        
        if(id[5]==1)
        s->bedzh(i-1,j+1) += dzp + s->bedzh(i,j) - s->bedzh(i-1,j+1) + tan(s->phi(i,j))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]) - tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);  
        
        if(id[6]==1)
        s->bedzh(i+1,j-1) += dzp + s->bedzh(i,j) - s->bedzh(i+1,j-1) + tan(s->phi(i,j))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]) - tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
        
        if(id[7]==1)
        s->bedzh(i+1,j+1) += dzp + s->bedzh(i,j) - s->bedzh(i+1,j+1) + tan(s->phi(i,j))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]) -tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);
        
        s->bedzh(i,j) += dzp;
        s->slideflag(i,j) += dzp;
        
}
