/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"sandslide_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sliceint.h"

sandslide_f::sandslide_f(lexer *p) : norm_vec(p), bedslope(p), fh(p)
{
    if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;

	fac1=p->S92*(1.0/6.0);
	fac2=p->S92*(1.0/12.0);
}

sandslide_f::~sandslide_f()
{
}

void sandslide_f::start(lexer *p, fdm * a, ghostcell *pgc)
{
    ALOOP
    {
    slope(p,a,pgc,teta,alpha,gamma,phi);
    a->test(i,j,k)=phi*180.0/PI;
    }
    

    for(int qn=0; qn<p->S91; ++qn)
    {
        
    SLICELOOP4
    fh(i,j)=a->bedzh(i,j);
    
    count=0;
    pgc->dgcslpol(p,a->bedzh,p->dgcsl4,p->dgcsl4_count,14);
    a->bedzh.ggcpol(p);
    
    topo_zh_update(p,a,pgc);

    SLICELOOP4
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    {
		slide(p,a,pgc);
    }
    
    SLICELOOP4
    a->bedzh(i,j)=fh(i,j);

    pgc->gcsl_start4(p,a->bedzh,1);

    count=pgc->globalimax(count);

    p->slidecells=count;
    
    if(p->slidecells==0)
    break;

    if(p->mpirank==0)
    cout<<"sandslide_f corrections: "<<p->slidecells<<endl;
    }
}

void sandslide_f::slide(lexer *p, fdm * a, ghostcell *pgc)
{
		k = a->bedk(i,j);
		
		slope(p,a,pgc,teta,alpha,gamma,phi);

			
        // 1
        dh = a->bedzh(i,j) - a->bedzh(i-1,j);
        dh_corr = dh + tan(p->S93*(PI/180.0))*p->DXP[IM1];
        
        maxdh = tan(phi)*p->DXP[IM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
			fh(i,j)-= fac1*dh_corr;
            fh(i-1,j)+= fac1*dh_corr;
            
		++count;
		}

        // 2
        dh = a->bedzh(i,j) - a->bedzh(i+1,j);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DXP[IP];
        
        maxdh = tan(phi)*p->DXP[IP];
		
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
			fh(i,j)-= fac1*dh_corr;
            fh(i+1,j)+= fac1*dh_corr;
			
        ++count;
        }

        // 3
        dh = a->bedzh(i,j) - a->bedzh(i,j-1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DYP[JM1];
        
        maxdh = tan(phi)*p->DYP[JM1];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{            
			fh(i,j)-= fac1*dh_corr;
            fh(i,j-1)+= fac1*dh_corr;
			
        ++count;
        }

        // 4
        dh = a->bedzh(i,j) - a->bedzh(i,j+1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*p->DXP[JP];
        
        maxdh = tan(phi)*p->DYP[JP];
        
        if(dh>maxdh && fabs(dh)<1.0e15)
		{
			fh(i,j)-= fac1*dh_corr;
            fh(i,j+1)+= fac1*dh_corr;

        ++count;
        }
		
		
        // 5
        dh = a->bedzh(i,j) - a->bedzh(i-1,j-1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);
        
        maxdhs = tan(phi)*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JM1]*p->DYP[JM1]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
        {
			fh(i,j)-= fac2*dh_corr;
            fh(i-1,j-1)+= fac2*dh_corr;
        
        ++count;
        }
    

        // 6
        dh = a->bedzh(i,j) - a->bedzh(i-1,j+1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);
        
        maxdhs = tan(phi)*sqrt(p->DXP[IM1]*p->DXP[IM1] + p->DYP[JP]*p->DYP[JP]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{            
			fh(i,j)-= fac2*dh_corr;
            fh(i-1,j+1)+= fac2*dh_corr;
			
        ++count;
        }

        // 7
        dh = a->bedzh(i,j) - a->bedzh(i+1,j-1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
        
        maxdhs = tan(phi)*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JM1]*p->DYP[JM1]);
        
        if(dh>maxdhs && fabs(dh)<1.0e15)
		{
			fh(i,j)-= fac2*dh_corr;
            fh(i+1,j-1)+= fac2*dh_corr;

        ++count;
        }
    
        // 8
        dh = a->bedzh(i,j) - a->bedzh(i+1,j+1);
		dh_corr = dh + tan(p->S93*(PI/180.0))*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);
        
        maxdhs = tan(phi)*sqrt(p->DXP[IP]*p->DXP[IP] + p->DYP[JP]*p->DYP[JP]);

        if(dh>maxdhs && fabs(dh)<1.0e15)
		{            
			fh(i,j)-= fac2*dh_corr;
            fh(i+1,j+1)+= fac2*dh_corr;
            
        ++count;
        }
        
}

void sandslide_f::topo_zh_update(lexer *p, fdm *a,ghostcell *pgc)
{
	pgc->gcsl_start4(p,a->bedzh,1);
	
    ALOOP
    {
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    a->topo(i,j,k)=-a->bedzh(i,j)+p->pos_z();
    }
	
	pgc->start4a(p,a->topo,150);
}

