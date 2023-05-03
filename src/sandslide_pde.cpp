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
#include"sandslide_pde.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"ghostcell.h"

sandslide_pde::sandslide_pde(lexer *p) : norm_vec(p), bedslope(p), fh(p), ci(p)
{
    if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;

	dxs=sqrt(2.0*p->DXM*p->DXM);
	fac1 = (1.0/6.0);
	fac2 = (1.0/12.0);
}

sandslide_pde::~sandslide_pde()
{
}

void sandslide_pde::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    
    SLICELOOP4
    {
    s->slideflag(i,j)=0.0;
    ci(i,j)=0.0;
    }
    
    // mainloop
    for(int qn=0; qn<p->S91; ++qn)
    {
        count=0;
        
        // fill
        SLICELOOP4
        {
        fh(i,j)=0.0;
        
        diff_update(p,pgc,s);
        }
        
        pgc->gcsl_start4(p,fh,1);
        pgc->gcsl_start4(p,ci,1);
        

        
        // slide loop
        SLICELOOP4
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
        {
            slide(p,pgc,s);
        }
        
        pgc->gcslparax_fh(p,fh,4);
        
        // fill back
        SLICELOOP4
        {
        s->slideflag(i,j)+=fh(i,j);
        s->bedzh(i,j)+=fh(i,j);
        }
        
        pgc->gcsl_start4(p,s->bedzh,1);

        count=pgc->globalimax(count);

        p->slidecells=count;
        
        if(p->slidecells==0)
        break;

        if(p->mpirank==0)
        cout<<"sandslide_ped corrections: "<<p->slidecells<<endl;
    }
}

void sandslide_pde::slide(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double dt = 0.1*p->DXM*p->DXM;
    double sqd = (1.0/(p->DXM*p->DXM));

    fh(i,j) =  dt*sqd*( (s->bedzh(i+1,j)-s->bedzh(i,j))*0.5*(ci(i+1,j)+ci(i,j)) 
                        -(s->bedzh(i,j)-s->bedzh(i-1,j))*0.5*(ci(i,j)+ci(i-1,j))
                                
                        +(s->bedzh(i,j+1)-s->bedzh(i,j))*0.5*(ci(i,j+1)+ci(i,j)) 
                        -(s->bedzh(i,j)-s->bedzh(i,j-1))*0.5*(ci(i,j)+ci(i,j-1)));
    
  
}

void sandslide_pde::diff_update(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double uvel,vvel;
    double nx,ny,nz,norm;
    double nx0,ny0;
    double nz0,bx0,by0,gamma;
    
    int kmem=0;
    double dH;
    
    k = s->bedk(i,j);
        

    dH = sqrt(pow((s->bedzh(i+1,j)-s->bedzh(i-1,j))/p->DXM,2.0) + pow((s->bedzh(i,j+1)-s->bedzh(i,j-1))/p->DXM,2.0));
        
        
    bx0 = (s->bedzh(i+1,j)-s->bedzh(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
    by0 = (s->bedzh(i,j+1)-s->bedzh(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
     
    gamma = atan(sqrt(bx0*bx0 + by0*by0));


            if(gamma>s->phi(i,j))
            {
            ci(i,j) = 1.0;
            
            ++count;
            }
            
            if(gamma<s->phi(i,j))
            ci(i,j) = 0.0;

}


