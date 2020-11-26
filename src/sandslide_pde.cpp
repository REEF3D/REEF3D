/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"sandslide_pde.h"
#include"lexer.h"
#include"fdm.h"
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

void sandslide_pde::start(lexer *p, fdm * a, ghostcell *pgc)
{
    /*for(int qn=0; qn<p->S91; ++qn)
    {
    count=0;
    
    ALOOP
    ci(i,j,k)=0.0;
    pgc->start4a(p,ci,150);
    
    
    topo_zh_update(p,a,pgc,zh);
    
    ILOOP
    JLOOP
    {
		diff_update(p,a,pgc,zh);
    }
    
    pgc->start4a(p,ci,150);

    ILOOP
    JLOOP
    {
		slide(p,a,pgc,zh);
    }
    
    ALOOP 
    zh(i,j,k)=fh(i,j,k);

    pgc->start4a(p,zh,150);

    count=pgc->globalimax(count);

    p->slidecells=count;
    
    if(p->slidecells==0)
    break;

    if(p->mpirank==0)
    cout<<"sandslide_pde corrections: "<<p->slidecells<<endl;
    }*/
}

void sandslide_pde::slide(lexer *p, fdm * a, ghostcell *pgc)
{
    double dt = 0.1*p->DXM*p->DXM;
    double sqd = (1.0/(p->DXM*p->DXM));

    /*
    KLOOP
    PBASECHECK
    {
    fh(i,j,k) = zh(i,j,k) + dt*sqd*( (zh(i+1,j,k)-zh(i,j,k))*0.5*(ci(i+1,j,k)+ci(i,j,k)) -(zh(i,j,k)-zh(i-1,j,k))*0.5*(ci(i,j,k)+ci(i-1,j,k))
                                    +(zh(i,j+1,k)-zh(i,j,k))*0.5*(ci(i,j+1,k)+ci(i,j,k)) -(zh(i,j,k)-zh(i,j-1,k))*0.5*(ci(i,j,k)+ci(i,j-1,k)));
    }
    */
  
}

void sandslide_pde::diff_update(lexer *p, fdm * a, ghostcell *pgc)
{/*
    int kmem=0;
    double dH;

        KLOOP
        PBASECHECK
        {
            if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
            kmem=k+1;
        }

		k = kmem;
		
		slope(p,a,pgc,zh,teta,alpha,gamma,phi);
        

        dH = sqrt(pow((zh(i+1,j,k)-zh(i-1,j,k))/p->DXM,2.0) + pow((zh(i,j+1,k)-zh(i,j-1,k))/p->DXM,2.0));

            if(fabs(dH)>=tan(phi))
            {
            //cout<<i<<" "<<j<<"  SAND "<<gamma*(180.0/PI)<<" "<<teta*(180.0/PI)<<"  "<<phi*(180.0/PI)<<endl;
            KLOOP
            PBASECHECK
            ci(i,j,k) = 1.0;
            
            ++count;
            }
            
            if(fabs(dH)<tan(phi))
            KLOOP
            PBASECHECK
            ci(i,j,k) = 0.0;*/

}

void sandslide_pde::topo_zh_update(lexer *p, fdm *a,ghostcell *pgc)
{/*
	pgc->start4a(p,zh,150);
	
    ALOOP
    {
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    a->topo(i,j,k)=-zh(i,j,k)+p->pos_z();
    }
	
	pgc->start4a(p,a->topo,150);*/
}


