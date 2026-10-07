/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
#include"sediment_mixture.h"
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
    
    SLICEBASELOOP
    {
    s->slide_fh(i,j)=0.0;
    ci(i,j)=0.0;
    }
    
    // smallest cell for the pseudo time step
    dxmin=1.0e20;
    SLICELOOP4
    {
    dxmin = MIN(dxmin,p->DXN[IP]);
    
    if(p->j_dir==1 && p->gknoy>1)
    dxmin = MIN(dxmin,p->DYN[JP]);
    }
    dxmin = pgc->globalmin(dxmin);
    
    // mainloop
    for(int qn=0; qn<p->S91; ++qn)
    {
        count=0;
        
        // fill
        SLICEBASELOOP
        {
        fh(i,j)=0.0;
        
        diff_update(p,pgc,s);
        }
        
        pgc->gcsl_start4(p,fh,1);
        
        if(s->pmix!=nullptr)
        s->pmix->slide_zero(p,pgc);
        pgc->gcsl_start4(p,ci,1);
        

        
        // slide loop: sediment cells in the S 77 window only
        SEDSLICELOOP
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
        {
            slide(p,pgc,s);
        }
        
        pgc->gcslparax_fh(p,fh,4);
        
        // fill back
        SEDSLICELOOP
        {
        s->slide_fh(i,j)+=fh(i,j);
        s->bedzh(i,j)+=fh(i,j);
        }
        
        pgc->gcsl_start4(p,s->bedzh,1);
        
        // multi-fraction bed: sorting of the slid material
        if(s->pmix!=nullptr)
        s->pmix->slide_finish(p,pgc,s);

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
    // finite-volume diffusion of the bed where it is steeper than the angle of repose:
    // face coefficients dt/(cell width * centre distance), so the face flux is the same seen from
    // both cells (conservative on non-uniform grids); faces to cells without an erodible bed
    // (structures, DFBED<0, solids, domain boundary) and to cells outside the S 77 window are closed.
    // pseudo time step from the smallest cell (dxmin, start())
    double dt = 0.1*dxmin*dxmin;
    double fcf[4];
    const int ni[4] = {1,-1,0,0};
    const int nj[4] = {0,0,1,-1};
    
    fcf[0] = dt/(p->DXN[IP]*p->DXP[IP]);
    fcf[1] = dt/(p->DXN[IP]*p->DXP[IM1]);
    fcf[2] = dt/(p->DYN[JP]*p->DYP[JP]);
    fcf[3] = dt/(p->DYN[JP]*p->DYP[JM1]);
    
    for(int f=0;f<4;++f)
    {
    int ii=i+ni[f];
    int jj=j+nj[f];
    
        if(!SLIDE_NB(ni[f],nj[f]) || p->flagslice4[(ii-p->imin)*p->jmax + jj-p->jmin]<0 || p->DFBED[(ii-p->imin)*p->jmax + jj-p->jmin]<0
          || p->XP[IP+ni[f]]<=p->S77_xs || p->XP[IP+ni[f]]>=p->S77_xe)
        fcf[f] = 0.0;
    }

    fh(i,j) = 0.0;
    
    for(int f=0;f<4;++f)
    fh(i,j) += fcf[f]*(s->bedzh(i+ni[f],j+nj[f])-s->bedzh(i,j))*0.5*(ci(i+ni[f],j+nj[f])+ci(i,j));
    
    // multi-fraction bed: face fluxes with upwind composition
    if(s->pmix!=nullptr)
    s->pmix->slide_pde(p,s,ci,i,j,fcf);
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


