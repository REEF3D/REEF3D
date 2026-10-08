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
for more details->

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm.h"
#include"bedshear.h"
#include"sediment_fdm.h"
#include<vector>

void sediment_part::ini_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    double h;
    ILOOP
    JLOOP
    {
        h = p->ZN[0+marge];
        
        KLOOP
        PBASECHECK
        if(a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>0.0)
            h = -(a->topo(i,j,k-1)*p->DZP[KM1])/(a->topo(i,j,k)-a->topo(i,j,k-1)) + p->pos_z()-p->DZP[KM1];

        s->bedzh(i,j)=h;
        s->bedzh0(i,j)=h;
    }

    pgc->gcsl_start4(p,s->bedzh0,50);

    pst->seed_particles(p,a,pgc,s);

    fill_PQ_cfd(p,a,pgc);

    SLICELOOP4
        s->bedk(i,j)=0;

    SLICELOOP4
    {
        KLOOP
        PBASECHECK
        if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
        s->bedk(i,j)=k+1;

        s->reduce(i,j)=0.3;
    }

    pbedshear->taubed(p,a,pgc,s);
    pgc->gcsl_start4(p,s->tau_eff,1);
    pbedshear->taucritbed(p,a,pgc,s);
    pgc->gcsl_start4(p,s->tau_crit,1);

    pst->ini_fields(p,a,pgc,s);
    pst->print_particles(p,s);
    
    // S 10 2 two-way (Q 50 1): the bed is no solid for the fluid, the fluid flows through
    // the parcels and feels their drag; topo stays as the bed level for output
    if(p->S10==2 && p->Q50==1 && p->Q11==2)
    {
        p->toporead=0;
        p->topoforcing=0;
        
        topo_flags_off(p,a,pgc);
    }

    pgc->gcdf_update(p,a);
    
    pst->periodic_flags(p);
}

// the cells that are solid only because of the topo become fluid cells: flag4 and the face
// flags next to them are set again; faces on the domain boundary take the face flag of a
// fluid cell further up in the same column, so periodic and parallel sides keep their flags
void sediment_part::topo_flags_off(lexer *p, fdm *a, ghostcell *pgc)
{
    auto id = [&](int ii, int jj, int kk) {return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;};
    auto inside = [&](int ii, int jj, int kk) {return ii>=0 && ii<p->knox && jj>=0 && jj<p->knoy && kk>=0 && kk<p->knoz;};
    
    double psi = p->j_dir==1 ? -p->X41*(1.0/3.0)*(p->DXN[0+marge]+p->DYN[0+marge]+p->DZN[0+marge]) : -p->X41*0.5*(p->DXN[0+marge]+p->DZN[0+marge]);
    
    std::vector<int> ci,cj,ck;
    
    BASELOOP
    if(p->flag4[IJK]<0 && (p->solidread==0 || a->solid(i,j,k)>=psi))
    {
        p->flag4[IJK]=10;
        ci.push_back(i);
        cj.push_back(j);
        ck.push_back(k);
    }
    
    int nchg = pgc->globalisum(int(ci.size()));
    
    if(p->mpirank==0)
    cout<<"CPM two-way: "<<nchg<<" topo cells become fluid cells"<<endl;
    
    pgc->flagx(p,p->flag4);
    
    // face flag: flag(face q,q+e) = flag4(q), -10 if flag4(q)>0 and flag4(q+e)<0
    auto face = [&](int *fl, int ii, int jj, int kk, int di, int dj, int dk)
    {
        if(!inside(ii,jj,kk))
        return;
        
        if(inside(ii+di,jj+dj,kk+dk))
        {
            int f4 = p->flag4[id(ii,jj,kk)];
            fl[id(ii,jj,kk)] = (f4>0 && p->flag4[id(ii+di,jj+dj,kk+dk)]<0) ? -10 : f4;
        }
        else
        {
            // boundary face: copy from a fluid cell of the same column (x,y faces) or row (z faces)
            if(dk==0)
            {
                for(int q=p->knoz-1;q>=0;--q)
                if(p->flag4[id(ii,jj,q)]>0 && q!=kk)
                {
                    fl[id(ii,jj,kk)] = fl[id(ii,jj,q)];
                    break;
                }
            }
            else
            {
                for(int q=0;q<p->knox;++q)
                if(p->flag4[id(q,jj,kk)]>0 && q!=ii)
                {
                    fl[id(ii,jj,kk)] = fl[id(q,jj,kk)];
                    break;
                }
            }
        }
    };
    
    for(size_t q=0;q<ci.size();++q)
    {
        i=ci[q]; j=cj[q]; k=ck[q];
        
        face(p->flag1,i,j,k,1,0,0);
        face(p->flag1,i-1,j,k,1,0,0);
        
        if(p->j_dir==1)
        {
        face(p->flag2,i,j,k,0,1,0);
        face(p->flag2,i,j-1,k,0,1,0);
        }
        else
        p->flag2[IJK] = p->flag4[IJK];
        
        face(p->flag3,i,j,k,0,0,1);
        face(p->flag3,i,j,k-1,0,0,1);
    }
    
    pgc->flagx(p,p->flag1);
    pgc->flagx(p,p->flag2);
    pgc->flagx(p,p->flag3);
}
