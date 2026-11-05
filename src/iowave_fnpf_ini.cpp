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

#include"iowave.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void iowave::ini_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    wavegen_precalc_ini(p,pgc);
    wavegen_precalc_relax_func_fnpf(p,pgc);
    
    // generation zone check: no cell inside it, or all of them dry, means no waves at all
    if(relax_gen && p->B90>0)
    {
        int ngen=0, ngenwet=0;
        
        SLICELOOP4
        if(wgflag(i,j)==1)
        {
        ++ngen;
        
        if(p->wet[IJ]>0)
        ++ngenwet;
        }
        
        ngen=pgc->globalisum(ngen);
        ngenwet=pgc->globalisum(ngenwet);
        
        if(p->mpirank==0)
        {
            if(ngen==0)
            cout<<endl<<"!!! WARNING iowave: the wave generation zone (B 108 / B 96) has no cell inside the domain, no waves are generated !!!"<<endl<<endl;
            
            else if(ngenwet==0)
            cout<<endl<<"!!! WARNING iowave: all "<<ngen<<" cells of the wave generation zone are dry, no waves are generated !!!"<<endl<<endl;
            
            else if(ngenwet<ngen)
            cout<<"iowave: wave generation zone "<<ngen<<" cells, "<<ngenwet<<" wet"<<endl;
        }
    }
    
    if(p->B89==1 && p->B98==2)
    wavegen_precalc_decomp_space_fnpf(p,pgc);
    
    if(p->B89==1 && p->B98==3)
    wavegen_precalc_decomp_space_dirichlet_fnpf(p,pgc);

    wavegen_precalc_fnpf(p,c,pgc);
    
    if(p->I30==0)
    {
	fivec_relax(p,pgc,c->Fi);
    fifsf_relax(p,pgc,c->Fifsf);
    eta_relax(p,pgc,c->eta);
    }  
    
    if(p->I30==1)
	full_initialize_fnpf(p,c,pgc);
    
    if(p->I30==2)
	full_initialize_fnpf_masked(p,c,pgc);
    
    SLICELOOP4
    c->WL(i,j) = MAX(c->wd_criterion, c->eta(i,j) + c->depth(i,j));
}

