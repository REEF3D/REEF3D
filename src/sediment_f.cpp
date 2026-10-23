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

#include"sediment_f.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"topo.h"
#include"bedshear.h"
#include"patchBC_interface.h"
#include"sediment_exner.h"
#include"sediment_mixture.h"
#include"nhflow_suspended.h"
#include"ckpt_pack.h"

sediment_f::sediment_f(lexer *p, ghostcell *pgc, turbulence *pturb, patchBC_interface *ppBC)
{

    pBC = ppBC;
    
    sediment_logic(p,pgc,pturb);

	
    
    volume_token=0;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
    psed = this;
}

sediment_f::sediment_f(lexer *p, ghostcell *pgc, turbulence *pturb, patchBC_interface *ppBC, int amr) : amrgrid(amr)
{
    pBC = ppBC;

    sediment_logic(p,pgc,pturb);

    volume_token=0;

    if(p->F50==1)
	gcval_eta = 51;

    if(p->F50==2)
	gcval_eta = 52;

    if(p->F50==3)
	gcval_eta = 53;

    if(p->F50==4)
	gcval_eta = 54;

    psed = this;
}

sediment_f::~sediment_f()
{
    // the grids of the mesh refinement come and go with the regrids (the level-0 object lives for the run)
    delete pbed;
    delete pcbed;
    delete pslide;
    delete prelax;
    delete pvrans;
    delete pslope;
    delete preduce;
    delete ptopo;
    delete pnhfsusp;
    delete pnhfsuspdisc;
    delete pnhfsuspdiff;
    delete pbedshear;
    delete pbeddir;
    delete pks;
    delete s;
}

void sediment_f::start_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo *preto, solver *psolv)
{
    sedcalc=0;
    
	if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psolv);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psolv);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime)
		{            
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psolv);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
    
    sedcalc=1;
	}
    
    if(sedcalc==0)
    {
    fill_bedk(p,a,pgc);
    waterlevel(p,a,pgc);
    pbedshear->taubed(p,a,pgc,s);
    }
}

void sediment_f::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
    if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47 && p->count>0))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm_nhflow(p,d,pgc,pflow);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm_nhflow(p,d,pgc,pflow);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime )
		{
		sediment_algorithm_nhflow(p,d,pgc,pflow);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
	}
}

void sediment_f::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &P, slice &Q)
{
    if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47 && p->count>0))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime )
		{
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
	}
}





// ---------------------------------------------------------------------------------------------
// checkpoints (P 48, I 45): the scalars of the bed update (sediment_exner), of the wave-averaged bed shear
// (S 16 9), of the multi-fraction bed (S 51) and of the NHFLOW suspended load; the slices and fields of
// all sediment objects are in the registry of the grid (CFD: cfd_ckpt; NHFLOW, SFLOW: under the tag "reg")
// ---------------------------------------------------------------------------------------------

void sediment_f::ckpt_state(lexer *p, std::vector<double> &v)
{
    ckpt_put w{v};
    w.d(1.0);   // version
    w.i(volume_token); w.d(volume0);
    std::vector<double> ex;
    exner()->ckpt_state(ex);
    w.vd(ex);
    w.d(*pbedshear->ckpt_tw_old()); w.d(*pbedshear->ckpt_tw_dt()); w.i(*pbedshear->ckpt_tw_count());
    const int nf = (s->pmix!=nullptr) ? s->pmix->nf : 0;
    w.i(nf);
    for(int q=0; q<nf; ++q)
    w.i(s->pmix->noneq_ini_k[q]);
    int *wl = (pnhfsusp!=nullptr) ? pnhfsusp->ckpt_wlini() : nullptr;
    w.i(wl!=nullptr ? *wl : 0);
}

bool sediment_f::ckpt_restore(lexer *p, ghostcell *pgc, const std::vector<double> &v)
{
    ckpt_get r{v};
    double version=0.0;
    r.d(version);
    if(!r.ok || int(version)!=1)
    return false;
    r.i(volume_token); r.d(volume0);
    std::vector<double> ex;
    r.vd(ex);
    size_t k=0;
    exner()->ckpt_restore(ex,k);
    r.d(*pbedshear->ckpt_tw_old()); r.d(*pbedshear->ckpt_tw_dt()); r.i(*pbedshear->ckpt_tw_count());
    int nf=0;
    r.i(nf);
    if(nf!=((s->pmix!=nullptr) ? s->pmix->nf : 0))
    return false;
    for(int q=0; q<nf; ++q)
    r.i(s->pmix->noneq_ini_k[q]);
    int wl=0;
    r.i(wl);
    if(pnhfsusp!=nullptr && pnhfsusp->ckpt_wlini()!=nullptr)
    *pnhfsusp->ckpt_wlini() = wl;
    return r.ok && k==ex.size();
}

// the rank: the cumulative limiter volume of the multi-fraction bed (its log)
void sediment_f::ckpt_rank_state(lexer *p, std::vector<double> &v)
{
    if(s->pmix!=nullptr)
    v.assign(s->pmix->ckpt_clip(),s->pmix->ckpt_clip()+s->pmix->nf);
}

bool sediment_f::ckpt_rank_restore(lexer *p, ghostcell *pgc, const std::vector<std::vector<double>> &v, bool same)
{
    if(same && s->pmix!=nullptr && v.size()==1 && (int)v[0].size()==s->pmix->nf)
    std::copy(v[0].begin(),v[0].end(),s->pmix->ckpt_clip());
    return true;
}

double* sediment_f::ckpt_wln()
{
    return (pnhfsusp!=nullptr) ? pnhfsusp->ckpt_wln() : nullptr;
}
