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

#include"fnpf_laplace_cds2.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"solver.h"
#include"fnpf_bed_update.h"

fnpf_laplace_cds2::fnpf_laplace_cds2(lexer *p) 
{
    pbed = new fnpf_bed_update(p);
    
    gcval=250;
    if(p->j_dir==0)
    gcval=150;
}

fnpf_laplace_cds2::~fnpf_laplace_cds2()
{
}

void fnpf_laplace_cds2::start(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, double *f, slice &Fifsf)
{
    p->poissoniter=0;
    p->poissontime=0.0;
    
    starttime=pgc->timer();
    
    // Fused single-pass assembly (previously two full LOOP sweeps).
    // Column-constant factors are hoisted; per-cell arithmetic is kept in
    // the original evaluation order, so M and rhs are bit-identical.
    const int sI = p->jmax*p->kmaxF;
    const int sJ = p->kmaxF;
    
    const int    *const __restrict flag7 = p->flag7;
    const double *const __restrict sigx  = p->sigx;
    const double *const __restrict sigy  = p->sigy;
    const double *const __restrict sigxx = p->sigxx;
    const double *const __restrict Uin   = c->Uin;
    
    double *const __restrict Mp = c->M.p.data();
    double *const __restrict Mn = c->M.n.data();
    double *const __restrict Ms = c->M.s.data();
    double *const __restrict Mw = c->M.w.data();
    double *const __restrict Me = c->M.e.data();
    double *const __restrict Mt = c->M.t.data();
    double *const __restrict Mb = c->M.b.data();
    double *const __restrict rhs = c->rhsvec.V.data();
    
    const double ydir = p->y_dir;
    const double xdir = p->x_dir;
    
    // The mixed derivatives 2*sigx*d2Fi/dxdsig (and y) are not in the 7-point
    // matrix but in rhs, evaluated with the current f. Without further
    // iterations they are those of the previous RK stage (or time step), a
    // lag of one stage: O(dt) and growing with the wave slope (sigx ~ Ex).
    // With A324>0 the assembly and solve are repeated (Picard) with the new
    // iterate until the change is below A325 times the change of the first
    // solve, at most A324 times. A324 0 is the previous behaviour.
    const int nouter = 1 + MAX(p->A324,0);
    double dfirst=0.0;
    
    if(nouter>1 && int(fold.size())!=p->imax*p->jmax*(p->kmax+2))
    fold.resize(p->imax*p->jmax*(p->kmax+2));
    
    for(int qo=0; qo<nouter; ++qo)
    {
    if(qo>0)
    {
    // bed ghost cells and halos from the new iterate
    pbed->bedbc_sig(p,c,pgc,f,pf);
    pgc->start7V(p,f,c->bc,gcval);
    
    starttime=pgc->timer();
    }
    
    if(nouter>1)
    FLOOP
    fold[FIJK] = f[FIJK];
    
	n=0;
    ILOOP
    JLOOP
    {
        k=0;
        const int c0 = FIJK;
        
        const int wP = p->wet[IJ];
        const int wS = p->wet[Im1J];
        const int wN = p->wet[Ip1J];
        const int wE = p->wet[IJm1];
        const int wW = p->wet[IJp1];
        const int bcS = c->bc(i-1,j);
        const int bcN = c->bc(i+1,j);
        
        const double dxn = 1.0/(p->DXP[IP]*p->DXN[IP]);
        const double dxs = 1.0/(p->DXP[IM1]*p->DXN[IP]);
        const double dyw = 1.0/(p->DYP[JP]*p->DYN[JP])*ydir;
        const double dye = 1.0/(p->DYP[JM1]*p->DYN[JP])*ydir;
        const double hxy = dxn + dxs + dyw + dye;
        
        const double dxc = p->DXP[IP]+p->DXP[IM1];
        const double dyc = p->DYP[JP]+p->DYP[JM1];
        
        const double sz = p->sigz[IJ];
        const double bx = c->Bx(i,j);
        const double by = c->By(i,j);
        
        KLOOP
        if(p->flag4[IJK]>0)
        {
            const int q = c0 + k;
            
            if(wP==1 && flag7[q]>0)
            {
            const double sx  = sigx[q];
            const double sy  = sigy[q];
            const double sxx = sigxx[q];
            const double s2  = sx*sx + sy*sy + sz*sz;
            const double dzc = p->DZN[KP]+p->DZN[KM1];
            
            double mp = hxy + (s2/(p->DZP[KM1]*p->DZN[KP])) + (s2/(p->DZP[KM1]*p->DZN[KM1]));
            double mn = -dxn;
            double ms = -dxs;
            double mw = -dyw;
            double me = -dye;
            double mt = -(s2/(p->DZP[KM1]*p->DZN[KP])  + sxx/dzc);
            double mb = -(s2/(p->DZP[KM1]*p->DZN[KM1]) - sxx/dzc);
            
            double rv = 2.0*sx*(f[q+sI+1] - f[q-sI+1] - f[q+sI-1] + f[q-sI-1])/(dxc*dzc)
                      + 2.0*sy*(f[q+sJ+1] - f[q-sJ+1] - f[q+sJ-1] + f[q-sJ-1])/(dyc*dzc)*ydir;
            
            // KBEDBC
            if(flag7[q-1]<0)
            {
            const double ab = mb;
            denom = sz + bx*sx + by*sy;
            
                mn +=  ab*2.0*p->DZN[KP]*bx/(denom*dxc);
                ms += -ab*2.0*p->DZN[KP]*bx/(denom*dxc);
                
                mw +=  ab*2.0*p->DZN[KP]*by/(denom*dyc);
                me += -ab*2.0*p->DZN[KP]*by/(denom*dyc);
                
                mt += ab;
                mb = 0.0;
            }
            
            // south
            const bool sdry = (flag7[q-sI]<0 || wS==0);
            
            if(p->B98<=2)
            if(sdry)
            {
            mp += ms;
            ms = 0.0;
            }
            
            if(p->B98>2)
            {
            if(sdry && bcS==0)
            {
            mp += ms;
            ms = 0.0;
            }
            
            if(flag7[q-sI]<0 && bcS==1  && p->A329==1)
            {
            rv += ms*Uin[q-sI]*p->DXP[IM1];
            mp += ms;
            ms = 0.0;
            }
            
            if(flag7[q-sI]<0 && bcS==1  && p->A329>=2)
            {
            denom = -1.5*p->XP[IM1] + 2.0*p->XP[IP] - 0.5*p->XP[IP1];
            
            rv += (2.0/3.0)*ms*Uin[q-sI]*denom;
            mp += (4.0/3.0)*ms;
            mn -= (1.0/3.0)*ms;
            ms = 0.0;
            }
            }
            
            // north
            const bool ndry = (flag7[q+sI]<0 || wN==0);
            
            if(p->B99<=2)
            if(ndry)
            {
            mp += mn;
            mn = 0.0;
            }
            
            if(p->B99>2)
            {
            if(ndry && bcN==0)
            {
            mp += mn;
            mn = 0.0;
            }
            
            if(flag7[q+sI]<0 && bcN==2  && p->A329==1)
            {
            rv -=  2.0*sx*(f[q+sI+1] - f[q-sI+1] - f[q+sI-1] + f[q-sI-1])
                        /(dxc*dzc)*xdir;
                        
            rv +=  2.0*sx*(Uin[q+sI+1] - Uin[q+sI-1])
                        /(dzc)*xdir;
                        
            rv -= mn*Uin[q+sI]*p->DXP[IP1];
            mp += mn;
            mn = 0.0;
            }
            
            if(flag7[q+sI]<0 && bcN==2  && p->A329>=2)
            {
            rv -=  2.0*sx*(f[q+sI+1] - f[q-sI+1] - f[q+sI-1] + f[q-sI-1])
                        /(dxc*dzc)*xdir;
                        
            rv +=  2.0*sx*(Uin[q+sI+1] - Uin[q+sI-1])
                        /(dzc)*xdir;
                        
            denom = -0.5*p->XP[IM1] + 2.0*p->XP[IP] - 1.5*p->XP[IP1];
            
            rv += (2.0/3.0)*mn*Uin[q+sI]*denom;
            mp += (4.0/3.0)*mn;
            ms -= (1.0/3.0)*mn;
            mn = 0.0;
            }
            }

            // east
            if(flag7[q-sJ]<0 || wE==0)
            {
            mp += me;
            me = 0.0;
            }

            // west
            if(flag7[q+sJ]<0 || wW==0)
            {
            mp += mw;
            mw = 0.0;
            }
            
            // FSFBC
            if(flag7[q+2]<0 && flag7[q+1]>0)
            {
            rv -= mt*f[q+2];
            mt = 0.0;
            }
            
            Mp[n] = mp;
            Mn[n] = mn;
            Ms[n] = ms;
            Mw[n] = mw;
            Me[n] = me;
            Mt[n] = mt;
            Mb[n] = mb;
            rhs[n] = rv;
            }
            
            else
            if(wP==0 || flag7[q]<0)
            {
            Mp[n] = 1.0;
            Mn[n] = 0.0;
            Ms[n] = 0.0;
            Mw[n] = 0.0;
            Me[n] = 0.0;
            Mt[n] = 0.0;
            Mb[n] = 0.0;
            rhs[n] = 0.0;
            }
            
        ++n;
        }
    }
    
    endtime=pgc->timer();
    //if(p->mpirank==0 && (p->count%p->P12==0))
	//cout<<"LAPLCE_time: "<<endtime-starttime<<endl;
    
    starttime=pgc->timer();
    psolv->startF(p,pgc,f,c->rhsvec,c->M,8);
    endtime=pgc->timer();
    
    p->poissoniter+=p->solveriter;
    p->poissontime+=endtime-starttime;
    
    if(nouter>1)
    {
    double dmax=0.0;
    
    FLOOP
    dmax = MAX(dmax,fabs(f[FIJK]-fold[FIJK]));
    
    dmax = pgc->globalmax(dmax);
    
    if(qo==0)
    dfirst = dmax;
    
    if(qo>0 && dmax<=p->A325*dfirst)
    break;
    }
    }
    
    p->laplacetime+=p->poissontime;
    
    
    
	if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"Fi_iter: "<<p->poissoniter;
    //cout<<" Final_residual: "<<p->final_res;
    cout<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
    }
}

