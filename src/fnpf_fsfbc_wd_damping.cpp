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

#include"fnpf_fsfbc_wd.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"solver2D.h"

void fnpf_fsfbc_wd::damping(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f, int gcval, double alpha)
{
    double starttime=pgc->timer();
    
    int check=0;
    double alpha_eta=1.0;
    
    // Fifsf damping
    if((gcval==60 || gcval==160) && (p->A357==1 || p->A357==2))
    {
    check=1;
    alpha_eta=1.0;
    }
    
    // eta damping
    if((gcval==55 || gcval==155) && (p->A357==1 || p->A357==3))
    {
    check=1;
    alpha_eta=p->A366;
    }
    
    // With vb==0 everywhere the implicit damping system is diagonal with
    // solution f, so the global 2D BiCGStab solve is skipped. This is the
    // normal state whenever no cell is breaking.
    if(p->A350>0 && check==1)
    {
        double vbmax=0.0;
        SLICELOOP4
        vbmax = MAX(vbmax,fabs(c->vb(i,j)));
        
        vbmax = pgc->globalmax(vbmax);
        
        if(vbmax<=0.0)
        check=0;
    }
    
    if(p->A350>0 && check==1)
    {
        n=0;
        SLICELOOP4
        {
            if(p->wet[IJ]==1 || p->A343>=1)
            {
             visc = alpha_eta*c->vb(i,j);
                
             c->N.p[n] =   visc/(p->DXP[IM1]*p->DXN[IP])*p->x_dir
                         + visc/(p->DXP[IP]*p->DXN[IP])*p->x_dir
                         + visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                         + visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                           
                           + 1.0/(alpha*p->dt);
            
             c->rvec.V[n] =   f(i,j)/(alpha*p->dt);
             
             c->N.s[n] = -visc/(p->DXP[IM1]*p->DXN[IP])*p->x_dir;
             c->N.n[n] = -visc/(p->DXP[IP]*p->DXN[IP])*p->x_dir;
             
             c->N.e[n] = -visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;
             c->N.w[n] = -visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
            }

            if(p->wet[IJ]==0 && p->A343>=1)
            {
             c->N.p[n] =  1.0;
            
             // A343 2/3: keep eta/Fifsf of temporarily dry cells (eta=0 would
             // put a full water column under a dry flag, released on rewetting)
             c->rvec.V[n] =   (p->A343>=2) ? f(i,j) : 0.0;
             
             c->N.s[n] = 0.0;
             c->N.n[n] = 0.0;
             
             c->N.e[n] = 0.0;
             c->N.w[n] = 0.0;
            }
     
         ++n;
        }
        
        
        // walls: zero gradient of the new value (the diagonal takes the coupling). Before, the coupling went to the
        // right-hand side with the old value f(i,j), which damped the change of the cells next to a wall less: in a
        // Ting & Kirby flume 4 cells across, the surf zone kept 7 % more wave height than in 2D. Dry neighbours keep
        // the lagged form, and so does mesh refinement (G 1): with the implicit wall, fnpf_3d_amr_break_romero (A 350 3
        // on patches near the coast) went unstable at t = 6.6 s
        const bool implicit_wall = (p->G1==0 && p->amr_patch==0);
        
        n=0;
        SLICELOOP4
        {
            if(p->wet[IJ]==1 || p->A343>=1)
            {
                if(p->flagslice4[Im1J]<0 && implicit_wall)
                {
                c->N.p[n] += c->N.s[n];
                c->N.s[n] = 0.0;
                }
                else
                if(p->flagslice4[Im1J]<0 || (p->wet[Im1J]==0 && p->A343>=1))
                {
                c->rvec.V[n] -= c->N.s[n]*f(i,j);
                c->N.s[n] = 0.0;
                }
                
                if(p->flagslice4[Ip1J]<0 && implicit_wall)
                {
                c->N.p[n] += c->N.n[n];
                c->N.n[n] = 0.0;
                }
                else
                if(p->flagslice4[Ip1J]<0 || (p->wet[Ip1J]==0 && p->A343>=1))
                {
                c->rvec.V[n] -= c->N.n[n]*f(i,j);
                c->N.n[n] = 0.0;
                }
                
                if(p->flagslice4[IJm1]<0 && implicit_wall)
                {
                c->N.p[n] += c->N.e[n];
                c->N.e[n] = 0.0;
                }
                else
                if(p->flagslice4[IJm1]<0 || (p->wet[IJm1]==0 && p->A343>=1))
                {
                c->rvec.V[n] -= c->N.e[n]*f(i,j);
                c->N.e[n] = 0.0;
                }
                
                if(p->flagslice4[IJp1]<0 && implicit_wall)
                {
                c->N.p[n] += c->N.w[n];
                c->N.w[n] = 0.0;
                }
                else
                if(p->flagslice4[IJp1]<0 || (p->wet[IJp1]==0 && p->A343>=1))
                {
                c->rvec.V[n] -= c->N.w[n]*f(i,j);
                c->N.w[n] = 0.0;
                }
            }
     
        ++n;
        }
        
        pgc->gcsl_start4(p,f,gcval);
        
        // the rows with viscosity only (Jacobi, same residual), BiCGStab over all columns if that does not converge
        if(!damping_local(p,c,pgc,f))
        psolv->start(p,pgc,f,c->N,c->xvec,c->rvec,4);
        
        double time=pgc->timer()-starttime;
        
        if(p->mpirank==0 && p->count%p->P12==0 && p->D21==1)
        cout<<"fsfbc_damping: "<<p->solveriter<<"  fsfbc_damping_time: "<<setprecision(3)<<time<<endl;
    }
}

// The damping system (1/(alpha dt) + L_vb) f_new = f/(alpha dt): only the rows with breaking or coastline
// viscosity vb > 0 (A 346) couple to their neighbours. Every other row is diagonal: a wet cell keeps f, a
// dry cell takes its right-hand side. BiCGStab iterated over all columns of the domain for the few rows of
// the breaking zones and the coastline band; here the diagonal rows are set directly and the coupled rows
// are solved by Jacobi iterations (the rows are diagonally dominant by 1/(alpha dt): convergence factor
// kappa/(1+kappa), kappa = 4 vb alpha dt/dx^2) until the residual is 1e-10 of the initial one. (The test of
// BiCGStab, |r|_2/cellnumtot2D < N 44, is met by the old f on large grids; BiCGStab still does one
// iteration, which nearly solves this system. Jacobi with that test would do none.)
// As in BiCGStab the iterate is a separate vector with the halo exchange of the parallel boundaries, and f
// gets the solution in its interior cells. Without convergence after N 46 iterations f is not changed
// and the caller solves with BiCGStab.
bool fnpf_fsfbc_wd::damping_local(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    const int nij = p->imax*p->jmax;
    
    if(int(damp_x.size())!=nij)
    damp_x.assign(nij,0.0);
    
    double *const x = damp_x.data();
    const int sI = p->jmax;
    
    // diagonal rows directly, coupled rows into the list
    damp_act.clear();
    
    n=0;
    SLICELOOP4
    {
        if(c->N.s[n]!=0.0 || c->N.n[n]!=0.0 || c->N.e[n]!=0.0 || c->N.w[n]!=0.0)
        {
            damp_act.push_back({IJ,n});
            x[IJ] = f(i,j);
        }
        else
        x[IJ] = c->rvec.V[n]/c->N.p[n];
        
        ++n;
    }
    
    const int na = int(damp_act.size());
    damp_r.resize(na);
    
    p->solveriter = 0;
    double residual0 = 0.0;
    
    while(true)
    {
        pgc->gcslparaxijk_single(p,x,4);
        
        double sum = 0.0;
        
        for(int a=0; a<na; ++a)
        {
            const int q = damp_act[a].q;
            const int m = damp_act[a].n;
            
            const double r = c->rvec.V[m] - (c->N.p[m]*x[q]
                                           + c->N.n[m]*x[q+sI]
                                           + c->N.s[m]*x[q-sI]
                                           + c->N.w[m]*x[q+1]
                                           + c->N.e[m]*x[q-1]);
            damp_r[a] = r;
            sum += r*r;
        }
        
        const double residual = sqrt(pgc->globalsum(sum));
        
        if(p->solveriter==0)
        residual0 = residual;
        
        if(residual<=1.0e-10*residual0)
        break;
        
        if(p->solveriter>=p->N46)
        return false;
        
        for(int a=0; a<na; ++a)
        x[damp_act[a].q] += damp_r[a]/c->N.p[damp_act[a].n];
        
        ++p->solveriter;
    }
    
    SLICELOOP4
    f(i,j) = x[IJ];
    
    return true;
}
