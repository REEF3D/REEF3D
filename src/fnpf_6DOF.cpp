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

#include"fnpf_6DOF.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"solver.h"
#include"fnpf_laplace.h"
#include"fnpf_fsf.h"
#include"fnpf_bed_update.h"

fnpf_6DOF::fnpf_6DOF(lexer *p, fdm_fnpf *c, ghostcell *pgc) : initialized(false), foot(p), psiD(p), zeroslice(p)
{
    if(p->mpirank==0)
    cout<<"6DOF FNPF startup ..."<<endl;
    
    // one body, as in NHFLOW
    nbody = 1;
    
    for(int nb=0; nb<nbody; ++nb)
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));
    
    gcval = (p->j_dir==0) ? 150 : 250;
    footcount = 0;
    
    const int size = p->imax*p->jmax*(p->kmax+2);
    
    p->Darray(psi0,size);
    p->Darray(zero,size);
    
    for(int m=0; m<6; ++m)
    {
    psi[m] = nullptr;
    
    // 2D: sway, roll and yaw do not exist
    if(p->j_dir==0 && (m==1 || m==3 || m==5))
    continue;
    
    p->Darray(psi[m],size);
    }
    
    mark = new int[size]();
    
    pbed = new fnpf_bed_update(p);
}

fnpf_6DOF::~fnpf_6DOF()
{
}

void fnpf_6DOF::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->initialize_fnpf(p,c,pgc);
    
    geometry(p,c,pgc);
    extrapolate(p,c,pgc,c->Fi);
    
    initialized = true;
}

void fnpf_6DOF::forces(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_laplace *plap, fnpf_fsf *pf, 
                       slice &Keta, slice &Kfi, int iter)
{
    // psi_0: phi_t at fixed z on the free surface, phi_t|z = dFifsf/dt - Fz*deta/dt
    SLICELOOP4
    psiD(i,j) = Kfi(i,j) - c->Fz(i,j)*Keta(i,j);
    
    pgc->gcsl_start4(p,psiD,50);
    
    zero_face(p,c);
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->face_data_fnpf(p,c,pgc,-2,c->FBF,c->FBu,c->FBv,c->FBw);
    exchange_face(p,c,pgc);
    
    solve_psi(p,c,pgc,psolv,plap,pf,psi0,psiD);
    
    // added mass: six unit-mode solves, once per time step
    const bool refresh = (iter==0);
    
    for(int nb=0; nb<nbody; ++nb)
    {
        const bool computeA = refresh && !fb_obj[nb]->fnpf_fixed(p);
        
        if(computeA)
        for(int m=0; m<6; ++m)
        if(psi[m]!=nullptr)
        {
            zero_face(p,c);
            fb_obj[nb]->face_data_fnpf(p,c,pgc,m,c->FBF,c->FBu,c->FBv,c->FBw);
            exchange_face(p,c,pgc);
            
            solve_psi(p,c,pgc,psolv,plap,pf,psi[m],zeroslice);
        }
        
        fb_obj[nb]->forces_fnpf(p,c,pgc,psi0,psi,computeA);
    }
}

void fnpf_6DOF::motion(lexer *p, fdm_fnpf *c, ghostcell *pgc, int iter)
{
    // last stage: 2 for RK3 (A 310 3), 3 for RK4 (A 310 4)
    const bool finalize = (iter == ((p->A310==3) ? 2 : 3));
    
    for(int nb=0; nb<nbody; ++nb)
    {
        fb_obj[nb]->solve_eqmotion_fnpf(p,pgc,iter,finalize);
        fb_obj[nb]->update_position_fnpf(p,pgc,finalize);
        
        if(finalize)
        fb_obj[nb]->print_fnpf(p,pgc,iter);
    }
}

void fnpf_6DOF::footprint(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf, int gcval_eta, int gcval_fifsf)
{
    // harmonic extension of eta and Fifsf over the footprint (Gauss-Seidel, warm started
    // from the previous stage), the surrounding free surface acts as Dirichlet data
    if(footcount==0)
    return;
    
    for(int it=0; it<20; ++it)
    {
        SLICELOOP4
        if(foot(i,j)>0.5)
        {
            double se = eta(i-1,j) + eta(i+1,j);
            double sf = Fifsf(i-1,j) + Fifsf(i+1,j);
            double nn = 2.0;
            
            if(p->j_dir==1)
            {
            se += eta(i,j-1) + eta(i,j+1);
            sf += Fifsf(i,j-1) + Fifsf(i,j+1);
            nn += 2.0;
            }
            
            eta(i,j)   = se/nn;
            Fifsf(i,j) = sf/nn;
        }
        
        pgc->gcsl_start4(p,eta,gcval_eta);
        pgc->gcsl_start4(p,Fifsf,gcval_fifsf);
    }
}

void fnpf_6DOF::geometry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const int size = p->imax*p->jmax*(p->kmax+2);
    
    for(int n=0; n<size; ++n)
    c->FBF[n] = 0.0;
    
    SLICELOOP4
    foot(i,j) = 0.0;
    
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->ray_cast_fnpf(p,c,pgc,c->FBF,foot);
    
    pgc->gcparax7(p,c->FBF,7);
    pgc->gcsl_start4(p,foot,50);
    
    int count=0;
    SLICELOOP4
    if(foot(i,j)>0.5)
    ++count;
    
    footcount = pgc->globalisum(count);
    
    // Neumann data for phi: rigid-body velocity
    zero_face(p,c);
    for(int nb=0; nb<nbody; ++nb)
    fb_obj[nb]->face_data_fnpf(p,c,pgc,-1,c->FBF,c->FBu,c->FBv,c->FBw);
    exchange_face(p,c,pgc);
}

void fnpf_6DOF::extrapolate(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f)
{
    // two layers of body nodes next to the fluid: average of the fluid (or already
    // filled) neighbours; feeds the lagged cross terms of the Laplace rhs and the
    // sampling of the hull pressure
    const int size = p->imax*p->jmax*(p->kmax+2);
    const int sI = p->jmax*p->kmaxF;
    const int sJ = p->kmaxF;
    const double *FBF = c->FBF;
    
    for(int n=0; n<size; ++n)
    mark[n]=0;
    
    for(int pass=0; pass<2; ++pass)
    {
        ILOOP
        JLOOP
        FKLOOP
        {
            const int q = FIJK;
            
            if(p->flag7[q]>0 && FBF[q]>0.5 && mark[q]==0)
            {
                int nbr[6] = {q-sI, q+sI, q-sJ, q+sJ, q-1, q+1};
                const int nn = (p->j_dir==1) ? 6 : 4;
                
                if(p->j_dir==0)
                {
                nbr[2] = q-1;
                nbr[3] = q+1;
                }
                
                double sum=0.0;
                int cnt=0;
                
                for(int m=0; m<nn; ++m)
                {
                    const int r = nbr[m];
                    
                    if(p->flag7[r]>0 && (FBF[r]<0.5 || (mark[r]>0 && mark[r]<=pass)))
                    {
                    sum += f[r];
                    ++cnt;
                    }
                }
                
                if(cnt>0)
                {
                f[q] = sum/double(cnt);
                mark[q] = pass+1;
                }
            }
        }
        
        pgc->gcparax7(p,f,7);
    }
}

void fnpf_6DOF::solve_psi(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_laplace *plap, fnpf_fsf *pf, 
                          double *f, slice &D)
{
    // Dirichlet on the free surface outside the footprint
    FFILOOP4
    if(foot(i,j)<0.5)
    {
        f[FIJK]   = D(i,j);
        f[FIJKp1] = D(i,j);
        f[FIJKp2] = D(i,j);  
        f[FIJKp3] = D(i,j);
    }
    
    pbed->bedbc_sig(p,c,pgc,f,pf);
    pgc->start7V(p,f,c->bc,gcval);
    
    // the A329 inflow flux belongs to phi only
    double *uin = c->Uin;
    c->Uin = zero;
    
    plap->start(p,c,pgc,psolv,pf,f,D);
    
    c->Uin = uin;
    
    pgc->start7V(p,f,c->bc,gcval);
    
    extrapolate(p,c,pgc,f);
}

void fnpf_6DOF::zero_face(lexer *p, fdm_fnpf *c)
{
    const int size = p->imax*p->jmax*(p->kmax+2);
    
    for(int n=0; n<size; ++n)
    {
    c->FBu[n] = 0.0;
    c->FBv[n] = 0.0;
    c->FBw[n] = 0.0;
    }
}

void fnpf_6DOF::exchange_face(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    pgc->gcparax7(p,c->FBu,7);
    pgc->gcparax7(p,c->FBv,7);
    pgc->gcparax7(p,c->FBw,7);
}
