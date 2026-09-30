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

#include"net_interface.h"
#include"net_membrane.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"nhflow_membrane_beta.h"
#include"vrans_definitions.h"
#include<mpi.h>
#include<fstream>
#include<sstream>
#include<string>

// ctrl.txt: X 330 1 (A 520 1 or 2). membrane.dat (read by rank 0, broadcast):
//
//   # comment
//   membrane box       x0 x1 y0 y1 z_bottom z_top     starts a new membrane (in 2D y0, y1 are ignored)
//   membrane cylinder  xc yc R z_bottom z_top
//   name        text                   optional label
//   resistance  R_n [R_t]              hydraulic resistance [m/s], leakage u_n = (dp/rho)/R_n; default 1e4, 0
//   thickness   delta                  half width of the smeared layer [m]; default 1.5 max(dx,dy,dz)
//   mesh        h                      target triangle edge length [m]; default min(dx,dy)
//   fill        dh                     initial inner water level above the outside level [m]; default 0
//   print       dt                     vtp output interval [s]; default none
//   floorpressure 0|1|3                static pressure below the floor: 0 uniform head difference,
//                                      1 local excess head (comparison only), 3 shape of the
//                                      time-averaged ramp; default 3 with A 520 1, 0 with A 520 2
//   tau         t                      averaging time of floorpressure 3 [s]; default 2
//   projections n                      projection passes per stage (default 1: Rhie-Chow continuity flux)
//   poisson     0|1                    membrane mobility in the pressure Poisson equation; default 1
//                                      (0 only to demonstrate the splitting leakage of the projection)
//
// Parameter lines apply to the most recent 'membrane' line.

void net_interface::membrane_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    string content;
    int len=0;

    if(p->mpirank==0)
    {
        ifstream f("membrane.dat");

        if(!f)
        len=-1;
        else
        {
            stringstream ss;
            ss<<f.rdbuf();
            content = ss.str();
            len = (int)content.size();
        }
    }

    MPI_Bcast(&len,1,MPI_INT,0,pgc->mpi_comm);

    if(len<0)
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330: membrane.dat not found in the case directory !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    content.resize(len);
    if(len>0)
    MPI_Bcast(&content[0],len,MPI_CHAR,0,pgc->mpi_comm);

    // parse
    vector<membrane_param> mp;
    istringstream is(content);
    string line;
    int lineno=0;
    bool error=false;

    while(getline(is,line))
    {
        ++lineno;

        const size_t hash = line.find('#');
        if(hash!=string::npos)
        line = line.substr(0,hash);

        istringstream ls(line);
        string key;

        if(!(ls>>key))
        continue;

        if(key=="membrane")
        {
            string shape;
            membrane_param m;
            ls>>shape;

            if(shape=="box")
            {
                m.shape=1;
                if(!(ls>>m.x0>>m.x1>>m.y0>>m.y1>>m.zb>>m.zt))
                error=true;
            }
            else if(shape=="cylinder")
            {
                m.shape=2;
                if(!(ls>>m.xc>>m.yc>>m.R>>m.zb>>m.zt))
                error=true;
            }
            else
            error=true;

            m.name = "membrane"+to_string(mp.size());
            mp.push_back(m);
        }
        else if(mp.empty())
        error=true;
        else if(key=="name")
        {
            if(!(ls>>mp.back().name))
            error=true;
        }
        else if(key=="resistance")
        {
            if(!(ls>>mp.back().Rn))
            error=true;
            ls>>mp.back().Rt;
        }
        else if(key=="thickness")
        {
            if(!(ls>>mp.back().delta))
            error=true;
        }
        else if(key=="mesh")
        {
            if(!(ls>>mp.back().h))
            error=true;
        }
        else if(key=="fill")
        {
            if(!(ls>>mp.back().fill))
            error=true;
        }
        else if(key=="floorpressure")
        {
            if(!(ls>>mp.back().floorp))
            error=true;
        }
        else if(key=="tau")
        {
            if(!(ls>>mp.back().tau) || mp.back().tau<=0.0)
            error=true;
        }
        else if(key=="projections")
        {
            if(!(ls>>mp.back().projections) || mp.back().projections<1)
            error=true;
        }
        else if(key=="poisson")
        {
            if(!(ls>>mp.back().poisson))
            error=true;
        }
        else if(key=="print")
        {
            if(!(ls>>mp.back().printdt))
            error=true;
        }
        else
        error=true;

        if(error)
        {
            if(p->mpirank==0)
            cout<<"\n!!! membrane.dat, line "<<lineno<<": cannot read '"<<line<<"' !!!\n"<<endl;
            MPI_Abort(pgc->mpi_comm,1);
        }
    }

    if(mp.empty())
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330: no membrane defined in membrane.dat !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    // mobility field for the pressure Poisson equation, 1 away from membranes
    if(d->MBETA==nullptr)
    p->Darray(d->MBETA,p->imax*p->jmax*(p->kmax+2));

    if(d->MCHI==nullptr)
    p->Darray(d->MCHI,p->imax*p->jmax*(p->kmax+2));

    if(d->MRCX==nullptr)
    p->Darray(d->MRCX,p->imax*p->jmax*(p->kmax+2));

    if(d->MRCY==nullptr)
    p->Darray(d->MRCY,p->imax*p->jmax*(p->kmax+2));

    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    {
    d->MBETA[qn]=1.0;
    d->MCHI[qn]=0.0;
    }

    for(size_t m=0; m<mp.size(); ++m)
    d->MPROJ = MAX(d->MPROJ, mp[m].projections);

    for(size_t m=0; m<mp.size(); ++m)
    {
        pmem.push_back(new net_membrane(m,mp[m]));
        pmem.back()->initialize_nhflow(p,d,pgc);
        pmem.back()->fill_nhflow(p,d,pgc);
    }
}

void net_interface::membrane_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha,
                                            double *UH, double *VH, double *WH, slice &WL)
{
    // 1. mobility field beta (also builds the membrane cell maps for this stage)
    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    d->MBETA[qn]=1.0;

    for(auto m : pmem)
    m->mobility_nhflow(p,d,pgc,alpha);

    pgc->start4V(p,d->MBETA,1);

    // 2. static overpressure of the bag below its floor (prescribed pressure, see net_membrane)
    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    d->MCHI[qn]=0.0;

    for(auto m : pmem)
    m->static_pressure_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);

    pgc->start4V(p,d->MCHI,1);

    // 3. implicit porous-jump forcing
    //    Incremental pressure scheme (A 520 2): the predictor already contains the old non-hydrostatic
    //    pressure gradient, -a C/rho G P^n. For the projection to stay consistent, P^n has to act with
    //    the same mobility as the pressure correction (beta_h horizontal, beta vertical), not with the
    //    implicit forcing factor (I + a A)^-1 in the layer or 1 next to it. It is taken out before the
    //    forcing and put back with the correction mobility afterwards; without this the pressure
    //    increments accumulate in the layer and the scheme diverges.
    const bool incremental = (p->A520==2);

    if(incremental)
    membrane_pgrad(p,d,alpha,UH,VH,WH,WL,1);

    for(auto m : pmem)
    m->forcing_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);

    if(incremental)
    membrane_pgrad(p,d,alpha,UH,VH,WH,WL,-1);
}

void net_interface::membrane_pgrad(lexer *p, fdm_nhf *d, double alpha, double *UH, double *VH, double *WH, slice &WL, int mode)
{
    // Incremental scheme (A 520 2): the predictor contains the old pressure gradient -a C/rho G P^n
    // (nhflow_pjm_corr::upgrad, vpgrad, wpgrad: plain wide gradient). Next to the membrane it has to act
    // through the same operator as the pressure correction, i.e. with the face mobilities
    // (nhflow_membrane_gradx/grady) and beta on the vertical and sigma terms, not with the implicit
    // forcing factor (I + a A)^-1 in the layer or with mobility 1 next to it. Then P^n and PCORR act
    // like one full pressure P^{n+1} and the scheme is the projection of A 520 1.
    //   mode  1: take the plain gradient out before the forcing
    //   mode -1: put it back with the mobilities after the forcing
    double gx,gy,gz,sx,sy,a;
    const double *P = d->P;

    LOOP
    WETDRYDEEP
    {
        if(!nhflow_membrane_active(p,d,i,j,k))
        continue;

        a = alpha*p->dt*CPORNH;
        const double bv = d->MBETA[IJK];
        const double dPk = (P[FIJKp1]-P[FIJK])/p->DZN[KP];

        sx = 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*dPk;
        sy = 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*dPk;

        if(mode==1)
        {
            gx = (0.5*(P[FIp1JKp1]+P[FIp1JK])-0.5*(P[FIm1JKp1]+P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]) + sx;
            gy = p->j_dir==1 ? (0.5*(P[FIJp1Kp1]+P[FIJp1K])-0.5*(P[FIJm1Kp1]+P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]) + sy : 0.0;
            gz = dPk;
        }
        else
        {
            gx = -(nhflow_membrane_gradx(p,d,P,i,j,k) + bv*sx);
            gy = p->j_dir==1 ? -(nhflow_membrane_grady(p,d,P,i,j,k) + bv*sy) : 0.0;
            gz = -bv*dPk;
        }

        // wpgrad acts on WH without the water depth: dW = -a C/rho dP/(DZN WL)
        gx *= 1.0/p->W1;
        gy *= 1.0/p->W1;
        gz *= 1.0/(p->W1*MAX(WL(i,j),1.0e-20));

        d->U[IJK] += a*gx;
        UH[IJK]   += a*gx*WL(i,j);

        if(p->j_dir==1)
        {
        d->V[IJK] += a*gy;
        VH[IJK]   += a*gy*WL(i,j);
        }

        d->W[IJK] += a*gz;
        WH[IJK]   += a*gz*WL(i,j);
    }
}

void net_interface::membrane_reaction_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, slice &WL, bool finalize)
{
    for(auto m : pmem)
    m->reaction_nhflow(p,d,pgc,alpha,WL,finalize);
}
