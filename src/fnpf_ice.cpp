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

#include"fnpf_ice.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ice_contact_nscd.h"

fnpf_ice::fnpf_ice(lexer *p, fdm_fnpf *c, ghostcell *pgc) : nfloe(0),nobst(0),etat(p),dtot(p),pcontact(nullptr),comm(pgc->mpi_comm),
                                                        stage(0),nstage(3),stagewarn(0),printtime(0.0),printcount(0)
{
    rhow  = p->W1;
    g     = fabs(p->W22);
    wd    = p->wd;
    alpha = MAX(p->A381,0.0);
    zeta  = MAX(p->A382,0.0);
    Cd    = MAX(p->A383,0.0);
    nsub  = MAX(p->A386,1);
    is2D  = (p->j_dir==0) ? 1 : 0;
    
    // taper width from the global mean cell size, identical on all ranks
    double dxs=0.0, dxn=0.0;
    SLICELOOP4
    {
    dxs += is2D ? p->DXN[IP] : 0.5*(p->DXN[IP] + p->DYN[JP]);
    dxn += 1.0;
    }
    dxs = pgc->globalsum(dxs);
    dxn = pgc->globalsum(dxn);
    dxmean = dxs/MAX(dxn,1.0);
    taper = MAX(p->A389,0.0)*dxmean;
    
    // breaking
    breakflag = MAX(0,MIN(3,p->A390));
    sigf   = p->A391;
    KIC    = p->A392_K;
    Csplit = p->A392_C;
    Dmin   = (p->A393>0.0) ? p->A393 : 4.0*dxmean;
    ndir   = is2D ? 1 : MAX(1,p->A395_dir);
    noff   = MAX(1,p->A395_off);
    nbreak = 0;
    nstage = (p->A310==4) ? 4 : 3;

    // the lid stiffens the surface under the floes: g_eff = g*(1+alpha)
    dtfac = 1.0/sqrt(1.0+alpha);

    SLICELOOP4
    {
    etat(i,j)=0.0;
    dtot(i,j)=0.0;
    }
    
    read(p,pgc);

    Yn.resize(floe.size());
    D1.resize(floe.size());
    D2.resize(floe.size());
    D3.resize(floe.size());

    double ymin = is2D ? 0.0 : p->global_ymin;
    double ymax = is2D ? 0.0 : p->global_ymax;

    pcontact = new ice_contact_nscd(p->A384_mu, p->A384_e, p->A388, p->A385, p->global_xmin, p->global_xmax, ymin, ymax);
}

fnpf_ice::~fnpf_ice()
{
    delete pcontact;

    if(logout.is_open())
    logout.close();

    if(obstout.is_open())
    obstout.close();
}

void fnpf_ice::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // floes float at their equilibrium draft rho_i*h/rho_w, the surface under the footprints is
    // depressed by the same amount: p = rho_i*g*h balances -rho_w*g*eta, Fi stays at rest
    for(auto &fl : floe)
    kinematics(fl);

    footprint(p,c);

    // an initial offset of the floe (disp in ice_floes.dat) moves the surface under it along,
    // as for a body and its wetted surface pushed down together in a decay test
    SLICELOOP4
    c->eta(i,j) -= dtot(i,j);
    
    for(auto &e : cell)
    {
    i=e.i;
    j=e.j;
    const fnpf_ice_floe &fl = floe[e.f];
    const double draft = fl.rho*fl.h/rhow;
    c->eta(i,j) += e.phi*(zbottom(fl,e.xc,e.yc) + draft - wd);
    }

    pgc->gcsl_start4(p,c->eta, p->j_dir==0 ? 155 : 55);

    print_ini(p);
    print(p,c,pgc);

    if(p->mpirank==0)
    {
    cout<<"FNPF ice: "<<nfloe<<" floes, "<<nobst<<" obstacles"<<endl;
    cout<<"FNPF ice: lid stiffness A381 = "<<alpha<<", damping ratio A382 = "<<zeta<<", drag Cd A383 = "<<Cd<<endl;
    cout<<"FNPF ice: time step factor 1/sqrt(1+A381) = "<<dtfac<<endl;
    cout<<"FNPF ice: footprint edge taper +-"<<taper<<" m"<<endl;
    
    if(breakflag>0)
    cout<<"FNPF ice: breaking "<<((breakflag&1)?"flexural sigma_f = ":"")<<((breakflag&1)?to_string(sigf):"")
        <<((breakflag&2)?"  splitting F = C*K_IC*h*sqrt(D), C*K_IC = ":"")<<((breakflag&2)?to_string(Csplit*KIC):"")
        <<"  D_min = "<<Dmin<<" m"<<endl;
    
    for(auto &fl : floe)
    if(fl.type==0)
    {
    const double draft = fl.rho*fl.h/rhow;
    if(draft > 0.5*MAX(taper,dxmean))
    {
    cout<<"FNPF ice: warning, floe "<<fl.id<<" draft "<<draft<<" m is steep over the edge taper, surface slope ~ "
        <<draft/(2.0*MAX(taper,dxmean))<<", consider a larger A 389"<<endl;
    break;
    }
    }
    }
}

void fnpf_ice::prestep(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    SLICELOOP4
    etat(i,j) = 0.0;

    for(size_t n=0; n<floe.size(); ++n)
    {
    kinematics(floe[n]);
    get_state(floe[n],Yn[n].data());
    }

    footprint(p,c);

    stage=0;
}

void fnpf_ice::store_etat(lexer *p, fdm_fnpf *c, slice &K)
{
    // deta/dt of the current stage, from the kinematic FSBC
    SLICELOOP4
    etat(i,j) = K(i,j);
}

void fnpf_ice::lid_forcing(lexer *p, fdm_fnpf *c, slice &K, slice &eta)
{
    ++stage;

    if(stage>nstage)
    {
        if(stagewarn==0 && p->mpirank==0)
        cout<<"FNPF ice: more dynamic FSBC calls per step than RK stages, floes are not advanced further"<<endl;
        stagewarn=1;
    }

    // floe geometry of this stage
    for(auto &fl : floe)
    if(fl.type==0)
    kinematics(fl);

    // lid pressure into the dynamic FSBC, stage loads on the floes
    stage_forces(p,c,K,eta);

    // advance the floes with the same RK stage as the free surface
    if(stage<=nstage)
    rk_stage(p);
}

void fnpf_ice::poststep(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const int checkbreak = (breakflag>0 && p->count%MAX(1,p->A394)==0);
    
    // bending moments from the last RK stage loads, before contact moves the floes
    if(checkbreak)
    breaking_moments(p);
    
    if(p->mpirank==0)
    contact(p);
    
    if(checkbreak && p->mpirank==0)
    breaking_decide(p);
    
    broadcast_state(p,pgc);
    
    if(checkbreak)
    breaking_apply(p,pgc);

    print(p,c,pgc);
}

void fnpf_ice::timestep(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->N48==0)
    {
        if(p->count==0 && p->mpirank==0 && alpha>0.0)
        cout<<"FNPF ice: fixed time step (N 48 0), make sure dt is reduced by 1/sqrt(1+A381) = "<<dtfac<<endl;
        return;
    }

    double dtice = p->dt*dtfac;

    double vmax=0.0;
    double omax=0.0;
    for(auto &fl : floe)
    if(fl.type==0)
    {
    vmax = MAX(vmax, sqrt(fl.v[0]*fl.v[0] + fl.v[1]*fl.v[1]) + fabs(fl.w[2])*fl.rbound);

    // explicit floe update: lid spring frequency of the footprint (heave, roll/pitch), damping 2*zeta*omega
    omax = MAX(omax, MAX(fl.omega, sqrt(fl.klid/(fl.rho*fl.h)))*MAX(1.0,2.0*zeta));
    }

    // floes should not move more than half a cell per step
    if(vmax>1.0e-10)
    dtice = MIN(dtice, 0.5*dxmean/vmax);

    if(omax>0.0)
    dtice = MIN(dtice, 1.0/omax);

    p->dt = pgc->globalmin(dtice);
}

void fnpf_ice::broadcast_state(lexer *p, ghostcell *pgc)
{
    if(p->mpi_size<=1)
    return;

    const int nv = 16;
    vector<double> buf(nv*floe.size(),0.0);

    if(p->mpirank==0)
    for(size_t n=0; n<floe.size(); ++n)
    {
        double *b = &buf[nv*n];
        const fnpf_ice_floe &fl = floe[n];
        get_state(fl,b);
        b[13] = fl.Fc[0];
        b[14] = fl.Fc[1];
        b[15] = double(fl.type);
    }

    MPI_Bcast(buf.data(), int(buf.size()), MPI_DOUBLE, 0, pgc->mpi_comm);

    if(p->mpirank>0)
    for(size_t n=0; n<floe.size(); ++n)
    {
        const double *b = &buf[nv*n];
        fnpf_ice_floe &fl = floe[n];
        set_state(fl,b);
        fl.Fc[0] = b[13];
        fl.Fc[1] = b[14];
        fl.type  = int(b[15]);
    }
}

void fnpf_ice::get_state(const fnpf_ice_floe &fl, double *y)
{
    for(int q=0;q<3;++q) y[q]    = fl.x[q];
    for(int q=0;q<4;++q) y[3+q]  = fl.q[q];
    for(int q=0;q<3;++q) y[7+q]  = fl.v[q];
    for(int q=0;q<3;++q) y[10+q] = fl.w[q];
}

void fnpf_ice::set_state(fnpf_ice_floe &fl, const double *y)
{
    for(int q=0;q<3;++q) fl.x[q] = y[q];
    for(int q=0;q<4;++q) fl.q[q] = y[3+q];
    for(int q=0;q<3;++q) fl.v[q] = y[7+q];
    for(int q=0;q<3;++q) fl.w[q] = y[10+q];
}
