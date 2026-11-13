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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_timestep.h"
#include<iomanip>
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_sz_weights.h"
#include<cmath>

fnpf_timestep::fnpf_timestep(lexer* p):epsi(1.0e-3),maxtimestep(p->N49),c0_orig(p->N47)
{
}

fnpf_timestep::~fnpf_timestep()
{
}

// N 48: the time step of FNPF
//
//  0  fixed       dt = N 49
//  1  (default)   RK3, RK4: dt = N47 dx / MAX(|u|, sqrt(g h_max)), the shallow-water wave speed of the deepest column
//                 theta-Robin (A 310 5): the step of an accuracy target, the phase-speed error N 55 [%] of the reference
//                 wave (period N 52, or of the wave input), at most 4 x the stability step of RK3 at N 47 0.75, with a
//                 coastline (A 343 >= 1) or a floating body (X 10) at most that step (theta_step). The semi-implicit
//                 scheme is stable for every gravity wave, so its step follows from accuracy, as in SFLOW and NHFLOW
//  2  stability   dt = N47 C_RK / (omega_d + pi |u|/dx + pi |v|/dy), omega_d = sqrt(g D(k)) of the grid cut-off
//                 mode with the discrete Dirichlet-to-Neumann operator D of the sigma grid (Laplace cds2 and the
//                 free-surface sz stencil of A311), C_RK = sqrt(3) (RK3), 2 sqrt(2) (RK4): N47 is the
//                 fraction of the linear stability limit of the RK scheme (N47 < 1), for any depth, 2D and 3D.
//                 theta-Robin: N47 times the stability step of RK3 (N47 > 1 allowed)
//
// rule(): 1, 2 as above, 3 the accuracy step of theta-Robin (N 48 1 with A 310 5; N 48 0 with A 310 5 for dt_crit).
//
// N 52 T_min > 0, N 53 n: accuracy cap dt <= T_min/n (RK3, RK4). Once the stability limit is no longer the
// constraint, the RK damping of the short waves sets the accuracy (RK3: about (pi/12)(omega dt)^3 per period).
//
// In potential flow no wave travels at sqrt(g h) on the grid scale: the fastest discrete mode is the shortest
// wave the grid resolves, k = pi/dx, with omega = sqrt(g k tanh(k h)): in deep water the step of N 48 2 grows with
// sqrt(h/dx) against that of N 48 1.
int fnpf_timestep::rule(lexer *p)
{
    if(p->A310==5)
    return p->N48==2 ? 2 : 3;
    
    return p->N48==2 ? 2 : 1;
}

std::vector<double> fnpf_timestep::tab_lq, fnpf_timestep::tab_g;
double fnpf_timestep::tab_lqmin=-4.0, fnpf_timestep::tab_lqmax=6.0, fnpf_timestep::tab_dlq=0.0, fnpf_timestep::gmax_ratio=0.0;
int fnpf_timestep::dtn_ok=0;

double fnpf_timestep::rk_limit(lexer *p)
{
    return (p->A310==4) ? 2.0*sqrt(2.0) : sqrt(3.0);
}

// N 48 2: omega of the grid cut-off mode, kh2 = kx^2 + ky^2 the symbol of the horizontal cds2 Laplacian
// at the Nyquist mode (4/dx^2 on a uniform grid), h the local water depth
double fnpf_timestep::omega_discrete(lexer *p, double kh2, double h)
{
    if(dtn_ok==0)
    dtn_table_ini(p);
    
    const double g = fabs(p->W22)>0.0 ? fabs(p->W22) : 9.81;
    h = MAX(h,1.0e-10);
    
    const double gq = MAX(dtn_lookup(sqrt(kh2)*h), 0.0);
    
    return sqrt(g*gq/h);
}

// theta-Robin (A 310 5), N 48 1: theta-Robin is stable for every gravity wave, so the step is set by accuracy.
// The relative phase-speed error of a wave of period T with omega dt = W is, for theta 1/2 and the Newmark beta
// of A 308 (beta 1/4 without: Crank-Nicolson), 1 - phi/W with 4 sin^2(phi/2) = W^2/(1 + beta W^2) (exact for
// the linear oscillator; CN: phi = 2 atan(W/2), about W^2/12). dt is the largest step whose error at the
// reference wave is at most N 55 per cent, at most 4 x the RK3 step of N 48 2 at N 47 0.75 (3 dt_crit, the
// largest step tested; the corrector, breaking and the relaxation zones). The reference period T: N 52 if set,
// otherwise that of the wave input (regular waves: T, irregular: Tp). Without one: 2 x the RK3 step. Reconstructed
// waves (B 92 51-53) have none here: their peak period (fnpf_wave_period.h) would let the step grow past the
// frequencies of a focused packet.
double fnpf_timestep::theta_tref(lexer *p)
{
    if(p->N52>0.0)
    return p->N52;
    
    if(p->B90>0)
    return (p->B92>11) ? p->wTp : p->wT;
    
    return 0.0;
}

double fnpf_timestep::theta_error(lexer *p, double dt, double cu, double T)
{
    const double W = 2.0*PI*dt/T;
    double beta = 0.25;
    
    if(p->A308>0.0 && cu>0.0)
    {
    const double om = rk_limit(p)*dt/cu;
    beta = MIN(0.25, MAX(1.0/12.0, 0.25 - 1.0/(p->A308*p->A308*om*om)));
    }
    
    const double s = MIN(1.0, W/(2.0*sqrt(1.0 + beta*W*W)));
    const double phi = 2.0*asin(s);
    
    return fabs(1.0 - phi/W);
}

double fnpf_timestep::theta_step(lexer *p, double cu)
{
    // at most 4 x the RK3 step; with a coastline (A 343 >= 1) at most the RK3 step at N 47 0.75 (the shallowest
    // gauges grow with the step, and the moving coastline A 343 2 fails at 2 x at the wetted front; larger steps
    // with a static coastline: N 48 2 with N 47 1.5-2.25). With resolved floating bodies (X 10 1) also at most the
    // RK3 step: the natural periods of the body set the step there.
    const double dtmax = (p->A343>=1 || p->X10>0) ? 0.75*cu : 3.0*cu;
    const double T = theta_tref(p);
    
    if(T<=0.0 || cu<=0.0)
    return MIN(1.5*cu, dtmax);
    
    const double eps = 0.01*MAX(p->N55,1.0e-6);
    
    if(theta_error(p,dtmax,cu,T)<=eps)
    return dtmax;
    
    double a=0.0, b=dtmax;
    for(int q=0;q<60;++q)
    {
    const double m = 0.5*(a+b);
    if(theta_error(p,m,cu,T)<=eps)
    a = m;
    else
    b = m;
    }
    
    return a;
}

void fnpf_timestep::accuracy_cap(lexer *p)
{
    // theta-Robin, N 48 1: N 52 is the reference period of the accuracy target (theta_step), not a cap with N 53
    if(p->N52>0.0 && rule(p)!=3)
    p->dt = MIN(p->dt, p->N52/MAX(p->N53,1.0));
}

// weights w[m] of f[knoz-m] such that sz = sum_m w[m] f[knoz-m] (sigma units), as the sz() of the free-surface
// discretisation chosen by A311: Fornberg weights on the sigma nodes (fnpf_sz_weights)
int fnpf_timestep::sz_weights(lexer *p, double *w)
{
    const int N = p->knoz;
    int np = 3;
    
    if(p->A311==2 || p->A311==3 || p->A311==4 || p->A311==5)
    np = 5;
    
    if(p->A311==6)
    np = 7;
    
    while(np>3 && N<np-1)
    np -= 2;
    
    fnpf_sz_weights szw(p);
    
    const double *ww = szw.weights(np);
    
    if(ww==nullptr && np>3)
    ww = szw.weights(np=3);
    
    if(ww==nullptr)
    ww = szw.weights(np=2);
    
    for(int m=0;m<np;++m)
    w[m] = ww[m];
    
    return np;
}

// dimensionless discrete DtN g(q) = H D for q = kH on the sigma grid (flat bed): the vertical rows of
// fnpf_laplace_cds2,  q^2 phi_k - [ (phi_k+1 - phi_k)/DZN_k - (phi_k - phi_k-1)/DZN_k-1 ]/DZP_k-1 = 0,
// the bed ghost as KBEDBC (mt += mb), phi_knoz = 1, then g = sz(phi)
double fnpf_timestep::dtn_discrete(lexer *p, double q)
{
    const int N = p->knoz;
    std::vector<double> lo(N), di(N), up(N), rhs(N), phi(N+1);
    double w[7];
    
    for(int kk=0;kk<N;++kk)
    {
        const int m = kk+marge;
        const double at = 1.0/(p->DZP[m-1]*p->DZN[m]);
        const double ab = 1.0/(p->DZP[m-1]*p->DZN[m-1]);
        
        di[kk]  = q*q + at + ab;
        up[kk]  = -at;
        lo[kk]  = -ab;
        rhs[kk] = 0.0;
        
        if(kk==0)
        {
        up[kk] += lo[kk];
        lo[kk]  = 0.0;
        }
        
        if(kk==N-1)
        {
        rhs[kk] = at;   // phi_N = 1
        up[kk]  = 0.0;
        }
    }
    
    for(int kk=1;kk<N;++kk)
    {
        const double f = lo[kk]/di[kk-1];
        di[kk]  -= f*up[kk-1];
        rhs[kk] -= f*rhs[kk-1];
    }
    
    phi[N] = 1.0;
    phi[N-1] = rhs[N-1]/di[N-1];
    for(int kk=N-2;kk>=0;--kk)
    phi[kk] = (rhs[kk] - up[kk]*phi[kk+1])/di[kk];
    
    const int np = sz_weights(p,w);
    
    double g=0.0;
    for(int m=0;m<np;++m)
    g += w[m]*phi[N-m];
    
    return g;
}

void fnpf_timestep::dtn_table_ini(lexer *p)
{
    const int nq = 801;
    tab_lqmin = -4.0;
    tab_lqmax =  6.0;
    tab_dlq   = (tab_lqmax-tab_lqmin)/double(nq-1);
    
    tab_lq.resize(nq);
    tab_g.resize(nq);
    
    int nneg=0;
    gmax_ratio=0.0;
    
    for(int n=0;n<nq;++n)
    {
        tab_lq[n] = tab_lqmin + n*tab_dlq;
        const double q = pow(10.0,tab_lq[n]);
        tab_g[n] = dtn_discrete(p,q);
        
        if(tab_g[n]<=0.0)
        ++nneg;
        
        gmax_ratio = MAX(gmax_ratio, tab_g[n]/(q*tanh(q)));
    }
    
    dtn_ok=1;
    
    if(p->mpirank==0)
    {
        cout<<"FNPF time step: discrete DtN of the sigma grid, knoz "<<p->knoz<<", A311 "<<p->A311<<endl;
        cout<<"        kH      g_disc   kH tanh(kH)"<<endl;
        
        const double qs[7]={0.1,1.0,3.0,10.0,30.0,100.0,1000.0};
        for(int n=0;n<7;++n)
        cout<<setw(10)<<qs[n]<<setw(12)<<setprecision(5)<<dtn_lookup(qs[n])<<setw(12)<<qs[n]*tanh(qs[n])<<endl;
        
        cout<<"max g_disc/(kH tanh kH): "<<setprecision(4)<<gmax_ratio<<endl;
        
        if(nneg>0)
        cout<<"WARNING: discrete DtN <= 0 for "<<nneg<<" table entries, check the sigma grid and A311"<<endl;
    }
}

double fnpf_timestep::dtn_lookup(double q)
{
    if(q<=0.0)
    return 0.0;
    
    const double lq = log10(q);
    
    // shallow limit g ~ q^2, deep limit bounded by linear growth
    if(lq<=tab_lqmin)
    return tab_g[0]*pow(q/pow(10.0,tab_lqmin),2.0);
    
    if(lq>=tab_lqmax)
    return tab_g.back()*q/pow(10.0,tab_lqmax);
    
    int n = int((lq-tab_lqmin)/tab_dlq);
    n = MIN(n, int(tab_g.size())-2);
    
    const double r = (lq - tab_lq[n])/tab_dlq;
    
    return (1.0-r)*tab_g[n] + r*tab_g[n+1];
}

void fnpf_timestep::start(fdm_fnpf *c, lexer *p,ghostcell *pgc)
{
    double depthmax=0.0;

    p->umax=p->vmax=p->wmax=p->viscmax=irsm=jrsm=krsm=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
	p->dt_old=p->dt;

	p->umax=p->vmax=p->wmax=p->viscmax=0.0;
	sqd=1.0/(p->DXM*p->DXM);
    
    // k
    k=p->knoz;

// maximum velocities
    SLICELOOP4
	depthmax=MAX(depthmax,c->depth(i,j));

	SLICELOOP4
    if(i+p->origin_i>=5 || p->j_dir==1)
    {
	// the surface node (k = knoz): of the 3D velocities, or the surface ones (fnpf_memory.h)
	p->umax=MAX(p->umax,fabs(c->U!=nullptr ? c->U[FIJK] : c->Us(i,j)));
    p->vmax=MAX(p->vmax,fabs(c->V!=nullptr ? c->V[FIJK] : c->Vs(i,j)));
    p->wmax=MAX(p->wmax,fabs(c->W!=nullptr ? c->W[FIJK] : c->Ws(i,j)));
    }
    
    SLICELOOP4
	p->viscmax=MAX(p->viscmax, c->vb(i,j));
    
    // one MPI_Allreduce instead of five
    double red[5] = {depthmax, p->umax, p->vmax, p->wmax, p->viscmax};
    pgc->globalmax(red,5);
    depthmax  = red[0];
    p->umax   = red[1];
    p->vmax   = red[2];
    p->wmax   = red[3];
    p->viscmax= red[4];

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    }
	
    cu=cv=cw=1.0e10;
    
    
    

    if(p->mpirank==0 && (p->count%p->P12==0) && p->viscmax>0.0)
	cout<<"viscmax: "<<p->viscmax<<endl;
    

    SLICELOOP4
    {
    if(p->j_dir==1 && p->knoy>1)
    dx = MIN(p->DXN[IP],p->DYN[JP]);
    
    if(p->j_dir==0 || p->knoy==1)
    dx = p->DXN[IP];
    
    if(rule(p)>=2)
    {
    // discrete: the cut-off mode of the local grid in the local water depth, Doppler shift by the surface velocity
    if(p->wet[IJ]==0)
    continue;
    
    const int ydir = (p->j_dir==1 && p->knoy>1) ? 1 : 0;
    const double kx2 = (2.0/p->DXP[IP] + 2.0/p->DXP[IM1])/p->DXN[IP];
    const double ky2 = ydir==1 ? (2.0/p->DYP[JP] + 2.0/p->DYP[JM1])/p->DYN[JP] : 0.0;
    
    const double rate = omega_discrete(p,kx2+ky2,c->WL(i,j)) + PI*p->umax/p->DXN[IP] + (ydir==1 ? PI*p->vmax/p->DYN[JP] : 0.0);
    
    cu = MIN(cu, rk_limit(p)/(rate+1.0e-20));
    continue;
    }
    
    cu = MIN(cu, 1.0/((fabs(MAX(p->umax, sqrt(9.81*depthmax)))/dx)));
    
    if(p->j_dir==1 )
    cv = MIN(cv, 1.0/((fabs(MAX(p->vmax, sqrt(9.81*depthmax)))/dx)));
    }

    if(p->j_dir==1 )
    cu = MIN(cu,cv);
    
    // minimum over all ranks
    cu = pgc->globalmin(cu);
    
    // no water and no velocity: the time step would be N47*1e10 and the run would end silently
    if(cu>=1.0e10)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! FNPF: no time step, the domain holds no water (max water depth "<<depthmax<<" m, umax "<<p->umax<<" m/s) !!!"<<endl<<endl;
        
        pgc->final(true);
    }

   	p->dt=p->N47*cu;
    p->dt_crit=cu;
    
    // theta-Robin, N 48 1: the step of the accuracy target instead of N 47 times the stability step
    if(rule(p)==3)
    {
        p->dt=theta_step(p,cu);
        
        static bool said=false;
        if(!said && p->mpirank==0)
        {
        const double T = theta_tref(p);
        cout<<"FNPF time step N 48 1 (theta-Robin): phase-speed error N 55 "<<p->N55<<" % of the reference wave T = "<<T<<" s ("
            <<(p->N52>0.0 ? "N 52" : (T>0.0 ? "wave input" : "none: 2 x the RK3 step"))<<"), dt = "<<p->dt<<" s = "<<p->dt/(0.75*cu)
            <<" x the RK3 step of N 48 2, N 47 0.75"<<(p->A343>=1 || p->X10>0 ? " (at most 1 x with a coastline, A 343 >= 1, or a floating body, X 10)" : (p->dt>=3.0*cu*(1.0-1.0e-9) ? " (the bound of 4 x)" : ""))<<endl;
        }
        said=true;
    }
    
    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"dt: "<<p->dt<<endl;

    
    if (p->N48==0) 
    p->dt=maxtimestep;
    
    else
    {
	p->dt=MIN(p->dt,maxtimestep);
    
    // N 52, N 53: accuracy cap
    accuracy_cap(p);
    }
    
    
}

void fnpf_timestep::ini(fdm_fnpf* c, lexer* p,ghostcell* pgc)
{
    double depthmax,wdt;
    
	p->umax = p->vmax = p->wmax = -1e19;
    depthmax = -1e19;
    
    cu=cv=1.0e10;
    
    
    SLICELOOP4
	depthmax=MAX(depthmax,p->wd - c->bed(i,j));
	
	depthmax=pgc->globalmax(depthmax);
	
    k=p->knoz;
    SLICELOOP4
	p->umax=MAX(p->umax,fabs((c->Fifsf[Ip1J]-c->Fifsf[Im1J])/(p->DXP[IP]+p->DXP[IM1])));

	p->umax=pgc->globalmax(p->umax);


	k=p->knoz;
    SLICELOOP4
	p->vmax=MAX(p->vmax,fabs((c->Fifsf[IJp1]-c->Fifsf[IJm1])/(p->DYP[JP]+p->DYP[JM1])));

	p->vmax=pgc->globalmax(p->vmax);


	// (only printed; without the 3D velocities the surface ones, fnpf_memory.h)
	if(c->W!=nullptr)
	FLOOP
	p->wmax=MAX(p->wmax,fabs(c->W[FIJK]));
	
	if(c->W==nullptr)
	SLICELOOP4
	p->wmax=MAX(p->wmax,fabs(c->Ws(i,j)));

	p->wmax=pgc->globalmax(p->wmax);
	
	p->umax=MAX(p->umax,2.0*p->ufbmax);
	p->umax=MAX(p->umax,2.0*p->vfbmax);
	p->umax=MAX(p->umax,2.0*p->wfbmax);
    
    p->umax=MAX(p->umax,2.0*p->X210_u);
	p->umax=MAX(p->umax,2.0*p->X210_v);
	p->umax=MAX(p->umax,2.0*p->X210_w);

	p->dt=p->DXM/(p->umax+epsi);
    

    p->umax+=10.0;

 
    FLOOP
    FPWDCHECK
    {
    if(p->j_dir==1 && p->knoy>1)
    dx = MIN(p->DXN[IP],p->DYN[JP]);
    
    if(p->j_dir==0 || p->knoy==1)
    dx = p->DXN[IP];
    
    
    if(rule(p)>=2)
    {
    const int ydir = (p->j_dir==1 && p->knoy>1) ? 1 : 0;
    const double kx2 = (2.0/p->DXP[IP] + 2.0/p->DXP[IM1])/p->DXN[IP];
    const double ky2 = ydir==1 ? (2.0/p->DYP[JP] + 2.0/p->DYP[JM1])/p->DYN[JP] : 0.0;
    
    const double rate = omega_discrete(p,kx2+ky2,p->wd-c->bed(i,j)) + PI*p->umax/p->DXN[IP] + (ydir==1 ? PI*p->vmax/p->DYN[JP] : 0.0);
    
    cu = MIN(cu, rk_limit(p)/(rate+1.0e-20));
    continue;
    }
    
    cu = MIN(cu, 1.0/((fabs(MAX(p->umax, sqrt(9.81*depthmax)))/dx)));
    cv = MIN(cv, 1.0/((fabs(MAX(p->vmax, sqrt(9.81*depthmax)))/dx)));
    }

	cu = MIN(cu,cv);
    
   	p->dt=p->N47*cu;
    
    p->dt=pgc->globalmin(p->dt);
    p->dt_crit=pgc->globalmin(cu);
    
    // theta-Robin, N 48 1: the step of the accuracy target instead of N 47 times the stability step
    if(rule(p)==3)
    p->dt=theta_step(p,p->dt_crit);
    
    if(p->N48>=1)
    accuracy_cap(p);
    
    // N 48 0: the fixed time step N 49 from the first step on, as in start()
    if(p->N48==0)
    p->dt=maxtimestep;
    
	p->dt_old=p->dt;
    
    if(p->mpirank==0 && p->N48>0 && rule(p)>1)
    {
        if(rule(p)==3)
        cout<<"FNPF time step N 48 1 (theta-Robin): the step of the phase-speed error N 55 = "<<p->N55<<" % of the reference wave"<<endl;
        
        if(rule(p)==2)
        {
        cout<<"FNPF time step N 48 2: dt omega_max = N47 C_RK, C_RK = "<<setprecision(4)<<rk_limit(p)<<", N47 "<<p->N47<<endl;
        
        if(p->N47>1.0 && p->A310!=5)
        cout<<"WARNING: N47 > 1 is beyond the linear stability limit of the RK scheme with N 48 2"<<endl;
        }
    }
    
    if(p->mpirank==0 && p->N52>0.0 && rule(p)!=3)
    cout<<"FNPF time step accuracy cap: dt <= N52/N53 = "<<p->N52<<"/"<<p->N53<<" = "<<p->N52/MAX(p->N53,1.0)<<" s"<<endl;
    
    if(p->B94==0)
	wdt=p->phimean;

	if(p->B94==1)
	wdt=p->B94_wdt;
    
    // all columns dry: no wet cell set the time step above, it would be N47*1e10 and the run
    // would end after one step without a message (e.g. a bed level above the still water level
    // everywhere from a broken geo.dat interpolation)
    int nwet=0, ncol=0;
    double bedmin=1.0e20, bedmax=-1.0e20;
    
    SLICELOOP4
    {
    ++ncol;
    
    if(p->wet[IJ]>0)
    ++nwet;
    
    bedmin=MIN(bedmin,c->bed(i,j));
    bedmax=MAX(bedmax,c->bed(i,j));
    }
    
    nwet=pgc->globalisum(nwet);
    ncol=pgc->globalisum(ncol);
    bedmin=pgc->globalmin(bedmin);
    bedmax=pgc->globalmax(bedmax);
    
    if(nwet==0)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! FNPF: all "<<ncol<<" columns are dry: bed level "<<bedmin<<" to "<<bedmax
            <<" m, still water level F 60 = "<<p->phimean<<" m, wetting criterion A 344 = "<<p->A344<<" m."
            <<" Check the bed level of the grid (DIVEMesh geo.dat, G 51/G 52, G 23) and F 60 !!!"<<endl<<endl;
        
        pgc->final(true);
    }
    
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    cout<<"depth_max: "<<setprecision(3)<<depthmax<<" |  wave depth: "<<wdt<<endl;
    }
}

