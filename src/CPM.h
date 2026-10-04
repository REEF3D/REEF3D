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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------

--------------------------------------------------------------------
CPM : Continuum Particle Method
--------------------------------------------------------------------*/

#ifndef CPM_H_
#define CPM_H_

#include"increment.h"
#include"part.h"
#include"slice4.h"
#include"field4a.h"
#include"boundarycheck.h"
#include"vtp3D.h"
#include<fstream>
#include<random>
#include<vector>

class lexer;
class fdm;
class ghostcell;
class sediment_fdm;
class turbulence;
class part;
class vrans;
class field;

using namespace std;

/*--------------------------------------------------------------------
MP-PIC with unresolved parcels and Eulerian (grid based) particle stresses

  dUp/dt = Dp (Uf - Up) - grad(p)/rho_p + g - grad(Ps)/(theta rho_p) + Ff

  Dp      : Andrews & O'Rourke (1996) drag, point implicit
  Ff      : Coulomb friction of the packed bed (Q 12 2, Q 13 1), implicit
  Ps      : particle normal stress
            Q 12 1 : Snider (2001)
            Q 12 2 : packed bed, see CPM_stress_packedbed.cpp
                     - effective stress of the contact network (Terzaghi)
                     - contact pressure against over-packing, Johnson & Jackson (1987)
                     - Coulomb friction with mu(I) rheology (Jop et al. 2006), Q 13 1
--------------------------------------------------------------------*/

class CPM : public increment, private vtp3D
{
public:
    CPM(lexer*, ghostcell*);
    virtual ~CPM() = default;

    void move(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*);
    
    void mppic_RK2(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*);
    void mppic_EE1(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*);
    
    void plain_RK2(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*);

    void update(lexer*, fdm*, ghostcell*, sediment_fdm*, field&, field&);
    void topo_update(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void bedzh_update(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void topo_column(lexer*, fdm*, ghostcell*);
    void topo_iso(lexer*, ghostcell*, field&);
    double ptopo(lexer*, fdm*, double, double, double);

    void timestep(lexer*, ghostcell*);

    void seed_particles(lexer*, fdm*, ghostcell*, sediment_fdm*);
    
    // hotstart and log
    void state_write(lexer*, int);
    void state_read(lexer*, ghostcell*, int);
    void sedlog(lexer*, ghostcell*);
    void ini_fields(lexer*, fdm*, ghostcell*, sediment_fdm*);
    
    // print
    void print_particles(lexer*,sediment_fdm*);
    void print_3D_CPM(lexer*, ghostcell*,  std::vector<char>&, size_t&);
    void name_ParaView_parallel_CPM(lexer*, ofstream&);
    void name_ParaView_CPM(lexer*, ostream&, int*, int &);
    void offset_ParaView_CPM(lexer*, int*, int &);
    
    
private:
    void advec_plain(lexer*, fdm*, part&, sediment_fdm*, turbulence*,
                        double*, double*, double*, double*, double*, double*,
                        double&, double&, double&, double);
                        
    void advec_mppic(lexer*, fdm*, part&, sediment_fdm*, turbulence*,
                        double*, double*, double*, double*, double*, double*,
                        double&, double&, double&, double);
                        
    void nearbed_velocity(lexer*, fdm*, double, double, double, double);
    double seed_diameter(lexer*, int);
    // two-way coupling
public:
    void coupling_update(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void fluid_forcing(lexer*, fdm*, ghostcell*, double, field&, field&, field&);
    
    // two-way coupling: mixture continuity div(u) = -div(theta u_p), source for the pressure Poisson equation
    inline static CPM *coupled = nullptr;
    void continuity_source(lexer*, ghostcell*);
    void coupling_reset(lexer*, ghostcell*);
    field4a Dsrc;
private:
    
    void limiter(lexer*, fdm*, ghostcell*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
    double occupancy_max(lexer*, ghostcell*);
    void substep_euler(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*, double);
    void substep_rk2(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*, double);
    double substep_size(lexer*, ghostcell*, double, int);
    void grid_update(lexer*, fdm*, ghostcell*, sediment_fdm*, double*, double*, double*, double*, double*, double*);

    // drag
    double drag_model(lexer *, double, double, double, double);

    void count_particles(lexer*, fdm*, ghostcell*, sediment_fdm*);

    // particle stress
    void stress_snider(lexer*, ghostcell*, sediment_fdm*);
    void stress_packedbed(lexer*, ghostcell*, sediment_fdm*);
    void stress_overburden(lexer*, ghostcell*, sediment_fdm*);
    void friction(lexer*, fdm*, double, double, double, double&, double&, double&, double, double);
    void gradient(lexer*, ghostcell*, field&, field&, field&, field&);
    double contact_pressure(double, double);
    double contact_pressure_deriv(double, double);
    void dilatancy(lexer*, ghostcell*);
    field4a T0e,Tiso;
    
    void stress_gradient(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void pressure_gradient(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void volfrac_update(lexer*, ghostcell*, sediment_fdm*, double*, double*, double*, double*, double*, double*);
    void smooth(lexer*, ghostcell*, field&, int);
    void kernel(lexer*, double, double, double);
    bool wallcell(lexer*, int, int, int);
    
    // periodic boundaries in x and y (DIVEMesh C 21, C 22)
    //   1 serial: CPM folds the ghost deposits and wraps the parcels
    //   2 parallel: the neighbour exchange sums the deposits, the parcels are wrapped on receipt
    int perx, pery;
    void pfold(lexer*, field&);
public:
    void periodic_flags(lexer*);
private:
    void periodic_wrap(lexer*, double*, double*, double*, double*);
    
    // turbulent dispersion (Q 52), random displacement with the eddy diffusivity of the fluid
    void dispersion_update(lexer*, fdm*, ghostcell*);
    void dispersion(lexer*, double, double, double, double&, double&, double&, double);
    field4a Kt,dKx,dKy,dKz;
    double Ktmax;
    std::mt19937_64 rng;
    std::normal_distribution<double> gauss;
    
    void wallbc(lexer*, ghostcell*, sediment_fdm*);

    void boundcheck(lexer*, int);

    part P;

    slice4 bedch;

    field4a Tau,Ts;
    field4a cellSum;
    field4a Us,Vs,Ws;
    field4a Pov;
    field4a Kc,KUx,KUy,KUz;
    field4a dSx,dSy,dSz;
    field4a Locc,Lout,Lin,Ltc,Lloc,LA,Lh0,Lh1,Lh2,Lh3;

    // relax
    void relax_ini(lexer*);
    void relax(lexer*, ghostcell*, sediment_fdm*);
    double rf(lexer*, double, double);
    double r1(lexer*, double, double);
    double distcalc(lexer*, double , double, double , double, double);
    
    double heaviside(double);
    double epsi,HS;

    void print_vtp(lexer*,sediment_fdm*);
    void pvtp(lexer*,int);

    boundarycheck boundaries;

    int printcount;
    double printtime;

    field4a dPx,dPy,dPz;
    field4a dTx,dTy,dTz;

    double *tan_betaQ73,*betaQ73,*dist_Q73;

    int timestep_ini = 0;
    int nsub, zsplit, nrej_step, nclip_step, nit_step=0;
    int restored, logini;
    int open_side[6];
    double outvol;
    ofstream logout;
    double dtsub,cmax;
    double hmin;
    
    // kernel
    int ki[2],kj[2],kk[2];
    double kw[2][2][2];
    
    // parameters
    double theta_max, theta_bed, theta_0;
    double Fr, eta0, eta1;
    double mu_s, mu_2, I0;
    double Vp;
    
    double Dpx,Dpy,Dpz;
    double dPx_val,dPy_val,dPz_val;
    double Bx,By,Bz;
    double uf,vf,wf;
    double liftx,lifty,liftz;
    // near-bed closure of the last call: weight of the exposed layer, log-law factor,
    // reference point of the fluid velocity, grid solid velocity (pore water)
    double nb_w,nb_fac,nb_xr,nb_yr,nb_zr,nb_ug,nb_vg,nb_wg;
    
    // Bagnold sheltering (Q 57): stress carried by the moving grains per bed column
    void bagnold_update(lexer*, fdm*, ghostcell*, sediment_fdm*);
    slice4 tauGf;
    int bag_count=-1;
    double shelter_factor(lexer*, sediment_fdm*, double, double, double, double, double, double);
    std::vector<double> tauG,tauB;
    double shelter=1.0;
    // exposure of the parcels of the bed (S 10 1): the top grain layer of each bed column
    void exposure_update(lexer*, fdm*);
    std::vector<double> expo;
    double expo_w=-1.0;
    // sub-grid bedload layer (Q 58, S 10 1), see CPM_bedload.cpp
    void bedload_columns(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void bedload_exchange(lexer*, fdm*, ghostcell*, sediment_fdm*, double);
    void bedload_move(lexer*, sediment_fdm*, int, double);
    bool bedload_grain(lexer*, int, int, double, double&, double&, double&);
    void bedload_column(lexer*, double, double, int&, int&);
    double settling_velocity(lexer*, double);
    void bedload_occupancy(lexer*);
    bool bedload_rest(lexer*, fdm*, int);
    double bedload_place(lexer*, double, double, double, int, int, double, double);
    slice4 blTx,blTy,blGx,blGy,blH,blC,blCs;
    int bl_npick=0, bl_ndep=0, bl_nsus=0;
    double Urel,Vrel,Wrel;
    double Tsval;
    double dTx_val,dTy_val,dTz_val;
    double DragCoeff,Fd;
    double F,G,H;
};

#endif
