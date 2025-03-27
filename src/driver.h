/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef DRIVER_H_
#define DRIVER_H_

class benchmark;
class concentration;
class convection;
class data;
class density;
class diffusion;
class fdm;
class fdm2D;
class fdm_fnpf;
class fdm_nhf;
class field;
class fnpf;
class fnpf_printer;
class fnpf_timestep;
class freesurface;
class fsi;
class ghostcell;
class grid;
class heat;
class initialize;
class ioflow;
class lexer;
class momentum;
class momentum_RKLS3_df;
class momentum_RKLS3_sf;
class multiphase;
class net;
class nhflow;
class nhflow_convection;
class nhflow_diffusion;
class nhflow_forcing;
class nhflow_fsf;
class nhflow_momentum;
class nhflow_potential;
class nhflow_pressure;
class nhflow_printer;
class nhflow_reconstruct;
class nhflow_scalar_convection;
class nhflow_signal_speed;
class nhflow_timestep;
class nhflow_turbulence;
class particle_corr;
class patchBC_interface;
class poisson;
class potential;
class pressure;
class printer;
class ptf;
class reini;
class reinitopo;
class sediment;
class sixdof;
class sflow;
class solver;
class timestep;
class turbulence;
class vrans;

#include"increment.h"
#include<fstream>
#include<vector>

using namespace std;

class driver : public increment
{
public:

    driver(int&,char**);
    virtual ~driver() = default;

private:
    
    void cfd_driver();
    void nhflow_driver();
    void fnpf_driver();
    void ptf_driver();
    void sflow_driver();
    
    void loop_cfd();
    void loop_cfd_df();
    void loop_cfd_sf();
    void loop_nhflow();
    void loop_ptf();
    void loop_fnpf();
    
	void logic_cfd();
    void logic_ptf();
    void logic_fnpf();
    void logic_nhflow();
    
    void patchBC_logic();
    
	void driver_ini_cfd();
    void driver_ini_nhflow();
    void driver_ini_fnpf();
    void driver_ini_ptf();
    void driver_ini_sflow();
    
	void log_ini();
	void mainlog(lexer*);
	void maxlog(lexer*);
	void solverlog(lexer*);
    
	void makegrid(lexer*,ghostcell*);
	void makegrid_cds();
    void makegrid2D(lexer*,ghostcell*);
    void makegrid2D_basic(lexer*,ghostcell*);
    void makegrid2D_cds(lexer*,ghostcell*,fdm2D*);
    void makegrid_sigma(lexer*,ghostcell*);
    void makegrid_sigma_cds(lexer*,ghostcell*);
    
	void vec_test(lexer*,fdm*,ghostcell*,field&);
	void func_test(lexer*,fdm*,ghostcell*,field&);
	double calc();
    
    void stop(lexer*,fdm*,ghostcell*);

    void assign_density();

    benchmark* pbench;
    concentration* pconc;
    convection* pconcdisc;
    convection* pconvec;
    convection* pfsfdisc;
    convection* pheatdisc;
    convection* pturbdisc;
    data* pdata;
    density* pd;
    diffusion* pdiff;
    diffusion* pturbdiff;
    diffusion* pconcdiff;
    diffusion* psuspdiff;
    freesurface* pfsf;
    fdm* a;
    fdm2D* b;
    fdm_fnpf* c;
    fdm_nhf* d;
    fnpf* ppfsg;
    fnpf_printer* pfprint;
    fnpf_timestep* pftstep;
    fsi* pfsi;
    ghostcell* pgc;
    grid* pgrid;
    heat* pheat;
    ioflow* pflow;
    initialize* pini;
    lexer* p;
    multiphase* pmp;
    momentum* pmom;
    momentum_RKLS3_df* pmom_df;
    momentum_RKLS3_sf* pmom_sf;
    std::vector<net*> pnet;
    nhflow* pnhf;
    nhflow_convection* pnhfconvec;
    nhflow_diffusion* pnhfdiff;
    nhflow_diffusion* pnhfturbdiff;
    nhflow_forcing* pnhfdf;
    nhflow_fsf* pnhfsf;
    nhflow_momentum* pnhfmom;
    nhflow_pressure* pnhpress;
    nhflow_potential* pnhfpot;
    nhflow_timestep* pnhfstep;
    nhflow_printer* pnhfprint;
    nhflow_reconstruct* precon;
    nhflow_scalar_convection* pnhfscalarconvec;
    nhflow_signal_speed* pss;
    nhflow_turbulence* pnhfturb;
    particle_corr* ppls;
    patchBC_interface* pBC;
    poisson* ppois;
    potential* potflow;
    pressure* ppress;
    printer* pprint;
    ptf* pptf;
    reini* preini;
    reinitopo* preso;
    reinitopo* preto;
    sediment* psed;
    sflow* psflow;
    sixdof* p6dof;
    solver* plapsolv;
    solver* ppoissonsolv;
    solver* psolv;
    turbulence* pturb;
    turbulence* pturbcfd;
    timestep* ptstep;
    vrans* pvrans;
    
    double starttime, endtime;
    std::ofstream mainlogout;
    std::ofstream maxlogout;
    std::ofstream solvlogout;
    
    double nom,val;
    char version[100];
};

#endif
