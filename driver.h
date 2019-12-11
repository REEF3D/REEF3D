/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"increment.h"

class field;
class printer;
class initialize;
class diffusion;
class fdm;
class fdm2D;
class fdm_fnpf;
class lexer;
class momentum;
class ioflow;
class pressure;
class poisson;
class convection;
class turbulence;
class solver;
class ghostcell;
class timestep;
class freesurface;
class reini;
class particlecorr;
class sediment;
class bedload;
class suspended;
class topo;
class reinitopo;
class potential;
class heat;
class benchmark;
class sixdof;
class data;
class concentration;
class fnpf_fg;
class fnpf_sg;
class onephase;
class nsewave;
class nhflow_fsf;
class sflow;
class fnpf_vtu3D;
class fnpf_timestep;
class grid;

#include<iostream>
#include<fstream>
#include<iomanip>


#ifndef DRIVER_H_
#define DRIVER_H_

using namespace std;

class driver : public increment
{
public:

	driver(int&,char**);
	virtual ~driver();
    
    void start();
    
    void cfd_driver();
	void nsewave_driver();
    void nhflow_driver();
    void pfsg_driver();
    void pffg_driver();
    void sf_driver();
    
	void loop_cfd_fsi(fdm*);
	void loop_cfd(fdm*);
    void loop_nsewave(fdm*);
    void loop_nhflow(fdm*);
    void loop_fnpf_fg(fdm*);
    void loop_fnpf_sg();
    void loop_sflow(fdm*);
    
	void logic();
    void logic_fnpf_fg();
    void logic_fnpf_sg();
    void logic_sflow();
    
	void driver_ini();
    void driver_ini_nhflow();
    void driver_ini_fnpf_sg();
    void driver_ini_fnpf_fg();
    void driver_ini_sflow();
    
	void log_ini();
	void mainlog(lexer*);
	void maxlog(lexer*);
	void solverlog(lexer*);
	void sedimentlog(lexer*);
    
	void makegrid(lexer*,ghostcell*);
	void makegrid_cds();
    void makegrid2D(lexer*,ghostcell*);
    void makegrid2D_cds(lexer*,ghostcell*,fdm2D*);
    void makegrid_fnpf(lexer*,ghostcell*);
    void makegrid_fnpf_cds(lexer*,ghostcell*);
    void makegrid_nhflow(lexer*,ghostcell*);    
    
	void fill_vel(lexer*,fdm*,ghostcell*);
	void vec_test(lexer*,fdm*,ghostcell*,field&);
	void func_test(lexer*,fdm*,ghostcell*,field&);
    void mgc_test(lexer*,fdm*,ghostcell*);
	double calc();
    
    void stop(lexer*,fdm*,ghostcell*);

	printer* pprint;
	initialize* pini;
	diffusion* pdiff;
	diffusion* pturbdiff;
	diffusion* pconcdiff;
	diffusion* psuspdiff;
	fdm* a;
    fdm2D* b;
    fdm_fnpf *c;
	lexer* p;
	momentum* pmom;
	ioflow* pflow;
	pressure* ppress;
	poisson* ppois;
	convection* pconvec;
	convection* pturbdisc;
	convection* pfsfdisc;
	convection* pconcdisc;
    convection* pheatdisc;
	turbulence* pturb;
	solver* psolv;
	solver* ppoissonsolv;
    solver* plapsolv;
	ghostcell* pgc;
	timestep* ptstep;
	freesurface* pfsf;
	reini* preini;
	particlecorr* ppart;
	sediment* psed;
	bedload* pbed;
	suspended* psusp;
	topo* ptopo;
	reinitopo* preto;
	heat* pheat;
	potential* potflow;
	benchmark* pbench;
	sixdof* p6dof;
	data *pdata;
	concentration *pconc;
    fnpf_sg *ppfsg;
    fnpf_fg *ppffg;
    onephase *poneph;
    nsewave *pnse;
    nhflow_fsf *pnhfsf;
    sflow *psflow;
    fnpf_vtu3D *pfprint; 
    fnpf_timestep* pftstep;
    grid *pgrid;


private:
    double starttime, endtime;
    ofstream mainlogout;
    ofstream maxlogout;
    ofstream solvlogout;
    ofstream sedlogout;
	
	double nom,val;
};

#endif
