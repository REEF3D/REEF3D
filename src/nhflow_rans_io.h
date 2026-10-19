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

#ifndef NHFLOW_RANS_IO_H_
#define NHFLOW_RANS_IO_H_

#include"nhflow_turbulence.h"
#include"nhflow_strain.h"
#include<fstream>
#include<sstream>

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

class nhflow_rans_io : public nhflow_turbulence, public nhflow_strain
{
public:
    nhflow_rans_io(lexer*,fdm_nhf*);
	virtual ~nhflow_rans_io();
    
    void print_2D(lexer*, fdm_nhf*, ghostcell*,ofstream&,int) override final;
    void print_3D(lexer*, fdm_nhf*, ghostcell*, std::vector<char>&, size_t&) override final;
    void ini(lexer*, fdm_nhf*, ghostcell*) override final;
    void inflow(lexer*, fdm_nhf*, ghostcell*);
    void inflow_profile(lexer*, double, double, double&, double&, double&);
    double kinval(int,int,int) override final;
    double epsval(int,int,int) override final;
	void gcupdate(lexer*, fdm_nhf*, ghostcell*) override final;
	double ccipol_kinval(lexer*,ghostcell*,double,double,double) override final;
	double ccipol_epsval(lexer*,ghostcell*,double,double,double) override final;
    double ccipol_a_kinval(lexer*,ghostcell*,double,double,double) override final;
	double ccipol_a_epsval(lexer*,ghostcell*,double,double,double) override final;
    void kinget(int,int,int,double) override final;
    void epsget(int,int,int,double) override final;

    void name_ParaView_parallel(lexer*, ofstream&) override final;
    void name_ParaView(lexer*,stringstream&, int*, int &) override final;
    void offset_ParaView(lexer*, int*, int &) override final;
    
    void name_pvtp(lexer*, fdm_nhf*, ghostcell*,ofstream&) override final;
    void name_vtp(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &) override final;
    void offset_ParaView_2D(lexer*, int*, int &) override final;
    
    void ckpt_fields(const std::function<void(const char*, char, double*)> &fd, const std::function<void(const char*, char, int*)> &fi) override
    {
        fd("turb_KIN",'C',KIN); fd("turb_EPS",'C',EPS); fi("turb_WETN",'S',WETN);
        fd("turb_PK",'C',PK); fd("turb_PK0",'C',PK0); fd("turb_PK_b",'C',PK_b);
    }
    double *KIN,*EPS;
    int *WALLF;
    int *WETN;   // wet state of the column after the last turbulence step (newly wetted columns)
    
    void wetting_ini(lexer*, fdm_nhf*, double*, double*, bool);
    void length_limit(lexer*, fdm_nhf*, bool);
    
    // SST/Bradshaw limiter (A 564 2), shared by k-epsilon and k-omega
    void walldist_h(lexer*, fdm_nhf*, ghostcell*);
    double sst_walldist(lexer*, fdm_nhf*);
    double sst_F2(lexer*, fdm_nhf*, double, double);
    double *WDH = nullptr;   // distance to the nearest side wall with a wall function (A 519 2), walldist_h
	
	double const ke_c_1e, ke_c_2e,ke_sigma_k,ke_sigma_e;
	double const kw_alpha, kw_beta,kw_sigma_k,kw_sigma_w;

private:
    void kepsini_default(lexer*,fdm_nhf*,ghostcell*);

	float ffn;
	int q,iin,ii,jj,kk;
	int gcval_kin,gcval_eps,gcval_edv;

	double M,I,tau,H,B,ks,shearvel,kinbed,epsbed,omegabed;
    double uvel,refwalldist,fc,ev_fac, beddist, dist;
	double kinw,epsw;
	double walld,ddn,depth;
    
};

#endif


