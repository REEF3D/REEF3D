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

#include"lexer.h"

lexer::lexer() : cmu(0.09), position(this), interpolation(this), coordinates(this)
{
    sigT=0.9;
    veclength=0;
    vec2Dlength=0;
    
	control::ini_default();

    solveriter=0;
	mpirank=0;

	simtime=0.0;
	poissontime=0.0;
	pressval=0;
    alpha=0.0;
    solidread=toporead=porousread=topoforcing=0;
    gridgeo=nullptr;
    geobed=nullptr;
    net_count=0;
    mooring_count=0;
    for(int n=0; n<8; ++n)
    gcbl_count[n]=-2;

    // iowave Riemann / Flather edges: set by iowave on the base grid only (AMR patch grids have none)
    open_xm=open_xp=open_ym=open_yp=0;

    // in- and outflow (set by the discharge of ioflow, used before it by the wave generation),
    // irregular wave parameters (set for irregular waves only, B 526-528 default to wTp otherwise),
    // coordinate transformation flag
    Ui=Ua=Uo=0.0;
    Hi=Ho=0.0;
    Qi=Qo=0.0;
    wT=wV=wH=wA=wL=wd=ww=wk=wC=0.0;
    wHs=wAs=wwp=ww_s=ww_e=wTp=wLp=0.0;
    wN=0;
    cms_flag=0;

    // free-surface levels at the in- and outflow: set by ini_phi (F 60-62), read by the outflow
    // pressure (B 77 1) also without a free surface
    fsfin=fsfout=0.0;
    printcount=printcount_sixdof=0;
    maxlength=0.0;      // largest domain extent: set by initialize (CFD), makegrid_sigma (NHFLOW, FNPF), sflow_maxcoor

    // extra ghost cells per side: never set, read by vecsize (gcbextra) and reduced over the ranks
    gcextra1=gcextra2=gcextra3=gcextra4=gcextra4a=gcextra6=0;

    // direct-forcing slice ghost cells: allocated by gridini (CFD, NHFLOW) or gcsldf_update, not for SFLOW
    gcsldfeta4_count=gcsldfbed4_count=0;
    gcsldfeta4=gcsldfbed4=nullptr;
}

lexer::lexer(const lexer &b, int) : control(b), position(this), interpolation(this), coordinates(this), cmu(0.09)
{
    sigT=0.9;
    veclength=0;
    vec2Dlength=0;

    solveriter=0;
	mpirank=b.mpirank;
    mpi_size=1;

	simtime=b.simtime;
	poissontime=0.0;
	pressval=0;
    alpha=0.0;
    solidread=toporead=porousread=topoforcing=0;
    gridgeo=nullptr;
    geobed=nullptr;
    net_count=0;
    mooring_count=0;

    count=b.count;
    dt=b.dt;
    dt_old=b.dt_old;
    umax=vmax=wmax=0.0;
    gctime=xtime=0.0;
    utime=vtime=wtime=0.0;
    recontime=fsftime=0.0;
    lsmtime=0.0;
    solver_error=0;
    wenofunc=nullptr;
    for(int n=0; n<8; ++n)
    gcbl_count[n]=-2;

    // iowave Riemann / Flather edges: set by iowave on the base grid only (AMR patch grids have none)
    open_xm=open_xp=open_ym=open_yp=0;

    // the state of the base grid that the patch grids (AMR, FAS) read but do not set themselves:
    // in- and outflow (discharge of ioflow, used before it by the wave generation), wave parameters
    // (wave_lib, B 526-528 default to wTp), free-surface levels at the in- and outflow (ini_phi,
    // F 60-62; the outflow pressure B 77 1), largest domain extent (reinitialisation)
    Ui=b.Ui; Ua=b.Ua; Uo=b.Uo;
    Hi=b.Hi; Ho=b.Ho;
    Qi=b.Qi; Qo=b.Qo;
    wT=b.wT; wV=b.wV; wH=b.wH; wA=b.wA; wL=b.wL; wd=b.wd; ww=b.ww; wk=b.wk; wC=b.wC;
    wHs=b.wHs; wAs=b.wAs; wwp=b.wwp; ww_s=b.ww_s; ww_e=b.ww_e; wTp=b.wTp; wLp=b.wLp;
    wN=b.wN;
    cms_flag=0;     // coordinate transformation: set by the patch builder (reefamr3d)
    fsfin=b.fsfin; fsfout=b.fsfout;
    maxlength=b.maxlength;
    printcount=printcount_sixdof=0;

    // extra ghost cells per side: never set, read by vecsize (gcbextra) and reduced over the ranks
    gcextra1=gcextra2=gcextra3=gcextra4=gcextra4a=gcextra6=0;

    // direct-forcing slice ghost cells: allocated by gridini (CFD, NHFLOW) or gcsldf_update, not for SFLOW
    gcsldfeta4_count=gcsldfbed4_count=0;
    gcsldfeta4=gcsldfbed4=nullptr;
}

lexer::~lexer()
{
    curv_logical_free();
}
