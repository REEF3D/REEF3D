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

#include"lexer.h"
#include<cstdlib>
#include<cstring>
#include<iostream>
#include<mpi.h>

//  The grids of a mesh refinement (CFD, NHFLOW, FNPF, SFLOW, SEASTATE patches; the FAS grids of
//  SEASTATE) are lexers built from the base grid (level 0) with lexer(const lexer&, int).  Besides
//  the keys (control) and its own geometry (the builders: reefamr3d, reefamr, seastate_fas2) a grid
//  has the runtime state of the lexer below, in four groups.  Every runtime member of the lexer
//  outside the geometry belongs to one of them; a new member that a solver reads on a patch goes
//  into the right list here, nowhere else.
//
//  RUN    level 0 sets it during the run, the grids read it and do not set it: in- and outflow
//         (ioflow discharge, hydrographs), free-surface levels at the in- and outflow, the gage
//         pressure, the wave parameters and wave times, the sediment times.  Copied at the
//         construction and again before every step of a grid (amr_follow).
//  ONCE   level 0 sets it at its initialisation, constant afterwards: the extent of the domain,
//         time step limits.  Copied at the construction.
//  BUILD  the start value of level 0, then the grid's own: the mean water level (SFLOW sets its
//         own), the interface width (the CFD builder sets its own), the sediment step.  Copied at
//         the construction.
//  OWN    the grid's own: maxima of the step, solver statistics, timers, volumes, the bodies'
//         working values.  The defaults of lexer_ini / control_calc at the construction.
//  The time (simtime, dt, dt_old, count) is copied at the construction and set by the stepping of
//  the refinement (with subcycling every level has its own step).
//
//  IO     the in- and outflow members of RUN (discharge, inflow and outflow velocities, the levels at the
//         in- and outflow): RUN, except on a block of a multi-block grid with open faces (G 54,
//         amr_own_io = 1), whose own in- and outflow (ioflow_amr_block) sets them for the block.
//
//  The contract check (environment REEF3D_AMR_CHECK=1, set by the regression suite): amr_follow
//  stops the run when a grid's own objects have changed a RUN member since the last amr_follow
//  (it belongs to BUILD or OWN then), or when level 0 has changed a ONCE member after the grid was
//  built (it belongs to RUN).  Without the environment a violation is printed once per member.

#define LEXER_AMR_IO_D(X) X(Ui) X(Ua) X(Uo) X(Qi) X(Qo) X(Hi) X(Ho) \
    X(fsfin) X(fsfout) X(fsfinval) X(fsfoutval)
#define LEXER_AMR_RUN_D(X) X(phiin) X(phiout) X(pressgage) \
    X(wavetime) X(wts) X(wte) X(sedtime) X(sedsimtime) X(sedwavetime) \
    X(wT) X(wV) X(wH) X(wA) X(wL) X(wd) X(ww) X(wk) X(wC) X(wHs) X(wAs) X(wwp) X(ww_s) X(ww_e) X(wTp) X(wLp)
#define LEXER_AMR_RUN_I(X) X(wN)

#define LEXER_AMR_ONCE_D(X) X(maxlength) X(xcoormax) X(xcoormin) X(ycoormax) X(ycoormin) X(zcoormax) X(zcoormin) \
    X(RK_alpha) X(mindt) X(maxdt) X(deltax)
#define LEXER_AMR_ONCE_I(X) X(velcorr) X(count_statestart)

#define LEXER_AMR_BUILD_D(X) X(phimean) X(psi) X(psi0) X(dtsed) X(bedmax) X(bedmin) X(maxbed) X(minbed)
#define LEXER_AMR_BUILD_I(X) X(reini_iter)

#define LEXER_AMR_TIME_D(X) X(simtime) X(dt) X(dt_old)
#define LEXER_AMR_TIME_I(X) X(count)

namespace
{
// the contract check: 1 stop, 0 print once
int amr_check_mode()
{
    static const int mode = (getenv("REEF3D_AMR_CHECK")!=nullptr && atoi(getenv("REEF3D_AMR_CHECK"))>0) ? 1 : 0;
    return mode;
}

// equal as stored (NaN equal to NaN: not a change)
template<class T> bool same(T a, T b)
{
    return std::memcmp(&a,&b,sizeof(T))==0;
}
}

// the runtime state of a grid built from b (lexer(const lexer&, int))
void lexer::amr_inherit(const lexer &b)
{
#define C(m) m = b.m;
    LEXER_AMR_IO_D(C) LEXER_AMR_RUN_D(C) LEXER_AMR_RUN_I(C)
    LEXER_AMR_ONCE_D(C) LEXER_AMR_ONCE_I(C)
    LEXER_AMR_BUILD_D(C) LEXER_AMR_BUILD_I(C)
    LEXER_AMR_TIME_D(C) LEXER_AMR_TIME_I(C)
#undef C

    // OWN: the defaults of lexer_ini and control_calc
    umax=vmax=wmax=0.0;
    viscmax=kinmax=epsmax=pressmax=omegamax=0.0;
    pressmin=1.0e9;
    ufbmax=vfbmax=wfbmax=0.0;
    ufb=vfb=wfb=pfb=qfb=rfb=0.0;
    ufbi=vfbi=wfbi=pfbi=qfbi=rfbi=0.0;
    xg=yg=zg=xgn=ygn=zgn=0.0;
    phi_fb=theta_fb=psi_fb=0.0;
    fbdt=fbmax=sfdt=sfmax=0.0;
    volume1=volume2=volume3=0.0;
    tank_vol=0.0;           // the water volume of the NHFLOW log (driver_log, level 0)
    slidecells=0.0;
    maxtopovel=0.0;
    final_res=0.0;
    pcnorm=ucnorm=vcnorm=wcnorm=0.0;
    alpha=0.0;
    solver_status=solver_error=0;
    solveriter=preconiter=inneriter=sediter=0;
    uiter=viter=witer=kiniter=epsiter=poissoniter=laplaceiter=lsmiter=suspiter=topoiter=heatiter=concentrationiter=0;
    FSI_count=0;
    presstime=veltime=reinitime=turbtime=plstime=itertime=0.0;
    wavecalctime=meantime=totaltime=gcmeantime=gctotaltime=Xmeantime=Xtotaltime=0.0;
    susptime=gctime=xtime=field4time=0.0;
    utime=vtime=wtime=recontime=fsftime=dftime=kintime=epstime=0.0;
    poissontime=laplacetime=matrixtime=ptime=sftime=fbtime=fsitime=0.0;
    lsmtime=heattime=concentrationtime=printouttime=0.0;
    printtime=sedprinttime=fsfprinttime=fsfsedprinttime=probeprinttime=stateprinttime=exportprinttime=0.0;
    printcount=printcount_sixdof=0;

    amr_last.clear();
}

// before every step of a grid: the RUN group from level 0 (with the contract check)
void lexer::amr_follow(const lexer &b)
{
    const int mode = amr_check_mode();
    static bool told[64] = {false};
    int n = 0;
    auto report = [&](const char *name, const char *what)
    {
        if(mode==0 && n<64 && told[n])
        return;
        if(n<64)
        told[n] = true;
        std::cerr<<"REEF3D AMR contract: "<<name<<" "<<what<<" (lexer_amr.cpp)"<<std::endl;
        if(mode==1)
        MPI_Abort(MPI_COMM_WORLD,3);
    };

    // RUN: written by the grid's own objects since the last amr_follow?
    const bool have = !amr_last.empty();
    if(!have)
    amr_last.resize(64);
    n = 0;
#define R(m) { if(have && !same<double>((double)m,amr_last[n])) report(#m,"was changed by the objects of a grid (not a RUN member)"); \
               m = b.m; amr_last[n] = (double)m; ++n; }
    // IO: the block's own with open faces (amr_own_io), else RUN
#define IO(m) { if(amr_own_io==1) ++n; else R(m) }
    LEXER_AMR_IO_D(IO) LEXER_AMR_RUN_D(R) LEXER_AMR_RUN_I(R)
#undef IO
#undef R

    // ONCE: changed on level 0 after the grid was built?
#define O(m) { if(!same(m,b.m)) report(#m,"differs from level 0 (level 0 changed it during the run: a RUN member)"); ++n; }
    LEXER_AMR_ONCE_D(O) LEXER_AMR_ONCE_I(O)
#undef O
}
