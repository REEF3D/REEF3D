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

#ifndef WAVE_LIB_FNPF_COUPLE_H_
#define WAVE_LIB_FNPF_COUPLE_H_

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include<mpi.h>
#include<vector>
#include<deque>

class fdm;
class sixdof;

using namespace std;

// B 92 62: waves from an FNPF model running in the same MPI job (reefcouple), receiver side
class wave_lib_fnpf_couple final : public wave_lib_precalc, public wave_lib_parameters
{
public:
    wave_lib_fnpf_couple(lexer*, ghostcell*);
    virtual ~wave_lib_fnpf_couple();
    
    double wave_u(lexer*,double,double,double) override final;
    double wave_v(lexer*,double,double,double) override final;
    double wave_w(lexer*,double,double,double) override final;
    double wave_eta(lexer*,double,double) override final;
    double wave_fi(lexer*,double,double,double) override final;
    
    void parameters(lexer*,ghostcell*) override final;
    void wave_prestep(lexer*,ghostcell*) override final;
    
    // cached evaluation in the relaxation zones: the horizontal stencil once
    void wave_cache_points(lexer*, const std::vector<double> &x, const std::vector<double> &y) override final;
    double wave_eta_c(lexer*, int) override final;
    double wave_fi_c(lexer*, int, double) override final;
    double wave_u_c(lexer*, int, double) override final;
    double wave_v_c(lexer*, int, double) override final;
    double wave_w_c(lexer*, int, double) override final;
    void wave_uvw_c(lexer*, int, double, double&, double&, double&) override final;
    
    // end of the receiver: stop to the source, drain its messages
    void finish();
    
    // feedback to the source at the start of the step (I 243 bodies, I 244 overlap box with the CFD fdm)
    void feedback(lexer*, ghostcell*, sixdof*, fdm*);
    
private:
    // horizontal stencil in the window: up to 4 columns with weights
    struct hstencil
    {
        int col[4];
        double w[4];
    };
    
    struct source
    {
        int rank;
        int i0,i1,j0,j1,ncol;
        bool ended=false;
        std::vector<double> buf;
    };
    
    struct level
    {
        double t;
        std::vector<double> E,U,V,W,F;
    };
    
    // overlap column (I 244): global FNPF indices, weight, sampling point in the receiver frame
    struct ocol
    {
        int gi,gj;
        double r,x,y;
        bool own;
        // FNPF 2D -> receiver 3D: the column stands for a crest line, averaged over the box
        std::vector<double> sx,sy;
        std::vector<char> sown;
    };
    
    void handshake(lexer*, ghostcell*);
    bool recv_level();
    void weights(double);
    void sync(lexer*);
    
    hstencil locate(double, double) const;
    void frame(double, double, double&, double&) const;
    double eta_h(const hstencil&) const;
    void vel_h(const hstencil&, double, double&, double&, double&) const;
    double fi_h(const hstencil&, double) const;
    void overlap_ini(lexer*);
    bool surface_cfd(lexer*, fdm*, double, double, double&, double&, double&, double&) const;
    void rotate(double, double, double&, double&) const;
    
    MPI_Comm ic;
    lexer *p;
    
    // FNPF grid
    int gx,gy,nz,jdirF,nfnpf;
    double wdF,t0F;
    std::vector<double> X,Y,ZN;
    
    // window of this rank in FNPF indices [wi0,wi1) x [wj0,wj1)
    int wi0,wi1,wj0,wj1,nwx,nwy,nwc;
    std::vector<double> B;        // bed of the window columns
    
    std::vector<source> S;
    std::deque<level> lev;
    bool ended=false;
    bool finished=false;
    
    // time interpolation: level index and weight
    std::vector<int> wl;
    std::vector<double> ww_;
    int order;
    double tlast;
    
    // frame: B 105 angle and origin, offsets I 231-I 233
    double cg,sg,dx0,dy0,dz0;
    
    std::vector<hstencil> cache;
    
    // Fi in the data (FNPF receiver)
    bool fi = false;
    
    // feedback
    bool fb_on = false;
    int nb_fb = 0;
    int fbsize = 0;
    std::vector<ocol> oc;
    std::deque<std::vector<double>> fbbuf;
    std::deque<MPI_Request> fbreq;
    
    int nwarn=0;
};

#endif
