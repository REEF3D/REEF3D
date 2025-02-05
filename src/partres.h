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
Authora: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef PARTRES2_H_
#define PARTRES2_H_

#include"increment.h"
#include"part.h"
#include"slice4.h"
#include"field4a.h"
#include"fieldint4a.h"
#include"boundarycheck.h"

class lexer;
class fdm;
class ghostcell;
class sediment_fdm;
class turbulence;
class part;
class vrans;

using namespace std;


class partres : public increment
{
public:
    partres(lexer*, ghostcell*);
    ~partres();
    
    void move_RK2(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*);
    
    void advec_plain(lexer*, fdm*, part&, sediment_fdm*, turbulence*, 
                        double*, double*, double*, double*, double*, double*, 
                        double&, double&, double&, double);
    void advec_mppic(lexer*, fdm*, part&, sediment_fdm*, turbulence*, 
                        double*, double*, double*, double*, double*, double*, 
                        double&, double&, double&, double);
    
    // drag
    double drag_model(lexer *, double, double, double, double);
    double drag_coefficient(double);
    
    void update(lexer*, fdm*, ghostcell*, sediment_fdm*, field&, field&);
    void count_particles(lexer*, fdm*, ghostcell*, sediment_fdm*);
    
    void stress_tensor(lexer*, ghostcell*, sediment_fdm*);
    void stress_gradient(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void pressure_gradient(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void cellSum_update(lexer*, ghostcell*, sediment_fdm*,int);
    void cellSum_full_update(lexer*, ghostcell*, sediment_fdm*,int);
    
    void bedchange(lexer*, fdm*, ghostcell*, sediment_fdm*,int);
    void bedchange_update(lexer*, ghostcell*, sediment_fdm*,int);
    
    void timestep(lexer*, ghostcell*);
    
    void seed_topo(lexer*, fdm*, ghostcell*, sediment_fdm*);
    
    void boundcheck(lexer*, fdm*, ghostcell*, sediment_fdm*, int);
    
    part P;
    
    slice4 bedch;
    
    field4a Tau,Ts;
    field4a cellSum;
    
    
    // relax
    void relax_ini(lexer*);
    void relax(lexer*, ghostcell*, sediment_fdm*);
    double rf(lexer*, double, double);     
    double r1(lexer*, double, double);
    double distcalc(lexer*, double , double, double , double, double);
    
    // print
    void print_particles(lexer*,sediment_fdm*);
    void print_vtp(lexer*,sediment_fdm *);
    void pvtp(lexer*);
    void piecename_pos(lexer*, int);
    void header_pos(lexer*);
        
        
private:
    
    boundarycheck boundaries;
    
    const int irand;
	const double drand;
    
    int num, printcount;
    double printtime;
    char name[100];
    char pname[100];
    
    int particle_count,empty_count,active_count;
    
    
    double topoval,solidval;
    
    double F,G,H;
    
    field4a dPx,dPy,dPz;
    field4a dTx,dTy,dTz;
    
    double dPx_val,dPy_val,dPz_val;
    double dTx_val,dTy_val,dTz_val;
    
    double Bx,By,Bz;
    double Urel,Vrel,Wrel;
    double uf,vf,wf;
    double Dpx,Dpy,Dpz;
    double Uabs_rel;
    double Fs,Fd,Ftot;
    
    double Tsval,T,Dp,Tf;
    double velDist;
    double Umax;
    double Sx,Sy,Sz;
    
    double Rep,Cd;
    
    double *tan_betaQ73,*betaQ73,*dist_Q73;
    
    double Ps,beta,epsilon,Tc;
    
    int timestep_ini;
    
    double DragCoeff;
    double F_tot;
};

#endif