/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef SEDPART_MOVEMENT_H_
#define SEDPART_MOVEMENT_H_

#include"increment.h"
#include"slice4.h"

#include <stdio.h>
#include <fstream>

class lexer;
class fdm;
class ghostcell;
class field4a;
class sediment_fdm;
class turbulence;

class sediment_part;
class particles_obj;

enum seedReturn:int
{
    REMOVE=-1,STOP,CONTINUE
};


class partres : public increment
/// Model for the movement of sediment particles following Tavouktsoglou et al. (2021)
/// @author Alexander Hanke
/// @date 2024
{
public:
        partres(lexer *);
        ~partres();

        void setup(lexer *, fdm &, double &);
        seedReturn seeding(lexer *, particles_obj &, size_t &, double, bool=false);
        
        void move_RK2(lexer *, fdm &, ghostcell&, particles_obj &, sediment_fdm &, turbulence &);
        void move_RK2_step1(lexer *, fdm &, ghostcell&, particles_obj &, sediment_fdm &, turbulence &, int &, int &);
        void move_RK2_step2(lexer *, fdm &, ghostcell&, particles_obj &, sediment_fdm &, turbulence &, int &, int &);
        void move_RK3(lexer *, fdm &, ghostcell&, particles_obj &, sediment_fdm &, turbulence &);
        
        void advec_plain(lexer *, fdm &, particles_obj &, size_t, sediment_fdm &, turbulence&, 
                        double*, double*, double*, double*, double*, double*, 
                        double&, double&, double&, double);
        void advec_pic(lexer *, fdm &, particles_obj &, size_t, sediment_fdm &, turbulence&, 
                        double*, double*, double*, double*, double*, double*, 
                        double&, double&, double&, double);
                        
        void sandslide(lexer *, fdm &, ghostcell&, particles_obj &, sediment_fdm &, turbulence &);
        
        void transfer(lexer *, particles_obj &, size_t &);
        void remove(lexer *, particles_obj &, size_t &);
        void make_moving(lexer *, fdm &, particles_obj &);
        void erosion(lexer *, fdm &, particles_obj &, sediment_fdm &);
        void deposition(lexer *, fdm &, particles_obj &, sediment_fdm &);
        void update(lexer *, fdm &, ghostcell &, particles_obj &);
        void debug(lexer *, fdm &, ghostcell &, particles_obj &, sediment_fdm &);
        double volume(lexer *, fdm &, particles_obj &);
        void writeState(lexer *, ofstream &);
        void readState(lexer *, ifstream &);
        void setupState(lexer *, fdm &, ghostcell &, particles_obj &);
        void setParticleMax(double);
private:
        double maxParticlesPerCell(lexer *, fdm &, double,bool=true,bool=false);
        void particleStressTensor(lexer *, fdm &, ghostcell &, particles_obj &);
        void particleStressTensorUpdateIJK(lexer *, fdm &, particles_obj &);
        void updateParticleStressTensor(lexer *, fdm &, particles_obj &, int, int, int);
        double theta_s(lexer *, fdm &, particles_obj &, int, int, int) const;
        double drag_model(lexer *, double, double, double, double, double) const;
        double settling_velocity(lexer *, double, double, double, double, double) const;
        void particlePerCell(lexer *, ghostcell &, particles_obj &);
        void timestep(lexer *, ghostcell &, particles_obj &);
        int activateNew(lexer *, fdm &, particles_obj &);
        void relative_velocity(lexer *, fdm &, particles_obj &, size_t, double &, double &, double &);
        double drag_coefficient(double) const;
        void addParticleForTransfer(lexer *, particles_obj &, size_t , particles_obj [6], int &);
        
    // relax
    void relax_ini(lexer*);
    void relax(lexer*, ghostcell*, sediment_fdm*);
    double rf(lexer*, double, double);     
    double r1(lexer*, double, double);
    double distcalc(lexer*, double , double, double , double, double);
        
    double *tan_betaQ73,*betaQ73,*dist_Q73;
	double val;
            
            /// @brief Sum of particles belonging to the stationary bed
            double *cellSumTopo;
            /// @brief Sum of particles belonging to the mobile bed
            double *cellSum;
            /// @brief Stress tensor for the particle-particle interaction
            double *stressTensor;
            /// @brief Number of particles per 2D bed column
            double *columnSum;
            /// @brief Relative density of fluid and particle
            const double drho;
            /// @brief Inverse of kinetik viscosity of the fluid
            const double invKinVis;
            /// @brief Stress tensor parameter
            const double Ps;
            /// @brief Stress tensor parameter
            const double beta;
            /// @brief Stress tensor parameter
            const double epsilon;
            /// @brief Critical solid volume fraction
            const double theta_crit;

            double dx;
            slice4 bedChange;

            double velDist=1.6;
            
    double Fd, Fs, Ft, F_tot, Re_p;
    double F,G,H;
    double maxcount;
    double Umax,Uabs;
    double fac;
    
    const int irand;
	const double drand;
};

#endif