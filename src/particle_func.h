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

#ifndef PARTICLEFUNC_H_
#define PARTICLEFUNC_H_

#include"increment.h"
#include"particles_obj.h"

class lexer;
class fdm;
class ghostcell;

class tracers_obj;
//class particles_obj;

/// Particle function class
/** A class containing all basic function to manipulate the position of tracers_objs. */
class particle_func: private increment
{
public:

protected:
    particle_func(lexer* p);
    virtual ~particle_func();
    
    // Parallelization
    int remove(lexer* p,tracers_obj* PP);
    int transfer(lexer* p,ghostcell* pgc, tracers_obj* PP,int);
    int transfer(lexer* p,ghostcell* pgc, particles_obj* PP,int);

    // Movement
    void advect(lexer* p, fdm* a, tracers_obj* PP,int=0,double=0,double=0,double=0);
    void advect(lexer* p, fdm* a, particles_obj* PP,int=0,double=0,double=0,double=0);
    void transport(lexer* p, fdm* a, particles_obj* PP,int=0);
    void make_stationary(lexer* p, fdm* a, tracers_obj* PP,int=0);
    void make_stationary(lexer* p, fdm* a, particles_obj* PP,int=0);
    void make_moving(lexer* p, fdm* a, particles_obj* PP);

    // Utility
    double reynolds(lexer* p, fdm* a, particles_obj* PP,int);
    double settling_vel(lexer* p, fdm* a, particles_obj* PP,int);
    double drag_coefficient(lexer* p, fdm* a, particles_obj* PP,int);
    double volume(particles_obj* PP,int);
    double maxParticlesPerCell(lexer* p, fdm* a, double,bool=true);
    int maxParticlesPerXY(lexer* p, fdm* a, double);
    void particlesPerCell(lexer* p, ghostcell* pgc, particles_obj* PP);
    void particleStressTensor(lexer* p, fdm* a, ghostcell* pgc, particles_obj* PP);
    void particleStressTensorUpdateIJK(lexer* p, fdm* a, particles_obj* PP);
    void updateParticleStressTensor(lexer* p, fdm* a, particles_obj* PP,int,int,int);
    double theta_s(lexer* p, fdm* a, particles_obj* PP,int,int,int);
    double drag_model(lexer* p, double,double,double,double,double) const;
    void debug(lexer* p, fdm* a, ghostcell* pgc, tracers_obj*);
    void fixPos(lexer* p, fdm* a, particles_obj* PP);

    // memory management
    void cleanup(lexer* p, fdm* a, tracers_obj* PP,int);
protected:
    /// @brief Inter-particle stresses per cell
    double* stressTensor;
    /// @brief Number of particles in a cell
    double* cellSum;
    /// @brief Volume change per column
    double* topoVolumeChange;

private:
    /// @brief Kinetic viscory of fluid\n Initialized using `lexer::W1` and `lexer::W2`
    const double kinVis;
    /// @brief Ratio of fluid and solid densities\n Initialized using `lexer::W1` and `lexer::S22`
    const double drho;
    /// @brief Constant for stress trensor calculation\n Initialized using `lexer::Q14` and given in Pascal
    const double Ps;
    /// @brief Constant for stress trensor calculation\n Initialized using `lexer::Q15` and should be in range of \f$2\leq\beta\leq5\f$
    const double beta;
    /// @brief Dampener for stress trensor calculation\n Used to dampen out sudden acceleration resulting from high packing densities, initialized using `lexer::Q16`
    const double epsilon;
    /// @brief Maximum solid volume fraction of a fully packed bed\n Usually between 60% and 65%
    const double theta_crit;
};

#endif