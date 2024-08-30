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

#ifndef SEDPART_H_
#define SEDPART_H_

#include "sediment.h"
#include "particle_func.h"
#include "increment.h"

#include "particles_obj.h"
#include "field4.h"
#include "bedshear_reduction.h"
#include "sediment_fdm.h"
#include "bedshear.h"
#include "bedslope.h"
#include "topo_relax.h"

namespace sediment_particle::movement
{
class base;
}

class lexer;
class fdm;
class ghostcell;
class ioflow;
class solver;
class reinitopo;
class fdm2D;
class slice;
class ofstrem;
class vrans;
class turbulence;

/// This class used particles on a Lagrangien framework and a VRANS sediment domain to simulate the influence of flow on the sediment
class sedpart : public sediment, private particle_func, private increment
{
public:

    sedpart(lexer* p,ghostcell* pgc ,turbulence* pturb);
    virtual ~sedpart();

    // CFD methods

    void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*) override;
    void ini_cfd(lexer*,fdm*,ghostcell*) override;
    void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*) override;

    // SFLOW methods
    
    void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&) override;
    void ini_sflow(lexer*, fdm2D*, ghostcell*) override;
    void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*) override;
    
    // 

    void erode(lexer*,fdm*);
    void deposit(lexer*,fdm*);

    // Statefile methods

    void write_state_particles(lexer *, ofstream&);
    void read_state_particles(lexer *, ifstream&);

    // Printing methods bedshear

    virtual void print_3D_bedshear(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedshear(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);

protected:

private:

    // 

    void fill_PQ_cfd(lexer*,fdm*,ghostcell*);

    // Seeding methods

    void seed_ini(lexer*,fdm*,ghostcell*);
    void seed(lexer*,fdm*);
    void posseed_box(lexer*,fdm*);
    void posseed_topo(lexer*,fdm*);
    void posseed_suspended(lexer*,fdm*);
    void point_source(lexer*,fdm*);
    void topo_influx(lexer*,fdm*);
    void seed_srand(lexer*);
    void seed_topo(lexer*,fdm*);
    void solid_influx(lexer*,fdm*);
    size_t set_active_topo(lexer*, fdm*);

    void debug(lexer*,fdm*,ghostcell*);

    // Printing methods
    
	void print_particles(lexer*);
	void print_vtp(lexer*);
	void pvtp_pos(lexer*);
    void header_pos(lexer*);
    void piecename_pos(lexer*,int);

public:

protected:

private:

    /// @brief Current capacity for particles
    int maxparticle;
    /// @brief Desired particles per cell
    int ppcell;
    /// @brief Number of particles over all partitions
    int gparticle_active;
    /// @brief Particles removed over all partitions
    int gremoved;
    /// @brief Particles exchanged between all partitions
    int gxchange;

    /// @brief Total volume of bed in cubic meter
    double volume;
    /// @brief Initial total volume of bed in cubic meter
    double volume0;

    /// @brief integer rand() scaler
    const int irand;
	/// @brief double rand() normalizer
	const double drand;

    /// @brief Initialization counter
    int inicount;

    /// @brief Marker for cells which should be seeded for a box
    field4 active_box;
    /// @brief Marker for cells which should be seeded with topography
	field4 active_topo;

    // SEDIMENT OBJECTS

    /// @brief Particle object
    particles_obj PP;
    /// @brief VRANS object
    vrans* pvrans;
    /// @brief FDM object for sediment
    sediment_fdm s;
    /// @brief Bed shear object
    bedshear pbedshear;
    /// @brief Movement object for particles
    sediment_particle::movement::base *movement;

    turbulence *pturb;

    topo_relax prelax;
    bedslope pslope;
    bedshear_reduction* preduce;

	// PRINT

	/// @brief Output file name
	char name[100];
    /// @brief Name of individual output files
    char pname[100];
    /// @brief Latest printed time
    double printtime;
    /// @brief Number of print iterations
    int printcount;
    /// @brief File numer
    int num;

    /// @brief Print out precision
    int prec;

    // DEFINITIONS

    #define PARTLOOP for(size_t n=0;n<PP.loopindex;n++)
    #define PARTICLELOOP for(size_t n=0;n<PP.loopindex;n++) if(PP.Flag[n]>INT32_MIN)
};

#endif