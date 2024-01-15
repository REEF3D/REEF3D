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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"particle.h"
#include"norm_vec.h"
#include"boundarycheck.h"
#include"field4.h"
#include "particles_obj.h"

#define PARTLOOP for(n=0;n<PP.loopindex;++n)

using namespace std;

#ifndef PARTICLE_F_H_
#define PARTICLE_F_H_

class particle_f : public particle_base, public norm_vec, public boundarycheck
{
public:
	particle_f(lexer*, fdm*, ghostcell*);
	virtual ~particle_f();
	virtual void start(lexer*,fdm*,ghostcell*,ioflow*);
    virtual void ini(lexer*,fdm*,ghostcell*,ioflow*);
    
	void advect(lexer*,fdm*,ghostcell*);
    
    void seed_ini(lexer*,fdm*,ghostcell*);
	void seed(lexer*,fdm*,ghostcell*);
    void posseed_box(lexer*,fdm*,ghostcell*);
    void posseed_topo(lexer*,fdm*,ghostcell*);
	void posseed_suspended(lexer*,fdm*,ghostcell*);
    
	void remove(lexer*,fdm*,ghostcell*);
	void particlex(lexer*, fdm*, ghostcell*);
    
    void allocate(lexer*,fdm*,ghostcell*);
	void print_particles(lexer*,fdm*,ghostcell*);
	void print_ascii(lexer*,fdm*,ghostcell*);
	
	field4 active_box;
	field4 active_topo;

	particles_obj PP;
	
	double **posxs;
	double **posxr;
	int *pxs;
	int *pxr;
    
	int n,i,j,k;
    int removed,xchange;
    
    int partnum, gpartnum;
    int ppcell;
	int gparticle_active;

    int gremoved,gxchange;
    
	const int irand;
	const double drand;
	double starttime;
	
	// PRINTVTU
	void print_vtu(lexer*,fdm*,ghostcell*);
	
	void pvtu_pos(fdm*,lexer*,ghostcell*);
    void header_pos(fdm*,lexer*,ghostcell*);
    void piecename_pos(fdm*,lexer*,ghostcell*, int);
	
	char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int gcval_phi;
    double printtime;
    int printcount;
};

#endif

