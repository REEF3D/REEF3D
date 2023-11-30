/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#define PC posmem[pcount] // ?
#define PARTLOOP for(n=maxparticle-1;n>=pcount;--n) // change to -> for(n=0;n<active;++n)

using namespace std;

#ifndef PARTICLE_F_H_
#define PARTICLE_F_H_

class particle_f : public particle_base, public norm_vec, public boundarycheck
{
public:
	particle_f(lexer*, fdm*, ghostcell*); // needed
	virtual ~particle_f(); // needed
	virtual void start(lexer*,fdm*,ghostcell*,ioflow*); // needed
    virtual void ini(lexer*,fdm*,ghostcell*,ioflow*); // needed
	virtual void setup(lexer*,fdm*,ghostcell*);  // ?
    
	void advect(lexer*,fdm*,ghostcell*,double**,int*,int);  // needed
    
    void seed_ini(lexer*,fdm*,ghostcell*);  // needed
	void seed(lexer*,fdm*,ghostcell*);  // needed
    void posseed(lexer*,fdm*,ghostcell*);  // needed -> rename or create sub functions for box, topo, etc
    void posseed_topo(lexer*,fdm*,ghostcell*);  // needed - To check
    
	void remove(lexer*,fdm*,ghostcell*);  // needed - at boundary
	void random_delete(lexer*,fdm*,ghostcell*);  // remove
	void parcount(lexer*,fdm*,ghostcell*); // ?
	void particlex(lexer*, fdm*, ghostcell*); // ?
	void xupdate(lexer*,fdm*,ghostcell*); // ?
    
    void allocate(lexer*,fdm*,ghostcell*); // needed
	void print_particles(lexer*,fdm*,ghostcell*); // needed
	void print_ascii(lexer*,fdm*,ghostcell*); // ? - no reference
	
	void setradius(lexer*,fdm*);  // ? - only ref setup
	void posradius(lexer*,fdm*,int); // ? - only ref setradius

	void normal(fdm*, double&,double&,double&,double&); // ? -  ref in posseed_topo
	void normreg(fdm*, int,int,int); // ? - only ref in normal
	
	field4 active; // ? - only use in particle_seed -> local declaration?
	field4 posnum; // ? - ref: count, delete, f, setup
	
	double **pos,**neg; // pos - needed // neg - remove
	double **posxs; // ? - allocate, particlex
	double **posxr; // ? - allocate, particlex
	int *pxs; // ? - allocate, particlex
	int *pxr; // ? - allocate, particlex
	int *posflag; // needed
	int *posmem; // ?
    
	int cellcount, ppcell; // cellcount vs ppcell?
    int pcount, pactive, partnum; // pcount vs pactive vs posactive vs partnum?
	int n,nn,q,qq,qn,count,check; // seed topo - make local
    double wx,wy,wz; // normal
    double di,dj,dk,dnorm; // normal
    double uvel,vvel,wvel; // remove
    int posactive,maxparticle; // ?
	int posactive_old, posbalance; // remove
	int gposactive_old, gposbalance; //remove
    int corrected,removed,xchange,reseeded; // remove: corrected, 
    double val1,val2; // normal - make local

    // old parameters
    double phix,phiy,phiz; // remove
    double xs,ys,zs; // remove
    double xp,yp,zp; // remove
    double xc,yc,zc; // remove
    double xcell,ycell,zcell; // remove
    double length,scalar,gamma; // remove
    double u1,u2,v1,v2,w1,w2; // make local
    double sp; // remove
    int ii,jj,kk; // needed

    int gnegactive,gposactive,gpcount,gncount,gcorrected,gremoved,greseeded,gxchange; // remove?
    

	double H,Hval,nvec[3],phival,lambda,cosinus; // remove H, Hval, cosinus - phival, lambda: seed_topo
	const double epsi,dx,dy,dz,rmin,rmax; // remove zero – rmin, rmax: seed_topo + radius + delete - epsi:?
	//int pnum;
	const int irand; // need
	const double drand; // need

	int maxpart, minpart; // remove
	
	
	// PRINTVTU
	void print_vtu(lexer*,fdm*,ghostcell*,double**,int*,int,int); // needed
	
	void pvtu_pos(fdm*,lexer*,ghostcell*); // needed
    void header_pos(fdm*,lexer*,ghostcell*); // needed
    void piecename_pos(fdm*,lexer*,ghostcell*, int); // needed
	
	char name[100],pname[100],epsvar[100]; // needed
    int iin,offset[100]; // needed
    float ffn; // needed print - maybe localise
    double printtime,printtime2; // printtime - unused -- remove printtime2
    int printcount; // needed - printing
};

#endif

