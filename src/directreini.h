/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"reini.h"
#include"gradient.h"
#include"resize.h"
#include"fieldint4.h"
#include"fieldint5.h"
#include"field4.h"

class picard;
class fieldint;
class vec;

using namespace std;

#ifndef DIRECTREINI_H_
#define DIRECTREINI_H_

class directreini : public reini, gradient, public resize_class
{
public:
	directreini(lexer* p, fdm *a);
	virtual ~directreini();
	virtual void start(fdm*,lexer*,field&, ghostcell*,ioflow*);
    virtual void startV(fdm*,lexer*,vec&,ghostcell*,ioflow*);
	
	virtual void vtp(lexer*,fdm*,ghostcell*);
	virtual void name_iter(lexer*,fdm*,ghostcell*);
	virtual void pvtp(lexer*,fdm*,ghostcell*);
	virtual void piecename(lexer*,fdm*,ghostcell*,int);

	double dstx, dsty, dstz, dnorm, sign;
	double sx,sy,sz,snorm,op;
	
	char name[200],pname[200],epsvar[200];
    int iin,offset[200];
    float ffn;

private:
    fieldint5 vertice, nodeflag;
	
	field4 d0;
	fieldint4 wallf;

    picard *ppicard;
    reini *ppreini;

	void triangulation(lexer*, fdm*, field&, fieldint&, fieldint&);
	void reconstruct(lexer*, fdm*, field&, fieldint&, fieldint&);
	void reini(lexer*, fdm*, ghostcell*, field&, fieldint&, fieldint&);
	void finalize(lexer *p, fdm*);
	void debug(lexer*,fdm*);
	
	void constraint(lexer *p, fdm*, ghostcell*, field&);
	void correction(lexer *p, fdm*, ghostcell*, field&);

	
	double determinant(double,double,double,double,double,double,double,double,double,
                        double,double,double,double,double,double,double);


	double sg(double);

	void addpoint(lexer*,fdm*,int,int);

	double starttime,endtime;

	int **tri, **facet, *confac, *numfac, **ijk, *reiniflag;
	double **ccpt, **pt, *ls, *lsfac,*lsvert;
	double  *ls1,*ls0, dV1,dV2,C1,C2,mi,eta;
	int numtri,numvert, numtri_mem, numvert_mem;
	int count,countM,n,nn,q;

	int ccptcount,facount,check;
	
	int polygon_sum,polygon_num,vertice_num;

	int gcval_phi,gcval_ro,gcval_iniphi,reiniter;
	const double epsi,zero;
	
	double H,H0,grad,dT,dirac;
	double lambda1,lambda2,dV,dval,Cs;
	void wallf_update(lexer*,fdm*,ghostcell*);
	
	double dx,dy,dz;
};

#endif

