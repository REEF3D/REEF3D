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

#include"nodefill.h"
#include"fieldint5.h"
#include"field5.h"

class lexer;
class fdm;
class ghostcell;

#ifndef FSF_VTP_H_
#define FSF_VTP_H_

using namespace std;

class fsf_vtp :  public nodefill
{

public:
	fsf_vtp(lexer*,fdm*,ghostcell*);
	virtual ~fsf_vtp();
	virtual void start(lexer*,fdm*,ghostcell*);

private:
	void triangulation(lexer*, fdm*, ghostcell*, field&);
	void reconstruct(lexer*, fdm*, field&);
	void addpoint(lexer*,fdm*,int,int);
	void finalize(lexer*,fdm*);
	
	int **tri, **facet, *confac, *numfac,*numpt;
	double **ccpt, **pt, *ls;
	double   dV1,dV2,C1,C2,mi;
	int numtri,numvert, numtri_mem, numvert_mem;
	int count,countM,n,nn,q;
	int ccptcount,facount,check;
	int polygon_sum,polygon_num,vertice_num;
	const double zero,interfac;
    double epsi;
	
	fieldint5 vertice, nodeflag;
	field5 eta;
	
	void print(lexer*,fdm*,ghostcell*);
    void pvtp(lexer*,fdm*,ghostcell*);
    void header(lexer*,fdm*,ghostcell*);
    void name_iter(lexer*,fdm*,ghostcell*);
    void name_time(lexer*,fdm*,ghostcell*);
    void piecename(lexer*,fdm*,ghostcell*, int);

    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int gcval_phi;
    double printtime,printtime2;
	int fsfprintcount;
	

};

#endif


