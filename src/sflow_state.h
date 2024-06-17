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

#ifndef SFLOW_STATE_H_
#define SFLOW_STATE_H_

#include"increment.h"
#include<fstream>

class lexer;
class fdm2D;
class ghostcell;

using namespace std;

class sflow_state : public increment
{

public:
	sflow_state(lexer*,fdm2D*,ghostcell*);
	virtual ~sflow_state();
	void write(lexer*,fdm2D*,ghostcell*);
    
    void ini_mainheader(lexer*,fdm2D*,ghostcell*);
    
    void write_result(lexer*,fdm2D*,ghostcell*);
    void write_mainheader(lexer*,fdm2D*,ghostcell*);
    void write_header(lexer*,fdm2D*,ghostcell*);
	
private:
    void filename_single(lexer*,fdm2D*,ghostcell*,int);
    void filename_continuous(lexer*,fdm2D*,ghostcell*);
    void filename_header(lexer*,fdm2D*,ghostcell*);
    
    char name[500];
    float ffn;
	int iin;
	double ddn;
	int printcount;
    int ini_token;
    int file_version,file_type;
    int qn;
    ofstream result;
    
    
    int is,ie,js,je;
    int is_global,ie_global,js_global,je_global;
    int is_global_root,ie_global_root,js_global_root,je_global_root;
    int is_flag,ie_flag,js_flag,je_flag;
    int flag;
    int *flag_all;
    int *is_flag_all,*ie_flag_all,*js_flag_all,*je_flag_all;
    int *is_global_all,*ie_global_all,*js_global_all,*je_global_all;
    
    
};

#endif
