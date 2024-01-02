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

class lexer;
class fdm;
class ghostcell;
class wave_lib;

using namespace std;

#ifndef WAVE_INTERFACE_H_
#define WAVE_INTERFACE_H_

class wave_interface 
{

public:
	wave_interface(lexer*, ghostcell*);
	virtual ~wave_interface();
	

	double wave_u(lexer*,ghostcell*,double,double,double);
	double wave_v(lexer*,ghostcell*,double,double,double);
    double wave_w(lexer*,ghostcell*,double,double,double);
    double wave_h(lexer*,ghostcell*,double,double,double);
    double wave_fi(lexer*,ghostcell*,double,double,double);
    double wave_eta(lexer*,ghostcell*,double,double);
    double wave_um(lexer*,ghostcell*,double,double);
    double wave_vm(lexer*,ghostcell*,double,double);
    
    
    // decomp
    double wave_u_space_sin(lexer*,ghostcell*,double,double,double,int);
    double wave_u_space_cos(lexer*,ghostcell*,double,double,double,int);
    double wave_u_time_sin(lexer*,ghostcell*,int);
    double wave_u_time_cos(lexer*,ghostcell*,int);
    
    double wave_v_space_sin(lexer*,ghostcell*,double,double,double,int);
    double wave_v_space_cos(lexer*,ghostcell*,double,double,double,int);
    double wave_v_time_sin(lexer*,ghostcell*,int);
    double wave_v_time_cos(lexer*,ghostcell*,int);
    
    double wave_w_space_sin(lexer*,ghostcell*,double,double,double,int);
    double wave_w_space_cos(lexer*,ghostcell*,double,double,double,int);
    double wave_w_time_sin(lexer*,ghostcell*,int);
    double wave_w_time_cos(lexer*,ghostcell*,int);
    
    double wave_eta_space_sin(lexer*,ghostcell*,double,double,int);
    double wave_eta_space_cos(lexer*,ghostcell*,double,double,int);
    double wave_eta_time_sin(lexer*,ghostcell*,int);
    double wave_eta_time_cos(lexer*,ghostcell*,int);
    
    double wave_fi_space_sin(lexer*,ghostcell*,double,double,double,int);
    double wave_fi_space_cos(lexer*,ghostcell*,double,double,double,int);
    double wave_fi_time_sin(lexer*,ghostcell*,int);
    double wave_fi_time_cos(lexer*,ghostcell*,int);
    
    void wave_prestep(lexer*,ghostcell*);


private:
    wave_lib *pwave;
    
    
    int n,m,count;
    int wtype;
    double wD;

	double starttime,endtime;
    static int printcheck;

};

#endif


