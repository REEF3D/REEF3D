/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef FNPF_FSF_H_
#define FNPF_FSF_H_

class lexer;
class ghostcell;
class fdm_fnpf;
class fnpf_laplace;
class slice;
class fnpf_convection;
class fnpf_ddx;
class fnpf_etadisc;

using namespace std;

class ioflow;

class fnpf_fsf 
{
public:
    // deleted through fnpf_fsf* (fnpf_amr): without a virtual destructor the derived destructors (fnpf_breaking: the
    // crest tracking) did not run
    virtual ~fnpf_fsf() {}
   
    virtual void fsfdisc(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&)=0;
    virtual void fsfdisc_ini(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&)=0;
    virtual void kfsfbc(lexer*,fdm_fnpf*,ghostcell*)=0;
    virtual void dfsfbc(lexer*,fdm_fnpf*,ghostcell*,slice&)=0;
    virtual void fsfwvel(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&)=0;
    virtual void wetdry(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&)=0;
    virtual void breaking(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&,slice&,double)=0;
    // the relaxation zones for the breaking models (A 350 4), once
    virtual void breaking_zones(lexer*,ioflow*) {}
    virtual void coastline_eta(lexer*,fdm_fnpf*,ghostcell*,slice&)=0;
    virtual void coastline_fi(lexer*,fdm_fnpf*,ghostcell*,slice&)=0;
    virtual void coastline_fi_ini(lexer*,fdm_fnpf*,ghostcell*,slice&)=0;
    virtual void coastline_vel(lexer*,fdm_fnpf*,ghostcell*,double*)=0;
    // the same on the surface velocities (no 3D U, V, W, fnpf_memory.h)
    virtual void coastline_vel_sl(lexer*,fdm_fnpf*,ghostcell*,slice&) {}
    virtual void damping(lexer*,fdm_fnpf*,ghostcell*,slice&,int,double)=0;
    
};

#endif
