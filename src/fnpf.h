/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

class lexer;
class fdm_fnpf;
class ghostcell;
class solver;
class convection;
class ioflow;
class reini;
class onephase;

using namespace std;

#ifndef FNPF_H_
#define FNPF_H_

class fnpf
{
public:
	virtual void start(lexer*, fdm_fnpf*, ghostcell*, solver*, convection*, ioflow*, reini*,onephase*)=0;
    virtual void ini(lexer*, fdm_fnpf*, ghostcell*, ioflow*, reini*, onephase*)=0;
    virtual void inidisc(lexer*, fdm_fnpf*, ghostcell*, ioflow*, solver*)=0;
    virtual void ini_wetdry(lexer*, fdm_fnpf*, ghostcell*)=0;

};

#endif
