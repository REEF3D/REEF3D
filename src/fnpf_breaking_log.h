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

#include"increment.h"
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;

using namespace std;

#ifndef FNPF_BREAKING_LOG_H_
#define FNPF_BREAKING_LOG_H_

class fnpf_breaking_log : public increment
{

public:
	fnpf_breaking_log(lexer*,fdm_fnpf*,ghostcell*);
	virtual ~fnpf_breaking_log();
	void write(lexer*,fdm_fnpf*,ghostcell*);

	
private:
    void filename(lexer*,fdm_fnpf*,ghostcell*);

    char name[500];
    float ffn;
	int iin;
	double ddn;
	int printcount;
    int ini_token;
    int file_version;
    ofstream result;
};

#endif
