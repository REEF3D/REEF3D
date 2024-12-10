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

#include"sediment_void.h"

sediment_void::sediment_void()
{

}

sediment_void::~sediment_void()
{

}

void sediment_void::start_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow,
                                    reinitopo *preto, solver *psolv)
{
}

void sediment_void::ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
}

void sediment_void::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
}

void sediment_void::ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void sediment_void::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow*, slice &P, slice &Q)
{
}

void sediment_void::ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
}
    
void sediment_void::update_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo*)
{
}

void sediment_void::update_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
}

void sediment_void::update_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{
}

void sediment_void::relax(lexer *p,ghostcell *pgc)
{
}

double sediment_void::qbeval(int ii, int jj)
{
    double val=0.0;

    return val;
}

void sediment_void::qbeget(int ii, int jj, double val)
{
}

double sediment_void::bedzhval(int ii, int jj)
{
    double val=0.0;

    return val;
}

void sediment_void::print_2D_bedload(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::print_3D_bedload(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu_bedload(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtp_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::print_2D_bedshear(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::print_3D_bedshear(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtp_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::print_2D_parameter1(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::print_3D_parameter1(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtp_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::print_2D_parameter2(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::print_3D_parameter2(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtp_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

double sediment_void::bedshear_point(lexer *p, ghostcell *pgc)
{
	return 0.0;
}

void sediment_void::start_susp(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, solver *psolv)
{
}

void sediment_void::ctimesave(lexer *p, fdm* a)
{

}
