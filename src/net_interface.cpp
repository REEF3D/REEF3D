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

#include"net_interface.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"

net_interface::net_interface(lexer *p, ghostcell *pgc)
{
    // no grid fields here: every 6DOF body (and sixdof_void) holds a net_interface, also without nets (X 320)
    // or membranes (X 330). The nets (initialize_cfd/_nhflow) and membranes (membrane_ini_nhflow) allocate
    // their own storage only when they are switched on.
}

net_interface::~net_interface()
{
}


void net_interface::ckpt_state(std::vector<double> &v) const
{
    v.push_back(double(pnet.size()));
    for(const net *n : pnet)
    {
        std::vector<double> s;
        n->ckpt_state(s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    }
}

bool net_interface::ckpt_restore(const std::vector<double> &v)
{
    if(v.empty() || size_t(v[0])!=pnet.size())
    return false;
    size_t k = 1;
    for(net *n : pnet)
    {
        if(k>=v.size())
        return false;
        const size_t len = size_t(v[k++]);
        if(k+len>v.size())
        return false;
        const std::vector<double> s(v.begin()+k,v.begin()+k+len);
        size_t j = 0;
        n->ckpt_restore(s,j);
        if(j!=len)
        return false;
        k += len;
    }
    return k==v.size();
}
