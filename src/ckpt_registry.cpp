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

#include"ckpt_registry.h"

std::map<unsigned long long,ckpt_registry::entry>& ckpt_registry::table()
{
    static std::map<unsigned long long,entry> t;
    return t;
}

std::mutex& ckpt_registry::lock()
{
    static std::mutex m;
    return m;
}

unsigned long long ckpt_registry::add(const void *owner, char kind, bool isint, void *V, bool local)
{
    static unsigned long long next = 0;
    std::lock_guard<std::mutex> g(lock());
    const unsigned long long id = ++next;
    table()[id] = entry{owner,kind,isint,V,local};
    return id;
}

unsigned long long ckpt_registry::add_raw(const void *owner, double *V, long long n, bool local)
{
    const unsigned long long id = add(owner,'R',false,V,local);
    std::lock_guard<std::mutex> g(lock());
    table()[id].n = n;
    return id;
}

void ckpt_registry::remove(unsigned long long id)
{
    if(id==0)
    return;
    std::lock_guard<std::mutex> g(lock());
    table().erase(id);
}

void ckpt_registry::visit(const void *owner, const std::function<void(const entry&)> &f)
{
    for(const auto &e : table())
    if(e.second.owner==owner)
    f(e.second);
}
