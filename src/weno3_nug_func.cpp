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

#include "weno3_nug_func.h"
#include "lexer.h"

weno3_nug_func::weno3_nug_func(lexer *p)
{
    ini(p);
}

void weno3_nug_func::ini(lexer* p)
{
    if(!iniflag)
    {
        qfx.resize(p->knox+2*marge);
        qfy.resize(p->knoy+2*marge);
        qfz.resize(p->knoz+2*marge);

        cfx.resize(p->knox+2*marge);
        cfy.resize(p->knoy+2*marge);
        cfz.resize(p->knoz+2*marge);

        isfx.resize(p->knox+2*marge);
        isfy.resize(p->knoy+2*marge);
        isfz.resize(p->knoz+2*marge);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag = true;
    }
}
