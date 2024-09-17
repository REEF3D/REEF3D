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

#include"part.h"
#include"lexer.h"
#include"ghostcell.h"

void part::xchange(lexer *p, ghostcell *pgc)
{
    xchange_count(p,pgc);
    pgc->gcpartnum(p,sendnum,recvnum);
    
    // check send/recv array size
    xchange_resize(p,pgc);
    
    // xchange part arrays
    xchange_fill(p,pgc,U);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,U);
    
    xchange_fill(p,pgc,V);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,V);
    
    xchange_fill(p,pgc,W);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,W);
    
    
    xchange_fill(p,pgc,URK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,URK1);
    
    xchange_fill(p,pgc,VRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,VRK1);
    
    xchange_fill(p,pgc,WRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,WRK1);
    
    
    xchange_fill(p,pgc,X);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,X);
    
    xchange_fill(p,pgc,Y);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Y);
    
    xchange_fill(p,pgc,Z);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Z);
    
    
    xchange_fill(p,pgc,XRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,XRK1);
    
    xchange_fill(p,pgc,YRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,YRK1);
    
    xchange_fill(p,pgc,ZRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,ZRK1);
    
    
    xchange_fill(p,pgc,Uf);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Uf);
    
    xchange_fill(p,pgc,Vf);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Vf);
    
    xchange_fill(p,pgc,Wf);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Wf);
    
    
    xchange_fill(p,pgc,D);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,D);
    
    xchange_fill(p,pgc,RO);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,RO);

}
