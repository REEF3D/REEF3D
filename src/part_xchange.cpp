/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

void part::xchange(lexer *p, ghostcell *pgc, slice &bedch, int mode)
{
    // count send/recv
    xchange_count(p,pgc,mode);
    pgc->gcpartnum(p,sendnum,recvnum);
    
    // check send/recv array size
    xchange_resize(p,pgc);
    
    // sendid
    xchange_sendid(p,pgc,mode);
    
    // xchange part arrays
    xchange_fill(p,pgc,mode,U);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,U);
    
    xchange_fill(p,pgc,mode,V);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,V);
    
    xchange_fill(p,pgc,mode,W);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,W);
    
    
    xchange_fill(p,pgc,mode,URK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,URK1);
    
    xchange_fill(p,pgc,mode,VRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,VRK1);
    
    xchange_fill(p,pgc,mode,WRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,WRK1);
    
    
    xchange_fill(p,pgc,mode,X);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,X);
    
    xchange_fill(p,pgc,mode,Y);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Y);
    
    xchange_fill(p,pgc,mode,Z);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Z);
    
    
    xchange_fill(p,pgc,mode,XRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,XRK1);
    
    xchange_fill(p,pgc,mode,YRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,YRK1);
    
    xchange_fill(p,pgc,mode,ZRK1);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,ZRK1);
    
    
    xchange_fill(p,pgc,mode,Uf);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Uf);
    
    xchange_fill(p,pgc,mode,Vf);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Vf);
    
    xchange_fill(p,pgc,mode,Wf);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Wf);
    
    
    xchange_fill(p,pgc,mode,D);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,D);
    
    xchange_fill(p,pgc,mode,RO);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,RO);
    
    xchange_fill(p,pgc,mode,Test);
    pgc->gcpartx(p,sendnum,recvnum,send,recv);
    xchange_fillback(p,pgc,Test);
    
    // Flag
    xchange_fillback_flag(p,pgc,bedch,mode);
    xchange_fill_flag(p,pgc,mode);
    
}
