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

#ifndef FNPF_CKPT_FIELDS_H_
#define FNPF_CKPT_FIELDS_H_

// the state of an FNPF grid at the end of a step (checkpoints, P 48, I 45): level 0 and every patch
// of the mesh refinement.  f(name, kind, pointer, may be absent): kind S slice, F node field (FIJK),
// C cell field (IJK); ints and doubles.

#include"lexer.h"
#include"fdm_fnpf.h"

template<class FD, class FI>
void fnpf_ckpt_fields(lexer *p, fdm_fnpf *c, FD &&fd, FI &&fi)
{
    // slices of fdm_fnpf
    fd("eta",'S',c->eta.V); fd("eta_n",'S',c->eta_n.V); fd("WL",'S',c->WL.V); fd("bed",'S',c->bed.V);
    fd("depth",'S',c->depth.V); fd("Fifsf",'S',c->Fifsf.V); fd("Fz",'S',c->Fz.V); fd("K",'S',c->K.V);
    fd("Fx",'S',c->Fx.V); fd("Fy",'S',c->Fy.V); fd("Ex",'S',c->Ex.V); fd("Ey",'S',c->Ey.V);
    fd("Exu",'S',c->Exu.V); fd("Eyu",'S',c->Eyu.V); fd("Exx",'S',c->Exx.V); fd("Eyy",'S',c->Eyy.V);
    fd("Bx",'S',c->Bx.V); fd("By",'S',c->By.V); fd("Bxx",'S',c->Bxx.V); fd("Byy",'S',c->Byy.V);
    fd("coastline",'S',c->coastline.V); fd("vb",'S',c->vb.V); fd("Hs",'S',c->Hs.V); fd("breaking_print",'S',c->breaking_print.V);
    fi("etaloc",'S',c->etaloc.V); fi("wet_n",'S',c->wet_n.V); fi("breaking",'S',c->breaking.V);
    fi("breaklog",'S',c->breaklog.V); fi("bc",'S',c->bc.V);

    // potential and velocities
    fd("Fi",'F',c->Fi); fd("Uin",'F',c->Uin); fd("U",'F',c->U); fd("V",'F',c->V); fd("W",'F',c->W);

    // resolved bodies (X 10): body nodes and the Neumann data of the body faces (fnpf_6DOF)
    fd("FBF",'F',c->FBF); fd("FBu",'F',c->FBu); fd("FBv",'F',c->FBv); fd("FBw",'F',c->FBw);

    // the sigma grid and the flags of the lexer that change in time (wetting and drying)
    fd("sig",'F',p->sig); fd("sigx",'F',p->sigx); fd("sigy",'F',p->sigy); fd("sigxx",'F',p->sigxx);
    fd("sigz",'S',p->sigz); fd("ZSN",'F',p->ZSN); fd("ZSP",'C',p->ZSP);
    fi("flagslice4",'S',p->flagslice4); fi("flag7",'F',p->flag7); fi("flag4",'C',p->flag4);
    fi("wet",'S',p->wet); fi("lwet_n",'S',p->wet_n); fi("deep",'S',p->deep);
}

#endif
