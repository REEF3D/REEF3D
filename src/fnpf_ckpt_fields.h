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
    fd("eta",'S',c->eta.data()); fd("eta_n",'S',c->eta_n.data()); fd("WL",'S',c->WL.data()); fd("bed",'S',c->bed.data());
    fd("depth",'S',c->depth.data()); fd("Fifsf",'S',c->Fifsf.data()); fd("Fz",'S',c->Fz.data()); fd("K",'S',c->K.data());
    fd("Fx",'S',c->Fx.data()); fd("Fy",'S',c->Fy.data()); fd("Ex",'S',c->Ex.data()); fd("Ey",'S',c->Ey.data());
    fd("Exu",'S',c->Exu.data()); fd("Eyu",'S',c->Eyu.data()); fd("Exx",'S',c->Exx.data()); fd("Eyy",'S',c->Eyy.data());
    fd("Bx",'S',c->Bx.data()); fd("By",'S',c->By.data()); fd("Bxx",'S',c->Bxx.data()); fd("Byy",'S',c->Byy.data());
    fd("coastline",'S',c->coastline.data()); fd("vb",'S',c->vb.data()); fd("Hs",'S',c->Hs.data()); fd("breaking_print",'S',c->breaking_print.data());
    fi("etaloc",'S',c->etaloc.data()); fi("wet_n",'S',c->wet_n.data()); fi("breaking",'S',c->breaking.data());
    fi("breaklog",'S',c->breaklog.data()); fi("bc",'S',c->bc.data());

    // potential and velocities
    fd("Fi",'F',c->Fi); fd("Uin",'F',c->Uin); fd("U",'F',c->U); fd("V",'F',c->V); fd("W",'F',c->W);

    // resolved bodies (X 10): body nodes and the Neumann data of the body faces (fnpf_6DOF)
    fd("FBF",'F',c->FBF); fd("FBu",'F',c->FBu); fd("FBv",'F',c->FBv); fd("FBw",'F',c->FBw);

    // the sigma grid and the flags of the lexer that change in time (wetting and drying)
    fd("sig",'F',p->sig); fd("sigx",'F',p->sigx); fd("sigy",'F',p->sigy); fd("sigxx",'F',p->sigxx);
    fd("sigz",'S',p->sigz); fd("ZSN",'F',p->ZSN); fd("ZSP",'C',p->ZSP);
    fi("flagslice4",'S',p->flagslice4); fi("flag7",'F',p->flag7); fi("flag4",'C',p->flag4.data());
    fi("wet",'S',p->wet); fi("lwet_n",'S',p->wet_n); fi("deep",'S',p->deep);
}

#endif
