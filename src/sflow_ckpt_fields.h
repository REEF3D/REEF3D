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

#ifndef SFLOW_CKPT_FIELDS_H_
#define SFLOW_CKPT_FIELDS_H_

// the state of an SFLOW grid at the end of a step (checkpoints, P 48, I 45): level 0 and every patch
// of the mesh refinement.  fd/fi(name, kind, pointer), kind S: slices (all slices of fdm2D, they are
// small) and the wet/dry flags of the lexer.  The momentum, pressure and turbulence objects add their
// own arrays (ckpt_fields).

#include"lexer.h"
#include"fdm2D.h"

template<class FD, class FI>
void sflow_ckpt_fields(lexer *p, fdm2D *b, FD &&fd, FI &&fi)
{
    fd("U",'S',b->U.data()); fd("V",'S',b->V.data()); fd("W",'S',b->W.data()); fd("UH",'S',b->UH.data()); fd("VH",'S',b->VH.data());
    fd("WH",'S',b->WH.data()); fd("UA",'S',b->UA.data()); fd("VA",'S',b->VA.data()); fd("MX",'S',b->MX.data());
    fd("MY",'S',b->MY.data()); fd("Fx",'S',b->Fx.data()); fd("Fs",'S',b->Fs.data()); fd("Fn",'S',b->Fn.data());
    fd("FEx",'S',b->FEx.data()); fd("Fy",'S',b->Fy.data()); fd("Fe",'S',b->Fe.data()); fd("Fw",'S',b->Fw.data());
    fd("FEy",'S',b->FEy.data()); fd("Ss",'S',b->Ss.data()); fd("Sn",'S',b->Sn.data()); fd("SSx",'S',b->SSx.data());
    fd("Se",'S',b->Se.data()); fd("Sw",'S',b->Sw.data()); fd("SSy",'S',b->SSy.data()); fd("Un",'S',b->Un.data());
    fd("Us",'S',b->Us.data()); fd("Ue",'S',b->Ue.data()); fd("Uw",'S',b->Uw.data()); fd("Vn",'S',b->Vn.data());
    fd("Vs",'S',b->Vs.data()); fd("Ve",'S',b->Ve.data()); fd("Vw",'S',b->Vw.data()); fd("Wn",'S',b->Wn.data());
    fd("Ws",'S',b->Ws.data()); fd("We",'S',b->We.data()); fd("Ww",'S',b->Ww.data()); fd("UHn",'S',b->UHn.data());
    fd("UHs",'S',b->UHs.data()); fd("UHe",'S',b->UHe.data()); fd("UHw",'S',b->UHw.data()); fd("VHn",'S',b->VHn.data());
    fd("VHs",'S',b->VHs.data()); fd("VHe",'S',b->VHe.data()); fd("VHw",'S',b->VHw.data()); fd("WHn",'S',b->WHn.data());
    fd("WHs",'S',b->WHs.data()); fd("WHe",'S',b->WHe.data()); fd("WHw",'S',b->WHw.data()); fd("QUn",'S',b->QUn.data());
    fd("QUs",'S',b->QUs.data()); fd("QUe",'S',b->QUe.data()); fd("QUw",'S',b->QUw.data()); fd("QVn",'S',b->QVn.data());
    fd("QVs",'S',b->QVs.data()); fd("QVe",'S',b->QVe.data()); fd("QVw",'S',b->QVw.data()); fd("ETAs",'S',b->ETAs.data());
    fd("ETAn",'S',b->ETAn.data()); fd("ETAe",'S',b->ETAe.data()); fd("ETAw",'S',b->ETAw.data()); fd("Ds",'S',b->Ds.data());
    fd("Dn",'S',b->Dn.data()); fd("De",'S',b->De.data()); fd("Dw",'S',b->Dw.data()); fd("dfx",'S',b->dfx.data());
    fd("dfy",'S',b->dfy.data()); fd("WL",'S',b->WL.data()); fd("eta",'S',b->eta.data()); fd("eta_n",'S',b->eta_n.data());
    fd("P",'S',b->P.data()); fd("Pn",'S',b->Pn.data()); fd("F",'S',b->F.data()); fd("Q",'S',b->Q.data()); fd("Qn",'S',b->Qn.data());
    fd("G",'S',b->G.data()); fd("H",'S',b->H.data()); fd("L",'S',b->L.data()); fd("Fext",'S',b->Fext.data());
    fd("Gext",'S',b->Gext.data()); fd("Hext",'S',b->Hext.data()); fd("ws",'S',b->ws.data()); fd("press",'S',b->press.data());
    fd("eddyv",'S',b->eddyv.data()); fd("kin",'S',b->kin.data()); fd("eps",'S',b->eps.data()); fd("hx",'S',b->hx.data());
    fd("hy",'S',b->hy.data()); fd("hp",'S',b->hp.data()); fd("dpx",'S',b->dpx.data()); fd("dpy",'S',b->dpy.data());
    fd("test",'S',b->test.data()); fd("Hs",'S',b->Hs.data()); fd("fs",'S',b->fs.data()); fd("fx",'S',b->fx.data());
    fd("fy",'S',b->fy.data()); fd("fz",'S',b->fz.data()); fd("bed",'S',b->bed.data()); fd("bed0",'S',b->bed0.data());
    fd("depth",'S',b->depth.data()); fd("solidbed",'S',b->solidbed.data()); fd("topobed",'S',b->topobed.data());
    fd("bednode",'S',b->bednode.data()); fd("breaking_print",'S',b->breaking_print.data()); fd("ks",'S',b->ks.data());
    fi("breaking",'S',b->breaking.data()); fi("wet1",'S',b->wet1.data());
    fi("deep1",'S',b->deep1.data()); fi("wet2",'S',b->wet2.data()); fi("deep2",'S',b->deep2.data());

    // the wet/dry flags of the lexer
    fi("wet",'S',p->wet); fi("lwet_n",'S',p->wet_n); fi("deep",'S',p->deep);
}

#endif
