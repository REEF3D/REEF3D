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
    fd("U",'S',b->U.V); fd("V",'S',b->V.V); fd("W",'S',b->W.V); fd("UH",'S',b->UH.V); fd("VH",'S',b->VH.V);
    fd("WH",'S',b->WH.V); fd("UA",'S',b->UA.V); fd("VA",'S',b->VA.V); fd("MX",'S',b->MX.V);
    fd("MY",'S',b->MY.V); fd("Fx",'S',b->Fx.V); fd("Fs",'S',b->Fs.V); fd("Fn",'S',b->Fn.V);
    fd("FEx",'S',b->FEx.V); fd("Fy",'S',b->Fy.V); fd("Fe",'S',b->Fe.V); fd("Fw",'S',b->Fw.V);
    fd("FEy",'S',b->FEy.V); fd("Ss",'S',b->Ss.V); fd("Sn",'S',b->Sn.V); fd("SSx",'S',b->SSx.V);
    fd("Se",'S',b->Se.V); fd("Sw",'S',b->Sw.V); fd("SSy",'S',b->SSy.V); fd("Un",'S',b->Un.V);
    fd("Us",'S',b->Us.V); fd("Ue",'S',b->Ue.V); fd("Uw",'S',b->Uw.V); fd("Vn",'S',b->Vn.V);
    fd("Vs",'S',b->Vs.V); fd("Ve",'S',b->Ve.V); fd("Vw",'S',b->Vw.V); fd("Wn",'S',b->Wn.V);
    fd("Ws",'S',b->Ws.V); fd("We",'S',b->We.V); fd("Ww",'S',b->Ww.V); fd("UHn",'S',b->UHn.V);
    fd("UHs",'S',b->UHs.V); fd("UHe",'S',b->UHe.V); fd("UHw",'S',b->UHw.V); fd("VHn",'S',b->VHn.V);
    fd("VHs",'S',b->VHs.V); fd("VHe",'S',b->VHe.V); fd("VHw",'S',b->VHw.V); fd("WHn",'S',b->WHn.V);
    fd("WHs",'S',b->WHs.V); fd("WHe",'S',b->WHe.V); fd("WHw",'S',b->WHw.V); fd("QUn",'S',b->QUn.V);
    fd("QUs",'S',b->QUs.V); fd("QUe",'S',b->QUe.V); fd("QUw",'S',b->QUw.V); fd("QVn",'S',b->QVn.V);
    fd("QVs",'S',b->QVs.V); fd("QVe",'S',b->QVe.V); fd("QVw",'S',b->QVw.V); fd("ETAs",'S',b->ETAs.V);
    fd("ETAn",'S',b->ETAn.V); fd("ETAe",'S',b->ETAe.V); fd("ETAw",'S',b->ETAw.V); fd("Ds",'S',b->Ds.V);
    fd("Dn",'S',b->Dn.V); fd("De",'S',b->De.V); fd("Dw",'S',b->Dw.V); fd("dfx",'S',b->dfx.V);
    fd("dfy",'S',b->dfy.V); fd("WL",'S',b->WL.V); fd("eta",'S',b->eta.V); fd("eta_n",'S',b->eta_n.V);
    fd("P",'S',b->P.V); fd("Pn",'S',b->Pn.V); fd("F",'S',b->F.V); fd("Q",'S',b->Q.V); fd("Qn",'S',b->Qn.V);
    fd("G",'S',b->G.V); fd("H",'S',b->H.V); fd("L",'S',b->L.V); fd("Fext",'S',b->Fext.V);
    fd("Gext",'S',b->Gext.V); fd("Hext",'S',b->Hext.V); fd("ws",'S',b->ws.V); fd("press",'S',b->press.V);
    fd("eddyv",'S',b->eddyv.V); fd("kin",'S',b->kin.V); fd("eps",'S',b->eps.V); fd("hx",'S',b->hx.V);
    fd("hy",'S',b->hy.V); fd("hp",'S',b->hp.V); fd("dpx",'S',b->dpx.V); fd("dpy",'S',b->dpy.V);
    fd("test",'S',b->test.V); fd("Hs",'S',b->Hs.V); fd("fs",'S',b->fs.V); fd("fx",'S',b->fx.V);
    fd("fy",'S',b->fy.V); fd("fz",'S',b->fz.V); fd("bed",'S',b->bed.V); fd("bed0",'S',b->bed0.V);
    fd("depth",'S',b->depth.V); fd("solidbed",'S',b->solidbed.V); fd("topobed",'S',b->topobed.V);
    fd("bednode",'S',b->bednode.V); fd("breaking_print",'S',b->breaking_print.V); fd("ks",'S',b->ks.V);
    fi("breaking",'S',b->breaking.V); fi("wet1",'S',b->wet1.V);
    fi("deep1",'S',b->deep1.V); fi("wet2",'S',b->wet2.V); fi("deep2",'S',b->deep2.V);

    // the wet/dry flags of the lexer
    fi("wet",'S',p->wet); fi("lwet_n",'S',p->wet_n); fi("deep",'S',p->deep);
}

#endif
