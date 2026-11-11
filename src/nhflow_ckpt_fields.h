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

#ifndef NHFLOW_CKPT_FIELDS_H_
#define NHFLOW_CKPT_FIELDS_H_

// the state of an NHFLOW grid at the end of a step (checkpoints, P 48, I 45): level 0 and every
// patch of the mesh refinement.  fd/fi(name, kind, pointer): kind S slice, F node field (FIJK: P,
// omegaF, the sigma grid), C cell field (IJK); nullptr arrays are skipped.  The turbulence model adds
// its own arrays (nhflow_turbulence::ckpt_fields).

#include"lexer.h"
#include"fdm_nhf.h"

template<class FD, class FI>
void nhflow_ckpt_fields(lexer *p, fdm_nhf *d, FD &&fd, FI &&fi)
{
    // slices
    fd("eta",'S',d->eta.data()); fd("eta_n",'S',d->eta_n.data()); fd("WL",'S',d->WL.data()); fd("detadt",'S',d->detadt.data());
    fd("detadt_n",'S',d->detadt_n.data()); fd("un",'S',d->un.data()); fd("vn",'S',d->vn.data()); fd("dudt",'S',d->dudt.data());
    fd("bed",'S',d->bed.data()); fd("depth",'S',d->depth.data()); fd("K",'S',d->K.data());
    fd("Ex",'S',d->Ex.data()); fd("Ey",'S',d->Ey.data()); fd("Exx",'S',d->Exx.data()); fd("Eyy",'S',d->Eyy.data());
    fd("Bx",'S',d->Bx.data()); fd("By",'S',d->By.data()); fd("Bxx",'S',d->Bxx.data()); fd("Byy",'S',d->Byy.data());
    fd("hx",'S',d->hx.data()); fd("hy",'S',d->hy.data()); fd("ks",'S',d->ks.data()); fd("coastline",'S',d->coastline.data());
    fd("vb",'S',d->vb.data()); fd("fs",'S',d->fs.data()); fd("breaking_print",'S',d->breaking_print.data()); fd("Hs",'S',d->Hs.data());
    fd("ETAs",'S',d->ETAs.data()); fd("ETAn",'S',d->ETAn.data()); fd("ETAe",'S',d->ETAe.data()); fd("ETAw",'S',d->ETAw.data());
    fd("Ds",'S',d->Ds.data()); fd("Dn",'S',d->Dn.data()); fd("De",'S',d->De.data()); fd("Dw",'S',d->Dw.data());
    fd("dfx",'S',d->dfx.data()); fd("dfy",'S',d->dfy.data());
    fi("etaloc",'S',d->etaloc.data()); fi("wet_n",'S',d->wet_n.data()); fi("breaking",'S',d->breaking.data());
    fi("breaklog",'S',d->breaklog.data()); fi("bc",'S',d->bc.data());

    // velocities, momenta, pressure, the cell fields of the step
    fd("U",'C',d->U); fd("V",'C',d->V); fd("W",'C',d->W); fd("omegaF",'F',d->omegaF);
    fd("UH",'C',d->UH); fd("VH",'C',d->VH); fd("WH",'C',d->WH); fd("P",'F',d->P);
    fd("RO",'C',d->RO); fd("VISC",'C',d->VISC); fd("EV",'C',d->EV); fd("EV0",'C',d->EV0); fd("EVH",'C',d->EVH);
    fd("F",'C',d->F); fd("G",'C',d->G); fd("H",'C',d->H); fd("L",'C',d->L);
    fd("Fext",'C',d->Fext); fd("Gext",'C',d->Gext); fd("Hext",'C',d->Hext);
    fd("POR",'C',d->POR); fd("PORPART",'C',d->PORPART); fd("PORDEM",'C',d->PORDEM);
    fd("KIN",'C',d->KIN); fd("CONC",'C',d->CONC); fd("SALT",'C',d->SALT); fd("TEMP",'C',d->TEMP); fd("SOLID",'C',d->SOLID); fd("FB",'C',d->FB); fd("FHB",'C',d->FHB);
    fd("PORSTRUC",'C',d->PORSTRUC); fd("Fx",'C',d->Fx); fd("Fy",'C',d->Fy); fd("Fz",'C',d->Fz);
    fd("FEx",'C',d->FEx); fd("FEy",'C',d->FEy);

    // the sigma grid and the flags of the lexer that change in time (wetting and drying)
    fd("sig",'F',p->sig); fd("sigx",'F',p->sigx); fd("sigy",'F',p->sigy); fd("sigt",'F',p->sigt); fd("sigxx",'F',p->sigxx);
    fd("sigz",'S',p->sigz); fd("pWL",'S',p->WL); fd("ZSN",'F',p->ZSN); fd("ZSP",'C',p->ZSP);
    fi("flagslice4",'S',p->flagslice4); fi("flag7",'F',p->flag7); fi("flag4",'C',p->flag4);
    fi("wet",'S',p->wet); fi("lwet_n",'S',p->wet_n); fi("deep",'S',p->deep);

    // the direct-forcing flags of the solids and floating bodies (nhflow_forcing, at the last forcing)
    fi("DF",'C',p->DF); fi("DFF",'F',p->DFF);
}

#endif
