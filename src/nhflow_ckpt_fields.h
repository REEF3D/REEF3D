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
    fd("eta",'S',d->eta.V); fd("eta_n",'S',d->eta_n.V); fd("WL",'S',d->WL.V); fd("detadt",'S',d->detadt.V);
    fd("detadt_n",'S',d->detadt_n.V); fd("un",'S',d->un.V); fd("vn",'S',d->vn.V); fd("dudt",'S',d->dudt.V);
    fd("bed",'S',d->bed.V); fd("depth",'S',d->depth.V); fd("K",'S',d->K.V);
    fd("Ex",'S',d->Ex.V); fd("Ey",'S',d->Ey.V); fd("Exx",'S',d->Exx.V); fd("Eyy",'S',d->Eyy.V);
    fd("Bx",'S',d->Bx.V); fd("By",'S',d->By.V); fd("Bxx",'S',d->Bxx.V); fd("Byy",'S',d->Byy.V);
    fd("hx",'S',d->hx.V); fd("hy",'S',d->hy.V); fd("ks",'S',d->ks.V); fd("coastline",'S',d->coastline.V);
    fd("vb",'S',d->vb.V); fd("fs",'S',d->fs.V); fd("breaking_print",'S',d->breaking_print.V); fd("Hs",'S',d->Hs.V);
    fd("ETAs",'S',d->ETAs.V); fd("ETAn",'S',d->ETAn.V); fd("ETAe",'S',d->ETAe.V); fd("ETAw",'S',d->ETAw.V);
    fd("Ds",'S',d->Ds.V); fd("Dn",'S',d->Dn.V); fd("De",'S',d->De.V); fd("Dw",'S',d->Dw.V);
    fd("dfx",'S',d->dfx.V); fd("dfy",'S',d->dfy.V);
    fi("etaloc",'S',d->etaloc.V); fi("wet_n",'S',d->wet_n.V); fi("breaking",'S',d->breaking.V);
    fi("breaklog",'S',d->breaklog.V); fi("bc",'S',d->bc.V);

    // velocities, momenta, pressure, the cell fields of the step
    fd("U",'C',d->U); fd("V",'C',d->V); fd("W",'C',d->W); fd("omegaF",'F',d->omegaF);
    fd("UH",'C',d->UH); fd("VH",'C',d->VH); fd("WH",'C',d->WH); fd("P",'F',d->P);
    fd("RO",'C',d->RO); fd("VISC",'C',d->VISC); fd("EV",'C',d->EV); fd("EV0",'C',d->EV0); fd("EVH",'C',d->EVH);
    fd("F",'C',d->F); fd("G",'C',d->G); fd("H",'C',d->H); fd("L",'C',d->L);
    fd("Fext",'C',d->Fext); fd("Gext",'C',d->Gext); fd("Hext",'C',d->Hext);
    fd("POR",'C',d->POR); fd("PORPART",'C',d->PORPART); fd("PORDEM",'C',d->PORDEM);
    fd("KIN",'C',d->KIN); fd("CONC",'C',d->CONC); fd("SALT",'C',d->SALT); fd("TEMP",'C',d->TEMP); fd("SOLID",'C',d->SOLID); fd("FB",'C',d->FB); fd("FHB",'C',d->FHB);
    fd("PORSTRUC",'C',d->PORSTRUC); fd("Fx",'C',d->Fx); fd("Fy",'C',d->Fy); fd("Fz",'C',d->Fz);
    fd("FEx",'C',d->FEx); fd("FEy",'C',d->FEy); fd("FSW",'C',d->FSW); fd("DWDT",'C',d->DWDT);

    // the sigma grid and the flags of the lexer that change in time (wetting and drying)
    fd("sig",'F',p->sig); fd("sigx",'F',p->sigx); fd("sigy",'F',p->sigy); fd("sigt",'F',p->sigt); fd("sigxx",'F',p->sigxx);
    fd("sigz",'S',p->sigz); fd("pWL",'S',p->WL); fd("ZSN",'F',p->ZSN); fd("ZSP",'C',p->ZSP);
    fi("flagslice4",'S',p->flagslice4); fi("flag7",'F',p->flag7); fi("flag4",'C',p->flag4);
    fi("wet",'S',p->wet); fi("lwet_n",'S',p->wet_n); fi("deep",'S',p->deep);

    // the direct-forcing flags of the solids and floating bodies (nhflow_forcing, at the last forcing)
    fi("DF",'C',p->DF); fi("DFF",'F',p->DFF);
}

#endif
