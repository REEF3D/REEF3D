/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"convection_void.h"
#include"fou.h"
#include"ifou.h"
#include"cds2.h"
#include"cds2_alt.h"
#include"cds4.h"
#include"icds2.h"
#include"icds4.h"
#include"quick.h"
#include"iquick.h"
#include"weno_hj.h"
#include"weno_hj_N.h"
#include"weno_hj_nug.h"
#include"weno_flux.h"
#include"weno_flux_nug.h"
#include"weno_flux_N.h"
#include"iweno_flux.h"
#include"iweno_hj.h"
#include"iweno_hj_nug.h"
#include"weno3_hj.h"
#include"weno3_flux.h"
#include"diff_void.h"
#include"ediff2.h"
#include"idiff2.h"
#include"idiff2_FS.h"

#include"hires.h"
#include"ihires.h"

#include"hric.h"
#include"hric_mod.h"
#include"cicsam.h"

#include"potential_v.h"
#include"potential_f.h"

