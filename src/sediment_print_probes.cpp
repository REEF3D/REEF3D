/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"bedshear.h"
#include"bedprobe_point.h"
#include"bedprobe_max.h"
#include"bedshear_probe.h"
#include"bedshear_max.h"
#include"bedprobe_line_x.h"
#include"bedprobe_line_y.h"

void sediment_f::print_probes(lexer *p, ghostcell *pgc, sediment_fdm *s, ioflow *pflow)
{

// sediment probes
	if(p->P121>0)
	pbedpt->bed_gauge(p,pgc,s);

	if(p->P122>0)
	pbedmax->bed_max(p,pgc,s);

	if(p->P123>0)
	pbedlinex->start(p,pgc,s,pflow);

	if(p->P124>0)
	pbedliney->start(p,pgc,s,pflow);

	if(p->P125>0)
	pbedshearprobe->bedshear_gauge(p,pgc,psed);

	if(p->P126>0)
	pbedshearmax->bedshear_maxval(p,pgc,psed);    
}
