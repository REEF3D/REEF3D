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


#define ILOOP	for(i=0; i<p->knox; ++i)
#define JLOOP	for(j=0; j<p->knoy; ++j)
#define KLOOP 	for(k=0; k<p->knoz; ++k)
#define PCHECK  if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
#define LOOP ILOOP JLOOP KLOOP PCHECK


#define ISTARTCHECK	if(si>=0)
#define IENDCHECK	if(si<p->knox)
#define JSTARTCHECK	if(sj>=0)
#define JENDCHECK	if(sj<p->knoy)
#define KSTARTCHECK	if(sk>=0)
#define KENDCHECK	if(sk<p->knoz)
#define SOLIDCHECK  if(p->flag4[(si-p->imin)*p->jmax*p->kmax + (sj-p->jmin)*p->kmax + sk-p->kmin]!=OBJ)
#define BFBCHK ISTARTCHECK IENDCHECK JSTARTCHECK JENDCHECK KSTARTCHECK KENDCHECK SOLIDCHECK


