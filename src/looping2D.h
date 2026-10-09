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

#ifndef LOOPING2D_H_
#define LOOPING2D_H_

// SLICE CONDITIONS
#define PSLICECHECK1  if(p->flagslice1[IJ]>0)
#define PSLICECHECK2  if(p->flagslice2[IJ]>0)
#define PSLICECHECK4  if(p->flagslice4[IJ]>0)
#define SSLICECHECK4  if(p->flagslice4[IJ]<0)

#define SEDSLICECHECK if(s->DFBED[IJ]>0)
#define SLICEFLEXCHECK  if(flagslice[IJ]>0)

#define WETDRY if(p->wet[IJ]==1)
#define WETDRYDEEP if(p->wet[IJ]==1 && p->deep[IJ]==1)

// SLICE BASE LOOPS
#define SLICEBASELOOP ILOOP JLOOP
#define JILOOP JLOOP ILOOP

// SLICE LOOPS
#define SLICELOOP1 IULOOP JLOOP  PSLICECHECK1
#define SLICELOOP2 ILOOP JVLOOP  PSLICECHECK2
#define SLICELOOP4 SLICEBASELOOP  PSLICECHECK4

#define SEDSLICELOOP SLICEBASELOOP  PSLICECHECK4 SEDSLICECHECK

#define TPSLICELOOP ITPLOOP JTPLOOP


#define SLICEFLEXLOOP IFLEXLOOP JFLEXLOOP SLICEFLEXCHECK

// GCBSL

#define GCSLB1 for(n=0;n<p->gcbsl1_count;++n)
#define GCSLB1CHECK if(p->gcbsl1[n][3]>0)
#define GCSL1LOOP GCSLB1 GCSLB1CHECK

#define QGCSLB1 for(q=0;q<p->gcbsl1_count;++q)
#define QGCSLB1CHECK if(p->gcbsl1[q][3]>0)

#define QQGCSLB1 for(qq=0;qq<p->gcbsl1_count;++qq)
#define QQGCSLB1CHECK if(p->gcbsl1[qq][3]>0)
#define QQGCSL1LOOP QQGCSLB1 QQGCSLB1CHECK


#define GCSLB2 for(n=0;n<p->gcbsl2_count;++n)
#define GCSLB2CHECK if(p->gcbsl2[n][3]>0)
#define GCSL2LOOP GCSLB2 GCSLB2CHECK

#define QGCSLB2 for(q=0;q<p->gcbsl2_count;++q)
#define QGCSLB2CHECK if(p->gcbsl2[q][3]>0)

#define QQGCSLB2 for(qq=0;qq<p->gcbsl2_count;++qq)
#define QQGCSLB2CHECK if(p->gcbsl2[qq][3]>0)
#define QQGCSL2LOOP QQGCSLB2 QQGCSLB2CHECK


#define GCSLB4 for(n=0;n<p->gcbsl4_count;++n)
#define GCSLB4CHECK if(p->gcbsl4[n][3]>0)
#define GCSL4LOOP GCSLB4 GCSLB4CHECK

#define QGCSLB4 for(q=0;q<p->gcbsl4_count;++q)
#define QGCSLB4CHECK if(p->gcbsl4[q][3]>0)

#define QQGCSLB4 for(qq=0;qq<p->gcbsl4_count;++qq)
#define QQGCSLB4CHECK if(p->gcbsl4[qq][3]>0)
#define QQGCSL4LOOP QQGCSLB4 QQGCSLB4CHECK


#define GCSLB4A for(n=0;n<p->gcbsl4a_count;++n)
#define GCSLB4ACHECK if(p->gcbsl4a[n][3]>0)

#define QGCSLB4A for(q=0;q<p->gcbsl4a_count;++q)
#define QGCSLB4ACHECK if(p->gcbsl4a[q][3]>0)

#define QQGCSLB4A for(qq=0;qq<p->gcbsl4a_count;++qq)
#define QQGCSLB4ACHECK if(p->gcbsl4a[qq][3]>0)
#define QQGCSL4ALOOP QQGCSLB4A QQGCSLB4ACHECK

#define GCSLDFETA4 for(n=0;n<p->gcsldfeta4_count;++n)
#define GCSLDFETA4CHECK if(p->gcsldfeta4[n][3]>0)
#define GCSLDFETA4LOOP GCSLDFETA4 GCSLDFETA4CHECK

#define GCSLDFBED4 for(n=0;n<p->gcsldfbed4_count;++n)
#define GCSLDFBED4CHECK if(p->gcsldfbed4[n][3]>0)
#define GCSLDFBED4LOOP GCSLDFBED4 GCSLDFBED4CHECK

#endif
