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

// LOOP

#define PSLICECHECK1  if(p->flagslice1[IJ]>0)
#define SLICELOOP1 IULOOP JULOOP  PSLICECHECK1

#define IULOOPSFLOW for(i=0; i<p->knox-p->ulastsflow; ++i)

#define PSLICECHECK2  if(p->flagslice2[IJ]>0)
#define SLICELOOP2 IVLOOP JVLOOP  PSLICECHECK2

#define PSLICECHECK4  if(p->flagslice4[IJ]>0)
#define SSLICECHECK4  if(p->flagslice4[IJ]<=0)
#define SLICELOOP4 ILOOP JLOOP  PSLICECHECK4

#define SLICEBASELOOP ILOOP JLOOP

#define SLICEREVLOOP4 IREVLOOP JREVLOOP  PSLICECHECK4

#define IREVLOOP	for(i=p->knox-1; i>=0; --i)
#define JREVLOOP	for(j=p->knoy-1; j>=0; --j)
#define SLICELOOPREV4 IREVLOOP JREVLOOP  PSLICECHECK4

#define PSLICECHECK4E  if(p->flagslice4[IJ]>=-1)
    
#define SLICELOOP4E ILOOP JLOOP  PSLICECHECK4

#define SLICEBASELOOP ILOOP JLOOP 


#define TPSLICECHECK  if(p->tpflagslice[IJ]>0)
#define TPSLICELOOP ITPLOOP JTPLOOP TPSLICECHECK

#define NSLICELOOP for(n=sizeS[0]; n<sizeS[1]; ++n)
#define NSLICELOOP1 for(n=p->sizeS1[0]; n<p->sizeS1[1]; ++n)
#define NSLICELOOP2 for(n=p->sizeS2[0]; n<p->sizeS2[1]; ++n)
#define NSLICELOOP4 for(n=p->sizeS4[0]; n<p->sizeS4[1]; ++n)
    
#define IFLEXLOOP	for(i=0; i<p->knox-ulast; ++i)
#define JFLEXLOOP	for(j=0; j<p->knoy-vlast; ++j)
#define KFLEXLOOP	for(k=0; k<p->knoz-wlast; ++k)
    
#define SLICEFLEXCHECK  if(flagslice[IJ]>0)
#define SLICEFLEXLOOP IFLEXLOOP JFLEXLOOP SLICEFLEXCHECK
    
#define SSLICECHECK4  if(p->flagslice4[IJ]<0)

#define WETDRY1 if(b->wet1(i,j)==1)
#define WETDRY2 if(b->wet2(i,j)==1)    
#define WETDRY if(p->wet[IJ]==1)
#define WETDRYDEEP if(p->wet[IJ]==1 && p->deep[IJ]==1)


// GCBSL

#define GCSLB1 for(n=0;n<p->gcbsl1_count;++n)
#define GCSLB1CHECK if(p->gcbsl1[n][3]>0)
#define GCSL1LOOP GCSLB1 GCSLB1CHECK

#define GGCSLB1 for(g=0;g<p->gcbsl1_count;++g)
#define GGCSLB1CHECK if(p->gcbsl1[g][3]>0)
#define GGCSL1LOOP GGCSLB1 GGCSLB1CHECK

#define QGCSLB1 for(q=0;q<p->gcbsl1_count;++q)
#define QGCSLB1CHECK if(p->gcbsl1[q][3]>0)
#define QGCSL1LOOP QGCSLB1 QGCSLB1CHECK

#define QQGCSLB1 for(qq=0;qq<p->gcbsl1_count;++qq)
#define QQGCSLB1CHECK if(p->gcbsl1[qq][3]>0)
#define QQGCSL1LOOP QQGCSLB1 QQGCSLB1CHECK


#define GCSLB2 for(n=0;n<p->gcbsl2_count;++n)
#define GCSLB2CHECK if(p->gcbsl2[n][3]>0)
#define GCSL2LOOP GCSLB2 GCSLB2CHECK

#define GGCSLB2 for(g=0;g<p->gcbsl2_count;++g)
#define GGCSLB2CHECK if(p->gcbsl2[g][3]>0)
#define GGCSL2LOOP GGCSLB2 GGCSLB2CHECK

#define QGCSLB2 for(q=0;q<p->gcbsl2_count;++q)
#define QGCSLB2CHECK if(p->gcbsl2[q][3]>0)
#define QGCSL2LOOP QGCSLB2 QGCSLB2CHECK

#define QQGCSLB2 for(qq=0;qq<p->gcbsl2_count;++qq)
#define QQGCSLB2CHECK if(p->gcbsl2[qq][3]>0)
#define QQGCSL2LOOP QQGCSLB2 QQGCSLB2CHECK


#define GCSLB3 for(n=0;n<p->gcbsl3_count;++n)
#define GCSLB3CHECK if(p->gcbsl3[n][3]>0)
#define GCSL3LOOP GCSLB3 GCSLB3CHECK

#define GGCSLB3 for(g=0;g<p->gcbsl3_count;++g)
#define GGCSLB3CHECK if(p->gcbsl3[g][3]>0)
#define GGCSL3LOOP GGCSLB3 GGCSLB3CHECK

#define QGCSLB3 for(q=0;q<p->gcbsl3_count;++q)
#define QGCSLB3CHECK if(p->gcbsl3[q][3]>0)
#define QGCSL3LOOP QGCSLB3 QGCSLB3CHECK

#define QQGCSLB3 for(qq=0;qq<p->gcbsl3_count;++qq)
#define QQGCSLB3CHECK if(p->gcbsl3[qq][3]>0)
#define QQGCSL3LOOP QQGCSLB3 QQGCSLB3CHECK


#define GCSLB4 for(n=0;n<p->gcbsl4_count;++n)
#define GCSLB4CHECK if(p->gcbsl4[n][3]>0)
#define GCSL4LOOP GCSLB4 GCSLB4CHECK

#define GGCSLB4 for(g=0;g<p->gcbsl4_count;++g)
#define GGCSLB4CHECK if(p->gcbsl4[g][3]>0)
#define GGCSL4LOOP GGCSLB4 GGCSLB4CHECK

#define QGCSLB4 for(q=0;q<p->gcbsl4_count;++q)
#define QGCSLB4CHECK if(p->gcbsl4[q][3]>0)
#define QGCSL4LOOP QGCSLB4 QGCSLB4CHECK

#define QQGCSLB4 for(qq=0;qq<p->gcbsl4_count;++qq)
#define QQGCSLB4CHECK if(p->gcbsl4[qq][3]>0)
#define QQGCSL4LOOP QQGCSLB4 QQGCSLB4CHECK


#define GCSLB4A for(n=0;n<p->gcbsl4a_count;++n)
#define GCSLB4ACHECK if(p->gcbsl4a[n][3]>0)
#define GCSL4ALOOP GCSLB4A GCSLB4ACHECK

#define GGCSLB4A for(g=0;g<p->gcbsl4a_count;++g)
#define GGCSLB4ACHECK if(p->gcbsl4a[g][3]>0)
#define GGCSL4ALOOP GGCSLB4A GGCSLB4ACHECK

#define QGCSLB4A for(q=0;q<p->gcbsl4a_count;++q)
#define QGCSLB4ACHECK if(p->gcbsl4a[q][3]>0)
#define QGCSL4ALOOP QGCSLB4A QGCSLB4ACHECK

#define QQGCSLB4A for(qq=0;qq<p->gcbsl4a_count;++qq)
#define QQGCSLB4ACHECK if(p->gcbsl4a[qq][3]>0)
#define QQGCSL4ALOOP QQGCSLB4A QQGCSLB4ACHECK

