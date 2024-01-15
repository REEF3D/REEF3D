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

#include"iterators.h"
#include"boundcheck.h"
#include"definitions.h"
#include"looping2D.h"
#include"iterators2D.h"
#include"iterators1D.h"

#define PI 3.14159265359
#define EE 2.71828182846

#define ILOOP	for(i=0; i<p->knox; ++i)
#define JLOOP	for(j=0; j<p->knoy; ++j)
#define KLOOP 	for(k=0; k<p->knoz; ++k)
#define IREVLOOP	for(i=p->knox-1; i>=0; --i)
#define JREVLOOP	for(j=p->knoy-1; j>=0; --j)
#define KREVLOOP 	for(k=p->knoz-1; k>=0; --k)
#define PCHECK  if(p->flag4[IJK]>0)
#define LOOP ILOOP JLOOP KLOOP PCHECK

#define FSCHECK  if(p->flag7[FIJK]<=0)
#define FPCHECK  if(p->flag7[FIJK]>0)
#define SCHECK  if(p->flag4[IJK]<0)
#define SFLUIDCHECK  if(p->flag4[IJK]<AIR)
#define USCHECK  if(p->flag1[IJK]<0)
#define VSCHECK  if(p->flag2[IJK]<0)
#define WSCHECK  if(p->flag3[IJK]<0)

#define KJILOOP KLOOP JLOOP ILOOP
#define UKJILOOP KLOOP JLOOP IULOOP
#define VKJILOOP KLOOP JVLOOP ILOOP
#define WKJILOOP KWLOOP JLOOP ILOOP
#define JILOOP JLOOP ILOOP

#define FKJILOOP FKLOOP FJLOOP FILOOP

#define WETDRYCHK if(p->wet[IJ]>0)
#define FILOOPWD ILOOP JLOOP KLOOP WETDRYCHK PCHECK

#define PLAINLOOP ILOOP JLOOP KLOOP

#define IMALOOP	for(i=-p->margin; i<p->knox+p->margin; ++i)
#define JMALOOP	for(j=-p->margin; j<p->knoy+p->margin; ++j)
#define KMALOOP 	for(k=-p->margin; k<p->knoz+p->margin; ++k)
#define MALOOP IMALOOP JMALOOP KMALOOP

#define IBLOOP	for(i=-1; i<p->knox+1; ++i)
#define JBLOOP	for(j=-1; j<p->knoy+1; ++j)
#define KBLOOP for(k=-1; k<p->knoz+1; ++k)
#define BLOOP IBLOOP JBLOOP KBLOOP

#define IBLOOP	for(i=-1; i<p->knox+1; ++i)
#define JBLOOP	for(j=-1; j<p->knoy+1; ++j)
#define KBLOOP for(k=-1; k<p->knoz+1; ++k)
#define BBASELOOP IBLOOP JBLOOP KBLOOP PBASECHECK

#define ITLOOP for(i=0; i<p->knox+1; ++i)
#define JTLOOP for(j=0; j<p->knoy+1; ++j)
#define KTLOOP for(k=0; k<p->knoz+1; ++k)
#define TCHECK if(p->flag4[IJK]>OBJ)
#define TLOOP ITLOOP JTLOOP KTLOOP

#define IFLEXLOOP	for(i=0; i<p->knox-ulast; ++i)
#define JFLEXLOOP	for(j=0; j<p->knoy-vlast; ++j)
#define KFLEXLOOP	for(k=0; k<p->knoz-wlast; ++k)
#define FLEXCHECK  if(flag[IJK]>0)
#define FLEXLOOP IFLEXLOOP JFLEXLOOP KFLEXLOOP FLEXCHECK

#define IULOOP	for(i=0; i<p->knox-p->ulast; ++i)
#define JULOOP	for(j=0; j<p->knoy; ++j)
#define KULOOP	for(k=0; k<p->knoz; ++k)
#define UCHECK  if(p->flag1[IJK]>0)
#define ULOOP IULOOP JULOOP KULOOP UCHECK

#define IVLOOP	for(i=0; i<p->knox; ++i)
#define JVLOOP	for(j=0; j<p->knoy-p->vlast; ++j)
#define KVLOOP	for(k=0; k<p->knoz; ++k)
#define VCHECK  if(p->flag2[IJK]>0)
#define VLOOP IVLOOP JVLOOP KVLOOP VCHECK

#define IWLOOP	for(i=0; i<p->knox; ++i)
#define JWLOOP	for(j=0; j<p->knoy; ++j)
#define KWLOOP	for(k=0; k<p->knoz-p->wlast; ++k)
#define WCHECK  if(p->flag3[IJK]>0)
#define WLOOP IWLOOP JWLOOP KWLOOP WCHECK

#define UBASECHECK  if(p->flag1[IJK]>OBJ)
#define VBASECHECK  if(p->flag2[IJK]>OBJ)
#define WBASECHECK  if(p->flag3[IJK]>OBJ)
#define PBASECHECK  if(p->flag4[IJK]>OBJ)

#define UBASELOOP IULOOP JULOOP KULOOP UBASECHECK
#define VBASELOOP IVLOOP JVLOOP KVLOOP VBASECHECK
#define WBASELOOP IWLOOP JWLOOP KWLOOP WBASECHECK
#define BASELOOP ILOOP JLOOP KLOOP PBASECHECK

#define URAWLOOP IULOOP JULOOP KULOOP
#define VRAWLOOP IVLOOP JVLOOP KVLOOP 
#define WRAWLOOP IWLOOP JWLOOP KWLOOP 
#define RAWLOOP ILOOP JLOOP KLOOP

#define ALOOP ILOOP JLOOP KLOOP PSOLIDCHECK

#define USOLIDCHECK  if(p->flag1[IJK]>SOLID)
#define VSOLIDCHECK  if(p->flag2[IJK]>SOLID)
#define WSOLIDCHECK  if(p->flag3[IJK]>SOLID)
#define PSOLIDCHECK  if(p->flag4[IJK]>SOLID)

#define USOLIDLOOP IULOOP JULOOP KULOOP USOLIDCHECK
#define VSOLIDLOOP IVLOOP JVLOOP KVLOOP VSOLIDCHECK
#define WSOLIDLOOP IWLOOP JWLOOP KWLOOP WSOLIDCHECK
#define SOLIDLOOP ILOOP JLOOP KLOOP PSOLIDCHECK    

#define UAIRCHECK  if(p->flag1[IJK]==AIR)
#define VAIRCHECK  if(p->flag2[IJK]==AIR)
#define WAIRCHECK  if(p->flag3[IJK]==AIR)
#define PAIRCHECK  if(p->flag4[IJK]==AIR)

#define UAIRLOOP IULOOP JULOOP KULOOP UAIRCHECK
#define VAIRLOOP IVLOOP JVLOOP KVLOOP VAIRCHECK
#define WAIRLOOP IWLOOP JWLOOP KWLOOP WAIRCHECK
#define AIRLOOP ILOOP JLOOP KLOOP PAIRCHECK    

#define UFLUIDCHECK  if(p->flag1[IJK]>=AIR)
#define VFLUIDCHECK  if(p->flag2[IJK]>=AIR)
#define WFLUIDCHECK  if(p->flag3[IJK]>=AIR)
#define PFLUIDCHECK  if(p->flag4[IJK]>=AIR)
    
#define PWDFLUIDCHECK  if(p->flag4[IJK]>=AIR && p->wet[IJ]>0)
#define FSWDCHECK  if(p->flag7[FIJK]<=0 || p->wet[IJ]==0)
#define FPWDCHECK  if(p->flag7[FIJK]>0  && p->wet[IJ]>0)

#define UFLUIDLOOP IULOOP JULOOP KULOOP UFLUIDCHECK
#define VFLUIDLOOP IVLOOP JVLOOP KVLOOP VFLUIDCHECK
#define WFLUIDLOOP IWLOOP JWLOOP KWLOOP WFLUIDCHECK
#define FLUIDLOOP ILOOP JLOOP KLOOP PFLUIDCHECK  

#define FILOOP	for(i=0; i<p->knox; ++i)
#define FJLOOP	for(j=0; j<p->knoy; ++j)
#define FKLOOP for(k=0; k<p->knoz+1; ++k)
#define FLOOP FILOOP FJLOOP FKLOOP FPCHECK  
#define FBASELOOP FILOOP FJLOOP FKLOOP 
    

#define ETALOC for(k=a->etaloc(i,j); k<a->etaloc(i,j)+1; ++k)
#define FILOOP4 ILOOP JLOOP ETALOC PFLUIDCHECK 

#define FETALOC for(k=c->etaloc(i,j); k<c->etaloc(i,j)+1; ++k)
#define FFILOOP4 ILOOP JLOOP FETALOC FPCHECK


#define ITPLOOP for(i=-1; i<p->knox; ++i)
#define JTPLOOP for(j=-1; j<p->knoy; ++j)
#define KTPLOOP for(k=-1; k<p->knoz; ++k)
#define TPCHECK  if(p->tpflag[IJK]>0)
#define TPLOOP ITPLOOP JTPLOOP KTPLOOP TPCHECK

#define NDBASECHECK  if(p->ndbaseflag[IJK]>OBJ)
#define NDBASELOOP ITPLOOP JTPLOOP KTPLOOP NDBASECHECK

#define MAX(aAa,bBb) ((aAa)>(bBb)?(aAa):(bBb))
#define MIN(aAa,bBb) ((aAa)<(bBb)?(aAa):(bBb))
#define SIGN(aAa) ((aAa)>= 0.0 ? 1.0 : -1.0)

#define MAX3(aAa,bBb,cCc) (((aAa)>(bBb)?(aAa):(bBb))>cCc?((aAa)>(bBb)?(aAa):(bBb)):cCc)
#define MIN3(aAa,bBb,cCc) (((aAa)<(bBb)?(aAa):(bBb))<cCc?((aAa)<(bBb)?(aAa):(bBb)):cCc)

#define GCB1 for(n=0;n<p->gcb1_count;++n)
#define GCB1CHECK if(p->gcb1[n][3]>0)
#define GC1LOOP GCB1 GCB1CHECK

#define GGCB1 for(g=0;g<p->gcb1_count;++g)
#define GGCB1CHECK if(p->gcb1[g][3]>0)
#define GGC1LOOP GGCB1 GGCB1CHECK

#define QGCB1 for(q=0;q<p->gcb1_count;++q)
#define QGCB1CHECK if(p->gcb1[q][3]>0)
#define QGC1LOOP QGCB1 QGCB1CHECK

#define QQGCB1 for(qq=0;qq<p->gcb1_count;++qq)
#define QQGCB1CHECK if(p->gcb1[qq][3]>0)
#define QQGC1LOOP QQGCB1 QQGCB1CHECK

#define GCB2 for(n=0;n<p->gcb2_count;++n)
#define GCB2CHECK if(p->gcb2[n][3]>0)
#define GC2LOOP GCB2 GCB2CHECK

#define GGCB2 for(g=0;g<p->gcb2_count;++g)
#define GGCB2CHECK if(p->gcb2[g][3]>0)
#define GGC2LOOP GGCB2 GGCB2CHECK

#define QGCB2 for(q=0;q<p->gcb2_count;++q)
#define QGCB2CHECK if(p->gcb2[q][3]>0)
#define QGC2LOOP QGCB2 QGCB2CHECK

#define QQGCB2 for(qq=0;qq<p->gcb2_count;++qq)
#define QQGCB2CHECK if(p->gcb2[qq][3]>0)
#define QQGC2LOOP QQGCB2 QQGCB2CHECK

#define GCB3 for(n=0;n<p->gcb3_count;++n)
#define GCB3CHECK if(p->gcb3[n][3]>0)
#define GC3LOOP GCB3 GCB3CHECK

#define GGCB3 for(g=0;g<p->gcb3_count;++g)
#define GGCB3CHECK if(p->gcb3[g][3]>0)
#define GGC3LOOP GGCB3 GGCB3CHECK

#define QGCB3 for(q=0;q<p->gcb3_count;++q)
#define QGCB3CHECK if(p->gcb3[q][3]>0)
#define QGC3LOOP QGCB3 QGCB3CHECK

#define QQGCB3 for(qq=0;qq<p->gcb3_count;++qq)
#define QQGCB3CHECK if(p->gcb3[qq][3]>0)
#define QQGC3LOOP QQGCB3 QQGCB3CHECK

#define GCB4 for(n=0;n<p->gcb4_count;++n)
#define GCB4CHECK if(p->gcb4[n][3]>0)
#define GC4LOOP GCB4 GCB4CHECK

#define GGCB4 for(g=0;g<p->gcb4_count;++g)
#define GGCB4CHECK if(p->gcb4[g][3]>0)
#define GGC4LOOP GGCB4 GGCB4CHECK

#define QGCB4 for(q=0;q<p->gcb4_count;++q)
#define QGCB4CHECK if(p->gcb4[q][3]>0)
#define QGC4LOOP QGCB4 QGCB4CHECK

#define QQGCB4 for(qq=0;qq<p->gcb4_count;++qq)
#define QQGCB4CHECK if(p->gcb4[qq][3]>0)
#define QQGC4LOOP QQGCB4 QQGCB4CHECK

//df
#define QGCDF4 for(q=0;q<p->gcdf4_count;++q)
#define QGCDF4CHECK if(p->gcdf4[q][3]>0)
#define QGCDF4LOOP QGCDF4 QGCDF4CHECK

#define GCDF4 for(n=0;n<p->gcdf4_count;++n)
#define GCDF4CHECK if(p->gcdf4[n][3]>0)
#define GCDF4LOOP GCDF4 GCDF4CHECK


#define GCB for(n=0;n<gcb_count;++n)
#define GCBCHECK if(gcb[n][3]>0)
#define GCLOOP GCB GCBCHECK

#define GGCB for(g=0;g<gcb_count;++g)
#define GGCBCHECK if(gcb[g][3]>0)
#define GGCLOOP GGCB GGCBCHECK

#define QGCB for(q=0;q<gcb_count;++q)
#define QGCBCHECK if(gcb[q][3]>0)
#define QGCLOOP QGCB QGCBCHECK

#define QQGCB for(qq=0;qq<gcb_count;++qq)
#define QQGCBCHECK if(gcb[qq][3]>0)
#define QQGCLOOP QQGCB QQGCBCHECK

#define GC4A  for(n=0;n<p->gcb4a_count;++n)
#define GCB4ACHECK if(p->gcb4a[n][3]>0)
#define GC4ALOOP  GC4A GCB4ACHECK

#define QGC4A  for(q=0;q<p->gcb4a_count;++q)
#define QGCB4ACHECK if(p->gcb4a[q][3]>0)
#define QGC4ALOOP  QGC4A QGCB4ACHECK

#define QQGC4A  for(qq=0;qq<p->gcb4a_count;++qq)
#define QQGCB4ACHECK if(p->gcb4a[qq][3]>0)
#define QQGC4ALOOP  QQGC4A QQGCB4ACHECK

#define GGC4A  for(g=0;g<p->gcb4a_count;++g)
#define GGCB4ACHECK if(p->gcb4a[g][3]>0)
#define GGC4ALOOP  GGC4A GGCB4ACHECK
        
#define GC6LOOP  for(n=0;n<p->gcb_fix;++n)
#define QGC6LOOP  for(q=0;q<p->gcb_fix;++q)
#define QQGC6LOOP  for(qq=0;qq<p->gcb_fix;++qq)
#define GGC6LOOP  for(g=0;g<p->gcb_fix;++g)


#define DT p->dt
#define NDT p->dt_old









