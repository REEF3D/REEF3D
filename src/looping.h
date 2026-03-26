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

#ifndef LOOPING_H_
#define LOOPING_H_

#include"definitions.h"

#include"looping2D.h"

#include"iterators1D.h"
#include"iterators2D.h"
#include"iterators3D.h"

// LOOPs
#define ILOOP for(i=0; i<p->knox; ++i)
#define JLOOP for(j=0; j<p->knoy; ++j)
#define KLOOP for(k=0; k<p->knoz; ++k)

#define ITLOOP for(i=0; i<p->knox+1; ++i)
#define JTLOOP for(j=0; j<p->knoy+1; ++j)
#define KTLOOP for(k=0; k<p->knoz+1; ++k)

#define ITPLOOP for(i=-1; i<p->knox; ++i)
#define JTPLOOP for(j=-1; j<p->knoy; ++j)
#define KTPLOOP for(k=-1; k<p->knoz; ++k)

#define IBLOOP for(i=-1; i<p->knox+1; ++i)
#define JBLOOP for(j=-1; j<p->knoy+1; ++j)
#define KBLOOP for(k=-1; k<p->knoz+1; ++k)

#define IMALOOP for(i=-p->margin; i<p->knox+p->margin; ++i)
#define JMALOOP for(j=-p->margin; j<p->knoy+p->margin; ++j)
#define KMALOOP for(k=-p->margin; k<p->knoz+p->margin; ++k)

#define IFLEXLOOP for(i=0; i<p->knox-ulast; ++i)
#define JFLEXLOOP for(j=0; j<p->knoy-vlast; ++j)
#define KFLEXLOOP for(k=0; k<p->knoz-wlast; ++k)

#define IULOOP for(i=0; i<p->knox-p->ulast; ++i)
#define JVLOOP for(j=0; j<p->knoy-p->vlast; ++j)
#define KWLOOP for(k=0; k<p->knoz-p->wlast; ++k)

#define FILOOP ILOOP
#define FJLOOP JLOOP
#define FKLOOP for(k=0; k<p->knoz+1; ++k)

#define ETALOC  for(k=a->etaloc(i,j); k<a->etaloc(i,j)+1; ++k)
#define FETALOC for(k=c->etaloc(i,j); k<c->etaloc(i,j)+1; ++k)

// CONDITIONS
#define FLEXCHECK   if(flag[IJK]>0)
#define UCHECK      if(p->flag1[IJK]>0)
#define UFLUIDCHECK if(p->flag1[IJK]>=AIR_FLAG)
#define USCHECK     if(p->flag1[IJK]<0)
#define VCHECK      if(p->flag2[IJK]>0)
#define VFLUIDCHECK if(p->flag2[IJK]>=AIR_FLAG)
#define VSCHECK     if(p->flag2[IJK]<0)
#define WCHECK      if(p->flag3[IJK]>0)
#define WFLUIDCHECK if(p->flag3[IJK]>=AIR_FLAG)
#define WSCHECK     if(p->flag3[IJK]<0)
#define PCHECK      if(p->flag4[IJK]>0)
#define SCHECK      if(p->flag4[IJK]<0)
#define PFLUIDCHECK if(p->flag4[IJK]>=AIR_FLAG)
#define PAIR_CHECK  if(p->flag4[IJK]==AIR_FLAG)
#define SFLUIDCHECK if(p->flag4[IJK]<AIR_FLAG)
#define PSOLIDCHECK if(p->flag4[IJK]>SOLID_FLAG)
#define PBASECHECK  if(p->flag4[IJK]>OBJ_FLAG)
#define FPCHECK     if(p->flag7[FIJK]>0)
#define FPWDCHECK   if(p->flag7[FIJK]>0 && p->wet[IJ]>0)
#define FSCHECK     if(p->flag7[FIJK]<=0)
#define FSWDCHECK   if(p->flag7[FIJK]<=0 || p->wet[IJ]==0)

// POROSITY
#define PORVAL1 (0.5*(a->porosity(i+1,j,k) + a->porosity(i,j,k)))
#define PORVAL2 (0.5*(a->porosity(i,j+1,k) + a->porosity(i,j,k)))
#define PORVAL3 (0.5*(a->porosity(i,j,k+1) + a->porosity(i,j,k)))
#define PORVAL4 a->porosity(i,j,k)

#define PORVAL4px a->porosity(i+1,j,k)
#define PORVAL4py a->porosity(i,j+1,k)
#define PORVAL4pz a->porosity(i,j,k+1)

#define CPOR4px   (1.0/(1.0+(p->B260*(PORVAL4px<1.0?1.0:0.0))))
#define CPOR4py   (1.0/(1.0+(p->B260*(PORVAL4py<1.0?1.0:0.0))))
#define CPOR4pz   (1.0/(1.0+(p->B260*(PORVAL4pz<1.0?1.0:0.0))))

#define CPOR1   (1.0/(1.0+(p->B260*(PORVAL1<1.0?1.0:0.0))))
#define CPOR2   (1.0/(1.0+(p->B260*(PORVAL2<1.0?1.0:0.0))))
#define CPOR3   (1.0/(1.0+(p->B260*(PORVAL3<1.0?1.0:0.0))))
#define CPOR4   (1.0/(1.0+(p->B260*(PORVAL4<1.0?1.0:0.0))))

#define PORVAL1m (0.5*(a->porosity(i,j,k) + a->porosity(i-1,j,k)))
#define PORVAL2m (0.5*(a->porosity(i,j,k) + a->porosity(i,j-1,k)))
#define PORVAL3m (0.5*(a->porosity(i,j,k) + a->porosity(i,j,k-1)))

#define CPOR1m   (1.0/(1.0+(p->B260*(PORVAL1m<1.0?1.0:0.0))))
#define CPOR2m   (1.0/(1.0+(p->B260*(PORVAL2m<1.0?1.0:0.0))))
#define CPOR3m   (1.0/(1.0+(p->B260*(PORVAL3m<1.0?1.0:0.0))))

#define PORVAL1p (0.5*(a->porosity(i+2,j,k) + a->porosity(i+1,j,k)))
#define PORVAL2p (0.5*(a->porosity(i,j+2,k) + a->porosity(i,j+1,k)))
#define PORVAL3p (0.5*(a->porosity(i,j,k+2) + a->porosity(i,j,k+1)))

#define CPOR1p   (1.0/(1.0+(p->B260*(PORVAL1p<1.0?1.0:0.0))))
#define CPOR2p   (1.0/(1.0+(p->B260*(PORVAL2p<1.0?1.0:0.0))))
#define CPOR3p   (1.0/(1.0+(p->B260*(PORVAL3p<1.0?1.0:0.0))))

#define PORVALNH d->POR[IJK]
#define PORVALNHm d->POR[IJK]
#define PORVALNHp d->POR[IJK]
#define CPORNH  (1.0/(1.0+(p->B260*(PORVALNH<1.0?1.0:0.0))))
#define CPORNHm (1.0/(1.0+(p->B260*(PORVALNHm<1.0?1.0:0.0))))
#define CPORNHp (1.0/(1.0+(p->B260*(PORVALNHp<1.0?1.0:0.0))))

// COMBINDED LOOPS
#define IJKLOOP ILOOP JLOOP KLOOP
#define KJILOOP KLOOP JLOOP ILOOP

// MAIN LOOPS
#define PLAINLOOP IJKLOOP
#define LOOP PLAINLOOP PCHECK
#define BASELOOP PLAINLOOP PBASECHECK
#define BASEREVLOOP KJILOOP PBASECHECK
#define TPLOOP KTPLOOP JTPLOOP ITPLOOP

#define ALOOP PLAINLOOP PSOLIDCHECK
#define AIRLOOP PLAINLOOP PAIR_CHECK

#define MALOOP IMALOOP JMALOOP KMALOOP

// BOUNDARY LOOPS
#define BLOOP IBLOOP JBLOOP KBLOOP
#define BBASELOOP IBLOOP JBLOOP KBLOOP PBASECHECK

// FLUID LOOPS
#define ULOOP IULOOP JLOOP KLOOP UCHECK
#define VLOOP ILOOP JVLOOP KLOOP VCHECK
#define WLOOP ILOOP JLOOP KWLOOP WCHECK

#define UFLUIDLOOP IULOOP JLOOP KLOOP UFLUIDCHECK
#define VFLUIDLOOP ILOOP JVLOOP KLOOP VFLUIDCHECK
#define WFLUIDLOOP ILOOP JLOOP KWLOOP WFLUIDCHECK
#define FLUIDLOOP PLAINLOOP PFLUIDCHECK

// SOLVER LOOPS
#define FLEXLOOP IFLEXLOOP JFLEXLOOP KFLEXLOOP FLEXCHECK

// FNPF LOOPS
#define FKJILOOP FKLOOP JLOOP ILOOP
#define FLOOP ILOOP JLOOP FKLOOP FPCHECK
#define FBASELOOP ILOOP JLOOP FKLOOP
#define FILOOP4 ILOOP JLOOP ETALOC PFLUIDCHECK
#define FFILOOP4 ILOOP JLOOP FETALOC FPCHECK

// FORCE LOOP
#define NDBASELOOP ITPLOOP JTPLOOP KTPLOOP

#define NETLOOP for (int n=0; n<p->net_count; ++n)

#define NLOOP1 for(n=p->sizeM1[0]; n<p->sizeM1[1]; ++n)
#define NLOOP2 for(n=p->sizeM2[0]; n<p->sizeM2[1]; ++n)
#define NLOOP3 for(n=p->sizeM3[0]; n<p->sizeM3[1]; ++n)
#define NLOOP4 for(n=p->sizeM4[0]; n<p->sizeM4[1]; ++n)
#define NLOOP4A for(n=p->sizeM4a[0]; n<p->sizeM4a[1]; ++n)
#define NLOOP6 for(n=p->sizeM6[0]; n<p->sizeM6[1]; ++n)
#define NLOOP9 for(n=p->sizeM9[0]; n<p->sizeM9[1]; ++n)
#define NLOOP for(n=sizeM[0]; n<sizeM[1]; ++n)
#define VECLOOP for(n=0; n<p->veclength; ++n)

//MAX, MIN, SIGN
#define MAX(aAa,bBb) ((aAa)>(bBb)?(aAa):(bBb))
#define MIN(aAa,bBb) ((aAa)<(bBb)?(aAa):(bBb))
#define SIGN(aAa) ((aAa)>= 0.0 ? 1.0 : -1.0)

#define MAX3(aAa,bBb,cCc) (((aAa)>(bBb)?(aAa):(bBb))>cCc?((aAa)>(bBb)?(aAa):(bBb)):cCc)
#define MIN3(aAa,bBb,cCc) (((aAa)<(bBb)?(aAa):(bBb))<cCc?((aAa)<(bBb)?(aAa):(bBb)):cCc)

//GCB
#define GCB1 for(n=0;n<p->gcb1_count;++n)
#define GCB1CHECK if(p->gcb1[n][3]>0)
#define GC1LOOP GCB1 GCB1CHECK

#define QGCB1 for(q=0;q<p->gcb1_count;++q)
#define QGCB1CHECK if(p->gcb1[q][3]>0)
#define QGC1LOOP QGCB1 QGCB1CHECK

#define QQGCB1 for(qq=0;qq<p->gcb1_count;++qq)
#define QQGCB1CHECK if(p->gcb1[qq][3]>0)
#define QQGC1LOOP QQGCB1 QQGCB1CHECK

#define GCB2 for(n=0;n<p->gcb2_count;++n)
#define GCB2CHECK if(p->gcb2[n][3]>0)
#define GC2LOOP GCB2 GCB2CHECK

#define QGCB2 for(q=0;q<p->gcb2_count;++q)
#define QGCB2CHECK if(p->gcb2[q][3]>0)
#define QGC2LOOP QGCB2 QGCB2CHECK

#define QQGCB2 for(qq=0;qq<p->gcb2_count;++qq)
#define QQGCB2CHECK if(p->gcb2[qq][3]>0)
#define QQGC2LOOP QQGCB2 QQGCB2CHECK

#define GCB3 for(n=0;n<p->gcb3_count;++n)
#define GCB3CHECK if(p->gcb3[n][3]>0)
#define GC3LOOP GCB3 GCB3CHECK

#define QGCB3 for(q=0;q<p->gcb3_count;++q)
#define QGCB3CHECK if(p->gcb3[q][3]>0)
#define QGC3LOOP QGCB3 QGCB3CHECK

#define QQGCB3 for(qq=0;qq<p->gcb3_count;++qq)
#define QQGCB3CHECK if(p->gcb3[qq][3]>0)
#define QQGC3LOOP QQGCB3 QQGCB3CHECK

#define GCB4 for(n=0;n<p->gcb4_count;++n)
#define GCB4CHECK if(p->gcb4[n][3]>0)
#define GC4LOOP GCB4 GCB4CHECK

#define QGCB4 for(q=0;q<p->gcb4_count;++q)
#define QGCB4CHECK if(p->gcb4[q][3]>0)
#define QGC4LOOP QGCB4 QGCB4CHECK

#define QQGCB4 for(qq=0;qq<p->gcb4_count;++qq)
#define QQGCB4CHECK if(p->gcb4[qq][3]>0)
#define QQGC4LOOP QQGCB4 QQGCB4CHECK

//df
#define QGCDF1 for(q=0;q<p->gcdf1_count;++q)
#define QGCDF1CHECK if(p->gcdf1[q][3]>0)
#define QGCDF1LOOP QGCDF1 QGCDF1CHECK

#define GCDF1 for(n=0;n<p->gcdf1_count;++n)
#define GCDF1CHECK if(p->gcdf1[n][3]>0)
#define GCDF1LOOP GCDF1 GCDF1CHECK

#define QGCDF2 for(q=0;q<p->gcdf2_count;++q)
#define QGCDF2CHECK if(p->gcdf2[q][3]>0)
#define QGCDF2LOOP QGCDF2 QGCDF2CHECK

#define GCDF2 for(n=0;n<p->gcdf2_count;++n)
#define GCDF2CHECK if(p->gcdf2[n][3]>0)
#define GCDF2LOOP GCDF2 GCDF2CHECK

#define QGCDF3 for(q=0;q<p->gcdf3_count;++q)
#define QGCDF3CHECK if(p->gcdf3[q][3]>0)
#define QGCDF3LOOP QGCDF3 QGCDF3CHECK

#define GCDF3 for(n=0;n<p->gcdf3_count;++n)
#define GCDF3CHECK if(p->gcdf3[n][3]>0)
#define GCDF3LOOP GCDF3 GCDF3CHECK

#define QGCDF4 for(q=0;q<p->gcdf4_count;++q)
#define QGCDF4CHECK if(p->gcdf4[q][3]>0)
#define QGCDF4LOOP QGCDF4 QGCDF4CHECK

#define GCDF4 for(n=0;n<p->gcdf4_count;++n)
#define GCDF4CHECK if(p->gcdf4[n][3]>0)
#define GCDF4LOOP GCDF4 GCDF4CHECK

#define GC4A for(n=0;n<p->gcb4a_count;++n)
#define GCB4ACHECK if(p->gcb4a[n][3]>0)
#define GC4ALOOP GC4A GCB4ACHECK

#define QGC4A for(q=0;q<p->gcb4a_count;++q)
#define QGCB4ACHECK if(p->gcb4a[q][3]>0)
#define QGC4ALOOP QGC4A QGCB4ACHECK

#define QQGC4A for(qq=0;qq<p->gcb4a_count;++qq)
#define QQGCB4ACHECK if(p->gcb4a[qq][3]>0)
#define QQGC4ALOOP QQGC4A QQGCB4ACHECK

#define GC6LOOP   for(n=0;n<p->gcb_fix;++n)
#define QGC6LOOP  for(q=0;q<p->gcb_fix;++q)
#define QQGC6LOOP for(qq=0;qq<p->gcb_fix;++qq)
#define GGC6LOOP  for(g=0;g<p->gcb_fix;++g)

#endif
