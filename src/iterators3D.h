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

#ifndef ITERATORS3_H_
#define ITERATORS3_H_

#define IJK  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin

#define Im1JK  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip1JK  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJm1K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin
#define IJp1K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin
#define IJKm1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1
#define IJKp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1

#define Ip1Jp1K  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin

#define Ip1JKp1  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1
#define Ip1JKm1  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1
#define Im1JKp1  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1
#define Im1JKm1  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1

#define IJp1Kp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1
#define IJp1Km1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin-1
#define IJm1Kp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+1
#define IJm1Km1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1

#define Im1Jm1Km1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1
#define Ip1Jm1Km1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1
#define Ip1Jp1Km1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin-1
#define Im1Jp1Km1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin-1

#define Im1Jm1Kp1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+1
#define Ip1Jm1Kp1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+1
#define Ip1Jp1Kp1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1
#define Im1Jp1Kp1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1

#define Ip1Jp2K  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+2)*p->kmax + k-p->kmin

#define Im1JKp2  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2
#define Ip1JKp2  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2

#define IJm1Kp2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+2
#define IJp1Kp2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+2


#define Im2JK  (i-p->imin-2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip2JK  (i-p->imin+2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJm2K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-2)*p->kmax + k-p->kmin
#define IJp2K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+2)*p->kmax + k-p->kmin
#define IJKm2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-2
#define IJKp2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2

#define Im3JK  (i-p->imin-3)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip3JK  (i-p->imin+3)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJm3K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-3)*p->kmax + k-p->kmin
#define IJp3K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+3)*p->kmax + k-p->kmin
#define IJKm3  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-3
#define IJKp3  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+3

#define Im4JK  (i-p->imin-4)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip4JK  (i-p->imin+4)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJKm4  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-4
#define IJKm5  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-5

#define Im5JK  (i-p->imin-5)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip5JK  (i-p->imin+5)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin

#define Im1Jm1K  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin
#define Im1Jp1K  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin
#define Ip1Jm1K  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin


#define FIJK  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin

#define FIm1JK  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIp1JK  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIJm1K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin
#define FIJp1K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin
#define FIJKm1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-1
#define FIJKp1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+1

#define FIm2JK  (i-p->imin-2)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIp2JK  (i-p->imin+2)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIJm2K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-2)*p->kmaxF + k-p->kmin
#define FIJp2K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+2)*p->kmaxF + k-p->kmin
#define FIJKm2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-2
#define FIJKp2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+2

#define FIm3JK  (i-p->imin-3)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIp3JK  (i-p->imin+3)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIJm3K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-3)*p->kmaxF + k-p->kmin
#define FIJp3K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+3)*p->kmaxF + k-p->kmin
#define FIJKm3  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-3
#define FIJKp3  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+3
#define FIJKm4  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-4
#define FIJKp4  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+4
#define FIJKm5  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-5
#define FIJKm6  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-6

#define FIp1JKp1  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+1
#define FIm1JKp1  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+1
#define FIp1JKm1  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-1
#define FIm1JKm1  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-1

#define FIJp1Kp1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin+1
#define FIJm1Kp1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin+1
#define FIJp1Km1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin-1
#define FIJm1Km1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin-1

#define FIp1Jp1K  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin
#define FIm1Jp1K  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin
#define FIp1Jm1K  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin
#define FIm1Jm1K  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin

#define FIp1Jp1Kp1  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin+1
#define FIm1Jm1Kp1  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin+1

#define FIp1Jp2K  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin+2)*p->kmaxF + k-p->kmin

#define FIm1JKp2  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+2
#define FIp1JKp2  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+2

#define FIJm1Kp2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin+2
#define FIJp1Kp2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin+2


#ifdef NHF_KCHECK
#include"nhf_kcheck.h"
#undef IJK
#undef Im1JK
#undef Ip1JK
#undef IJm1K
#undef IJp1K
#undef IJKm1
#undef IJKp1
#undef Ip1Jp1K
#undef Ip1JKp1
#undef Ip1JKm1
#undef Im1JKp1
#undef Im1JKm1
#undef IJp1Kp1
#undef IJp1Km1
#undef IJm1Kp1
#undef IJm1Km1
#undef Im1Jm1Km1
#undef Ip1Jm1Km1
#undef Ip1Jp1Km1
#undef Im1Jp1Km1
#undef Im1Jm1Kp1
#undef Ip1Jm1Kp1
#undef Ip1Jp1Kp1
#undef Im1Jp1Kp1
#undef Ip1Jp2K
#undef Im1JKp2
#undef Ip1JKp2
#undef IJm1Kp2
#undef IJp1Kp2
#undef Im2JK
#undef Ip2JK
#undef IJm2K
#undef IJp2K
#undef IJKm2
#undef IJKp2
#undef Im3JK
#undef Ip3JK
#undef IJm3K
#undef IJp3K
#undef IJKm3
#undef IJKp3
#undef Im4JK
#undef Ip4JK
#undef IJKm4
#undef IJKm5
#undef Im5JK
#undef Ip5JK
#undef Im1Jm1K
#undef Im1Jp1K
#undef Ip1Jm1K
#undef FIJK
#undef FIm1JK
#undef FIp1JK
#undef FIJm1K
#undef FIJp1K
#undef FIJKm1
#undef FIJKp1
#undef FIm2JK
#undef FIp2JK
#undef FIJm2K
#undef FIJp2K
#undef FIJKm2
#undef FIJKp2
#undef FIm3JK
#undef FIp3JK
#undef FIJm3K
#undef FIJp3K
#undef FIJKm3
#undef FIJKp3
#undef FIJKm4
#undef FIJKp4
#undef FIJKm5
#undef FIJKm6
#undef FIp1JKp1
#undef FIm1JKp1
#undef FIp1JKm1
#undef FIm1JKm1
#undef FIJp1Kp1
#undef FIJm1Kp1
#undef FIJp1Km1
#undef FIJm1Km1
#undef FIp1Jp1K
#undef FIm1Jp1K
#undef FIp1Jm1K
#undef FIm1Jm1K
#undef FIp1Jp1Kp1
#undef FIm1Jm1Kp1
#undef FIp1Jp2K
#undef FIm1JKp2
#undef FIp1JKp2
#undef FIJm1Kp2
#undef FIJp1Kp2
#define IJK (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Im1JK (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Ip1JK (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJm1K (i-p->imin)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJp1K (i-p->imin)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJKm1 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define IJKp1 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define Ip1Jp1K (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Ip1JKp1 (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define Ip1JKm1 (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define Im1JKp1 (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define Im1JKm1 (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define IJp1Kp1 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define IJp1Km1 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define IJm1Kp1 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define IJm1Km1 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define Im1Jm1Km1 (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define Ip1Jm1Km1 (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define Ip1Jp1Km1 (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define Im1Jp1Km1 (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-1,0,__FILE__,__LINE__)-p->kmin
#define Im1Jm1Kp1 (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define Ip1Jm1Kp1 (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define Ip1Jp1Kp1 (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define Im1Jp1Kp1 (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+1,0,__FILE__,__LINE__)-p->kmin
#define Ip1Jp2K (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j+2,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Im1JKp2 (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+2,0,__FILE__,__LINE__)-p->kmin
#define Ip1JKp2 (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+2,0,__FILE__,__LINE__)-p->kmin
#define IJm1Kp2 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+2,0,__FILE__,__LINE__)-p->kmin
#define IJp1Kp2 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+2,0,__FILE__,__LINE__)-p->kmin
#define Im2JK (i-p->imin-2)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Ip2JK (i-p->imin+2)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJm2K (i-p->imin)*p->jmax*p->kmax + (jchk(p,j-2,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJp2K (i-p->imin)*p->jmax*p->kmax + (jchk(p,j+2,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJKm2 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-2,0,__FILE__,__LINE__)-p->kmin
#define IJKp2 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+2,0,__FILE__,__LINE__)-p->kmin
#define Im3JK (i-p->imin-3)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Ip3JK (i-p->imin+3)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJm3K (i-p->imin)*p->jmax*p->kmax + (jchk(p,j-3,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJp3K (i-p->imin)*p->jmax*p->kmax + (jchk(p,j+3,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJKm3 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-3,0,__FILE__,__LINE__)-p->kmin
#define IJKp3 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k+3,0,__FILE__,__LINE__)-p->kmin
#define Im4JK (i-p->imin-4)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Ip4JK (i-p->imin+4)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define IJKm4 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-4,0,__FILE__,__LINE__)-p->kmin
#define IJKm5 (i-p->imin)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k-5,0,__FILE__,__LINE__)-p->kmin
#define Im5JK (i-p->imin-5)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Ip5JK (i-p->imin+5)*p->jmax*p->kmax + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Im1Jm1K (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Im1Jp1K (i-p->imin-1)*p->jmax*p->kmax + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define Ip1Jm1K (i-p->imin+1)*p->jmax*p->kmax + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmax + kchk(p,k,0,__FILE__,__LINE__)-p->kmin
#define FIJK (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIm1JK (i-p->imin-1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIp1JK (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJm1K (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJp1K (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJKm1 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-1,1,__FILE__,__LINE__)-p->kmin
#define FIJKp1 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+1,1,__FILE__,__LINE__)-p->kmin
#define FIm2JK (i-p->imin-2)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIp2JK (i-p->imin+2)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJm2K (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j-2,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJp2K (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j+2,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJKm2 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-2,1,__FILE__,__LINE__)-p->kmin
#define FIJKp2 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+2,1,__FILE__,__LINE__)-p->kmin
#define FIm3JK (i-p->imin-3)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIp3JK (i-p->imin+3)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJm3K (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j-3,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJp3K (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j+3,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIJKm3 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-3,1,__FILE__,__LINE__)-p->kmin
#define FIJKp3 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+3,1,__FILE__,__LINE__)-p->kmin
#define FIJKm4 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-4,1,__FILE__,__LINE__)-p->kmin
#define FIJKp4 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+4,1,__FILE__,__LINE__)-p->kmin
#define FIJKm5 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-5,1,__FILE__,__LINE__)-p->kmin
#define FIJKm6 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-6,1,__FILE__,__LINE__)-p->kmin
#define FIp1JKp1 (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+1,1,__FILE__,__LINE__)-p->kmin
#define FIm1JKp1 (i-p->imin-1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+1,1,__FILE__,__LINE__)-p->kmin
#define FIp1JKm1 (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-1,1,__FILE__,__LINE__)-p->kmin
#define FIm1JKm1 (i-p->imin-1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-1,1,__FILE__,__LINE__)-p->kmin
#define FIJp1Kp1 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+1,1,__FILE__,__LINE__)-p->kmin
#define FIJm1Kp1 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+1,1,__FILE__,__LINE__)-p->kmin
#define FIJp1Km1 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-1,1,__FILE__,__LINE__)-p->kmin
#define FIJm1Km1 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k-1,1,__FILE__,__LINE__)-p->kmin
#define FIp1Jp1K (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIm1Jp1K (i-p->imin-1)*p->jmax*p->kmaxF + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIp1Jm1K (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIm1Jm1K (i-p->imin-1)*p->jmax*p->kmaxF + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIp1Jp1Kp1 (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+1,1,__FILE__,__LINE__)-p->kmin
#define FIm1Jm1Kp1 (i-p->imin-1)*p->jmax*p->kmaxF + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+1,1,__FILE__,__LINE__)-p->kmin
#define FIp1Jp2K (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j+2,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k,1,__FILE__,__LINE__)-p->kmin
#define FIm1JKp2 (i-p->imin-1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+2,1,__FILE__,__LINE__)-p->kmin
#define FIp1JKp2 (i-p->imin+1)*p->jmax*p->kmaxF + (jchk(p,j,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+2,1,__FILE__,__LINE__)-p->kmin
#define FIJm1Kp2 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j-1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+2,1,__FILE__,__LINE__)-p->kmin
#define FIJp1Kp2 (i-p->imin)*p->jmax*p->kmaxF + (jchk(p,j+1,__FILE__,__LINE__)-p->jmin)*p->kmaxF + kchk(p,k+2,1,__FILE__,__LINE__)-p->kmin
#endif
#endif
