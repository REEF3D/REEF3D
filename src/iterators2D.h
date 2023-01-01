/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is frb->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Frb->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sb->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#define IJ  (i-p->imin)*p->jmax + (j-p->jmin) 

#define Im1J  (i-p->imin-1)*p->jmax + (j-p->jmin)
#define Ip1J  (i-p->imin+1)*p->jmax + (j-p->jmin)
#define IJm1  (i-p->imin)*p->jmax + (j-p->jmin-1)
#define IJp1  (i-p->imin)*p->jmax + (j-p->jmin+1)

#define Im1Jm1  (i-p->imin-1)*p->jmax + (j-p->jmin-1)
#define Ip1Jp1  (i-p->imin+1)*p->jmax + (j-p->jmin+1)

#define Im1Jp1  (i-p->imin-1)*p->jmax + (j-p->jmin+1)
#define Ip1Jm1  (i-p->imin+1)*p->jmax + (j-p->jmin-1)


#define Im2J  (i-p->imin-2)*p->jmax + (j-p->jmin)
#define Ip2J  (i-p->imin+2)*p->jmax + (j-p->jmin)
#define IJm2  (i-p->imin)*p->jmax + (j-p->jmin-2)
#define IJp2  (i-p->imin)*p->jmax + (j-p->jmin+2)

#define Im3J  (i-p->imin-3)*p->jmax + (j-p->jmin)
#define Ip3J  (i-p->imin+3)*p->jmax + (j-p->jmin)
#define IJm3  (i-p->imin)*p->jmax + (j-p->jmin-3)
#define IJp3  (i-p->imin)*p->jmax + (j-p->jmin+3)


#define I_J_1 n
#define I_J_2 n
#define I_J_3 n
#define I_J_4 n
#define I_J_4a n
#define I_J n

// ----------------------------------------------

#define Ip1_J C.n[n]
#define Im1_J C.s[n]
#define I_Jp1 C.w[n]
#define I_Jm1 C.e[n]

#define Ip2_J C.n[C.n[n]]
#define Im2_J C.s[C.s[n]]
#define I_Jp2 C.w[C.w[n]]
#define I_Jm2 C.e[C.e[n]]

#define Ip3_J C.n[C.n[C.n[n]]]
#define Im3_J C.s[C.s[C.s[n]]]
#define I_Jp3 C.w[C.w[C.w[n]]]
#define I_Jm3 C.e[C.e[C.e[n]]]
// 1 ----------------------------------------------

#define Ip1_J_1 bb->C1.n[n]
#define Im1_J_1 bb->C1.s[n]
#define I_Jp1_1 bb->C1.w[n]
#define I_Jm1_1 bb->C1.e[n]

#define Ip2_J_1 bb->C1.n[bb->C1.n[n]]
#define Im2_J_1 bb->C1.s[bb->C1.s[n]]
#define I_Jp2_1 bb->C1.w[bb->C1.w[n]]
#define I_Jm2_1 bb->C1.e[bb->C1.e[n]]

#define Ip3_J_1 bb->C1.n[bb->C1.n[bb->C1.n[n]]]
#define Im3_J_1 bb->C1.s[bb->C1.s[bb->C1.s[n]]]
#define I_Jp3_1 bb->C1.w[bb->C1.w[bb->C1.w[n]]]
#define I_Jm3_1 bb->C1.e[bb->C1.e[bb->C1.e[n]]]

#define Ip1_Jp1_1 bb->C1.n[bb->C1.w[n]]
#define Ip1_Jp1_1 bb->C1.n[bb->C1.t[n]]
#define Ip1_Jm1_1 bb->C1.n[bb->C1.e[n]]
#define Ip1_Jm1_1 bb->C1.n[bb->C1.b[n]]
#define Im1_Jp1_1 bb->C1.s[bb->C1.w[n]]
#define Im1_Jp1_1 bb->C1.s[bb->C1.t[n]]
#define Im1_Jm1_1 bb->C1.s[bb->C1.e[n]]
#define Im1_Jm1_1 bb->C1.s[bb->C1.b[n]]
#define I_Jp1p1_1 bb->C1.w[bb->C1.t[n]]
#define I_Jp1m1_1 bb->C1.w[bb->C1.b[n]]
#define I_Jm1p1_1 bb->C1.e[bb->C1.t[n]]
#define I_Jm1m1_1 bb->C1.e[bb->C1.b[n]]

// 2 ----------------------------------------------

#define Ip1_J_2 bb->C2.n[n]
#define Im1_J_2 bb->C2.s[n]
#define I_Jp1_2 bb->C2.w[n]
#define I_Jm1_2 bb->C2.e[n]

#define Ip2_J_2 bb->C2.n[bb->C2.n[n]]
#define Im2_J_2 bb->C2.s[bb->C2.s[n]]
#define I_Jp2_2 bb->C2.w[bb->C2.w[n]]
#define I_Jm2_2 bb->C2.e[bb->C2.e[n]]

#define Ip3_J_2 bb->C2.n[bb->C2.n[bb->C2.n[n]]]
#define Im3_J_2 bb->C2.s[bb->C2.s[bb->C2.s[n]]]
#define I_Jp3_2 bb->C2.w[bb->C2.w[bb->C2.w[n]]]
#define I_Jm3_2 bb->C2.e[bb->C2.e[bb->C2.e[n]]]

#define Ip1_Jp1_2 bb->C2.n[bb->C2.w[n]]
#define Ip1_Jp1_2 bb->C2.n[bb->C2.t[n]]
#define Ip1_Jm1_2 bb->C2.n[bb->C2.e[n]]
#define Ip1_Jm1_2 bb->C2.n[bb->C2.b[n]]
#define Im1_Jp1_2 bb->C2.s[bb->C2.w[n]]
#define Im1_Jp1_2 bb->C2.s[bb->C2.t[n]]
#define Im1_Jm1_2 bb->C2.s[bb->C2.e[n]]
#define Im1_Jm1_2 bb->C2.s[bb->C2.b[n]]
#define I_Jp1p1_2 bb->C2.w[bb->C2.t[n]]
#define I_Jp1m1_2 bb->C2.w[bb->C2.b[n]]
#define I_Jm1p1_2 bb->C2.e[bb->C2.t[n]]
#define I_Jm1m1_2 bb->C2.e[bb->C2.b[n]]

// 4 ----------------------------------------------

#define Ip1_J_4 bb->C4.n[n]
#define Im1_J_4 bb->C4.s[n]
#define I_Jp1_4 bb->C4.w[n]
#define I_Jm1_4 bb->C4.e[n]

#define Ip2_J_4 bb->C4.n[bb->C4.n[n]]
#define Im2_J_4 bb->C4.s[bb->C4.s[n]]
#define I_Jp2_4 bb->C4.w[bb->C4.w[n]]
#define I_Jm2_4 bb->C4.e[bb->C4.e[n]]

#define Ip3_J_4 bb->C4.n[bb->C4.n[bb->C4.n[n]]]
#define Im3_J_4 bb->C4.s[bb->C4.s[bb->C4.s[n]]]
#define I_Jp3_4 bb->C4.w[bb->C4.w[bb->C4.w[n]]]
#define I_Jm3_4 bb->C4.e[bb->C4.e[bb->C4.e[n]]]

#define Ip1_Jp1_4 bb->C4.n[bb->C4.w[n]]
#define Ip1_Jp1_4 bb->C4.n[bb->C4.t[n]]
#define Ip1_Jm1_4 bb->C4.n[bb->C4.e[n]]
#define Ip1_Jm1_4 bb->C4.n[bb->C4.b[n]]
#define Im1_Jp1_4 bb->C4.s[bb->C4.w[n]]
#define Im1_Jp1_4 bb->C4.s[bb->C4.t[n]]
#define Im1_Jm1_4 bb->C4.s[bb->C4.e[n]]
#define Im1_Jm1_4 bb->C4.s[bb->C4.b[n]]
#define I_Jp1p1_4 bb->C4.w[bb->C4.t[n]]
#define I_Jp1m1_4 bb->C4.w[bb->C4.b[n]]
#define I_Jm1p1_4 bb->C4.e[bb->C4.t[n]]
#define I_Jm1m1_4 bb->C4.e[bb->C4.b[n]]
