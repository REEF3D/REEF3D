/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#define Im1Jm2  (i-p->imin-1)*p->jmax + (j-p->jmin-2)
#define Ip1Jp2  (i-p->imin+1)*p->jmax + (j-p->jmin+2)
#define Im1Jp2  (i-p->imin-1)*p->jmax + (j-p->jmin+2)
#define Ip1Jm2  (i-p->imin+1)*p->jmax + (j-p->jmin-2)

#define Im2Jm1  (i-p->imin-2)*p->jmax + (j-p->jmin-1)
#define Ip2Jp1  (i-p->imin+2)*p->jmax + (j-p->jmin+1)
#define Im2Jp1  (i-p->imin-2)*p->jmax + (j-p->jmin+1)
#define Ip2Jm1  (i-p->imin+2)*p->jmax + (j-p->jmin-1)

#define Im2Jm2  (i-p->imin-2)*p->jmax + (j-p->jmin-2)
#define Ip2Jp2  (i-p->imin+2)*p->jmax + (j-p->jmin+2)
#define Im2Jp2  (i-p->imin-2)*p->jmax + (j-p->jmin+2)
#define Ip2Jm2  (i-p->imin+2)*p->jmax + (j-p->jmin-2)


#define Im2J  (i-p->imin-2)*p->jmax + (j-p->jmin)
#define Ip2J  (i-p->imin+2)*p->jmax + (j-p->jmin)
#define IJm2  (i-p->imin)*p->jmax + (j-p->jmin-2)
#define IJp2  (i-p->imin)*p->jmax + (j-p->jmin+2)

#define Im3J  (i-p->imin-3)*p->jmax + (j-p->jmin)
#define Ip3J  (i-p->imin+3)*p->jmax + (j-p->jmin)
#define IJm3  (i-p->imin)*p->jmax + (j-p->jmin-3)
#define IJp3  (i-p->imin)*p->jmax + (j-p->jmin+3)

