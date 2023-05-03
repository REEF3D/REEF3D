/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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

#include"vrans_net.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_net::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    // Disable common porosity
	ALOOP
	{
        a->porosity(i,j,k) = 1.0;
	}
	pgc->start4(p,a->porosity,1);
    p->B260 = 0.0;
    
    // Parallelisation ini
	p->Darray(xstart, p->mpi_size);
	p->Darray(xend, p->mpi_size);
	p->Darray(ystart, p->mpi_size);
	p->Darray(yend, p->mpi_size);
	p->Darray(zstart, p->mpi_size);
	p->Darray(zend, p->mpi_size);
	
	xstart[p->mpirank] = p->originx;
	ystart[p->mpirank] = p->originy;
	zstart[p->mpirank] = p->originz;
	xend[p->mpirank] = p->endx;
	yend[p->mpirank] = p->endy;
	zend[p->mpirank] = p->endz;
	
	for (int ii = 0; ii < p->mpi_size; ii++)
	{
		MPI_Bcast(&xstart[ii],1,MPI_DOUBLE,ii,pgc->mpi_comm);
		MPI_Bcast(&xend[ii],1,MPI_DOUBLE,ii,pgc->mpi_comm);
		MPI_Bcast(&ystart[ii],1,MPI_DOUBLE,ii,pgc->mpi_comm);
		MPI_Bcast(&yend[ii],1,MPI_DOUBLE,ii,pgc->mpi_comm);
		MPI_Bcast(&zstart[ii],1,MPI_DOUBLE,ii,pgc->mpi_comm);
		MPI_Bcast(&zend[ii],1,MPI_DOUBLE,ii,pgc->mpi_comm);
    }          
}

