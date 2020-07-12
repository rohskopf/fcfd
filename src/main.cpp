#include <stdlib.h>
#include <iostream>
#include "geo.h"
#include "mpi.h"

using namespace GEO_NS;

int main(int argc, char **argv)
{

  /* Initialize MPI */
  MPI_Init(&argc,&argv);

  /* Begin a GEO instance */
  GEO *geo = new GEO(argc, argv);

  /* Delete the memory */
  delete geo;

  /* Close MPI */
  int MPI_Comm_free(MPI_Comm *comm);
  MPI_Finalize();

  return EXIT_SUCCESS;
}

