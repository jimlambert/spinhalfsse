#ifndef GRAPH_INITS_H
#define GRAPH_INITS_H

#include "lattice.h"

namespace LatticeUtls
{

/*
 *  Accepts number of sites and a bool which is true for periodic boundary
 *  conditions and false for open.
 *
 *  Returns graph corresponding to nearest neighbour bonds 
 */
graph_t initchain(const int&, const bool&);


/*
 *  Accepts x and y dimensions of square lattice and boundary conditions in x
 *  and y.
 *
 *  Returns graph with single edge list (all edges treated equally).
 */
graph_t initisosqlatt
(
  const int& lx, 
  const int& ly, 
  const bool bcx, 
  const bool bcy
);
 

/*
 *  Accepts x and y dimensions of square lattice and boundary conditions in x
 *  and y.
 *
 *  Returns graph with x and y edge lists for anisotropic model.
 */
graph_t initsqlatt
(
  const int& lx, 
  const int& ly, 
  const bool bcx, 
  const bool bcy
);


/*
 *  Accepts number of unit cells in x and y, and boundary conditions for btoh
 *  (True for periodic)
 *
 *  Returns graph with single edge list for all bonds.
 */
graph_t inithoneycomb
(
  const int& nx, 
  const int& ny, 
  const bool bcx, 
  const bool bcy
);


/*
 *  Accepts x and y dimensions of square lattice and boundary conditions in x
 *  and y.
 *
 *  Returns graph with single edge list (all edges treated equally), and the
 *  dangling z-spins for the Kondo necklace model.
 */
graph_t initIsoSqKondoNecklace
(
  const int& nx,
  const int& ny,
  const bool bcx,
  const bool bcy 
);


/*
 *  Accepts x and y dimensions of square lattice and boundary conditions in x
 *  and y.
 *
 *  Graph has three edge lists for x bonds, y bonds, and dangling bonds
 */
graph_t initIsoSqKondoNecklace
(
  const int& nx,
  const int& ny,
  const bool bcx,
  const bool bcy 
);


graph_t initSqKondoNecklace
(
  const int& nx,
  const int& ny,
  const bool bcx,
  const bool bcy 
);


}

#endif // GRAPH_INITS_H
