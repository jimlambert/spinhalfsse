#ifndef UPDATE_FUNCTIONS_H
#define UPDATE_FUNCTIONS_H

#include "randomhelpers.h"
#include "configuration.h"
#include "spinHalfVertices.h"
#include "initfunctions.h"

namespace spinhalfsse
{


typedef 
struct LoopCoord
{
  int vtindex;
  int leg;
} loopcrd_t;

typedef std::vector<int> freesplst_t;

/*
 *  Accepts a configuration with an updated list of operators and a list of
 *  spins unattached to any operator.
 *
 *  The spins that are unattached to any operator are randomly changed. The
 *  spins that are connected to other operators are updated according to the new
 *  operator list structure.
 */
void spinupdt(Configuration&, const freesplst_t&);


void xordupdt(Configuration&);


/*
 *  Performs diagonal update on multiple types of bond weight list. 
 */
void diagupdt
(
  Configuration&, 
  LatticeUtls::Lattice&, 
  Parameters&, 
  wgtlsts_t&,
  RandomHelpers&
);


void oneloopupdt
(
  Configuration&,
  sgnlst_t&,
  trnlstset_t&,
  prblstset_t&,
  RandomHelpers& 
);


void oneloopupdt
(
  Configuration&,
  trnlstset_t&,
  prblstsets_t&,
  RandomHelpers& 
);



}

#endif // UPDATE_FUNCTONS_H
