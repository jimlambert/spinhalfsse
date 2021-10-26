#ifndef INIT_FUNCTIONS_H
#define INIT_FUNCTIONS_H

#include <array>

#include "parameters.h"
#include "configuration.h"
#include "lattice.h"
#include "spinHalfVertices.h" // definition of NVERTS

namespace spinhalfsse
{

typedef std::array<int, NVERTS> sgnlst_t;
typedef std::array<double,NVERTS> wgtlst_t;
typedef std::vector<wgtlst_t> wgtlsts_t;

typedef std::array<std::array<int,5>, 4> trnmat_t;
typedef std::array<std::array<double,5>, 4> prbmat_t;

typedef std::array<trnmat_t,NVERTS> trnlst_t;
typedef std::array<prbmat_t,NVERTS> prblst_t;

typedef std::array<trnlst_t,2> trnlstset_t;
typedef std::array<prblst_t,2> prblstset_t;

typedef std::vector<trnlstset_t> trnlstsets_t;
typedef std::vector<prblstset_t> prblstsets_t;

typedef std::vector<int> linklst_t;


typedef
struct LinksFrsts
{
  linklst_t linklst;
  std::vector<int> vfrst;
} lnkfrst_t;


/*
 *  Returns a 4x5 matrix with all zero entries.
 */
prbmat_t zeroprbmat();


/*
 *  Accepts model parameters.
 *
 *  Returns list of weights for the vertices in VERTEXTBL (see spin1vertices.h).
 */
wgtlst_t initwgts(const Parameters&);


/*
 *  Accepts wgtlst_t.
 *
 *  Returns list of sgns for the weights.
 */
sgnlst_t initsgns(const wgtlst_t&);


/*
 *  Accepts magnitude of spin flip at loop entrance.
 *
 *  Returns list of transformed vertices
 */
trnlst_t inittrnlst(const int&);


/*
 *  Accepts a transfer list and a list of weights.
 *
 *  Returns probability for each exit leg given an entrace leg. 
 */
prblst_t initprblst(const trnlst_t&, const wgtlst_t&);


/*
 *  Accepts a list of probability matrices.
 *
 *  Returns a list of cumulative probabilities.
 */
prblst_t initcpdlst(const prblst_t&);


/*
 *  Accepts a configuration.
 *
 *  Returns the linked list corresponding to the operators in the configuration.
 */
lnkfrst_t getlinklst(Configuration&);

}


#endif // INIT_FUNCTIONS_H
