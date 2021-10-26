#ifndef IO_FUNCTIONS_H
#define IO_FUNCTIONS_H

#include <string>

#include "parameters.h"
#include "spinHalfVertices.h"
#include "initfunctions.h"
#include "measutls.h"

namespace spinhalfsse
{

/*
 *  I/O functions for common structures
 */
void write(const trnmat_t&);

void write(const trnlst_t&);

void write(const prbmat_t&);

void write(const prblst_t&);

void write(const wgtlst_t&);

void write(const sgnlst_t&);

void write(const LatticeUtls::bond_t&);

void write(const LatticeUtls::bondlst_t&);

void write(const LatticeUtls::graph_t&);

void write(const linklst_t&);

void write(const phasearr_t&);

void writevertex(const int&);


/*
 *  ARGS
 *  ----
 *
 *    const std::string&   -   process name
 *    const int&           -   max step
 *    const int&           -   current step
 *
 *  Displays progress of <process name> as percentage. When process reachers
 *  termination condition(<max step> == <current step>) prints appropriate
 *  message and indents line. 
 */
void progress(const std::string&, const int&, const int&);


/*
 *  Accepts a parameter file name.
 *
 *  Returns parameters extracted from file
 */
paramfile_t extractparams(char*);


/*
 *  Accepts a file name of the format <name>.<ext>. 
 *
 *  Returns a string of the format <name>.out.
 */ 
std::string convertfilename(char* fname);


} // namespace spin1sse

#endif // IO_FUNCTIONS_H
