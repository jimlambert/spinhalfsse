#ifndef MEAS_UTLS_H
#define MEAS_UTLS_H

#include <vector>
#include <complex>

#include "configuration.h"
#include "parameters.h"
#include "initfunctions.h"

namespace spinhalfsse
{

typedef double meas_t;
typedef std::vector<meas_t> measarr_t;

typedef std::complex<double> cmeas_t;
typedef std::vector<cmeas_t> cmeasarr_t;

typedef std::complex<double> phase_t;
typedef std::vector<phase_t> phasearr_t;
typedef std::vector<phasearr_t> phasearrs_t;

typedef std::vector<int> posarr_t;
typedef std::vector<posarr_t> posarrs_t;


/*
 *  Accepts an initial site, start, and a final site, end.
 *
 *  Returns a posarr_t containing integers [start,end-1].
 */
posarr_t measutls_initposarr(const int&, const int&);


/*
 *  Accepts an  integer L corresponding to the length of the system,
 *  an integer k corresponding to a position in the BZ and 
 *  a list of positions. 
 *
 *  Returns the phases corresponding to those positions. 
 */
phasearr_t measutls_initphasearr(const int&, const int&, const posarr_t&);


/*
 *  Accepts a configuration, an array of position indices, and an array 
 *  of phases.
 *
 *  Returns the average of the sz components at the spin indices multiplied by the
 *  phases at the current propagation index. 
 *  For example, if the phases are trivial this is the average
 *  magnetization. If the phases are characterized by k=\pi the answer is the
 *  average staggered magnetization. 
 */
cmeas_t measutls_szmag
(
  const Configuration&,
  const posarr_t&, 
  const phasearr_t&
);


/*
 *  The same as the above function but only returns the real part.
 */
meas_t measutls_realszmag
(
  const Configuration&,
  const posarr_t&, 
  const phasearr_t&
);


/*
 *  Measures the real value of the Dz operator at the given momenta
 */
meas_t measutls_realdzmag
(
  const Configuration&,
  const posarr_t&, 
  const phasearr_t&
);


struct cmeasset
{
  cmeas_t avemag{0.0,0.0};
  cmeas_t varmag{0.0,0.0};  
  cmeas_t suscep{0.0,0.0};
};


struct measset
{
  meas_t avemag = 0.0;
  meas_t varmag = 0.0;
  meas_t suscep = 0.0;
};


typedef std::vector<cmeasset> cmeassets_t;
typedef std::vector<measset> meassets_t;


cmeassets_t measutls_szmagconf
(
  Configuration&,
  const Parameters&,
  const posarr_t&,
  const phasearrs_t&  
);


meassets_t measutls_realszmagconf
(
  Configuration&,
  const Parameters&,
  const posarr_t&,
  const phasearrs_t&  
);


meassets_t measutls_realdzmagconf
(
  Configuration&,
  const Parameters&,
  const posarr_t&,
  const phasearrs_t&  
);


meassets_t measutls_realszmagconf
(
  Configuration&,
  const Parameters&,
  const posarrs_t&,
  const phasearrs_t&  
);


meassets_t measutls_realdzmagconf
(
  Configuration&,
  const Parameters&,
  const posarrs_t&,
  const phasearrs_t&  
);


/*
 *  Measures the different between the operators,
 *  N+ = S+S- and N-= S-S+ divided by the system size
 *  L on bonds of type <bondType>
 */
int currDiff
(
  Configuration&,
  const int&
);


}

#endif // MEAS_UTLS_H
