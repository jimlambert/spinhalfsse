#include <assert.h>
#include <iostream>
#include <cmath>
#include <complex>

#include "measutls.h"

using namespace std;

namespace spinhalfsse
{


posarr_t measutls_initposarr(const int& start, const int& end)
{
  posarr_t posarr;

  for(int i=start; i<end; i++)
    posarr.push_back(i);

  return posarr;
}


phasearr_t measutls_initphasearr
(
  const int& L, 
  const int& k, 
  const posarr_t& posarr
)
{
  phasearr_t phasearr;

  for(size_t i=0; i<posarr.size(); i++)
  {
    phase_t arg = (double)2*1i*(double)k*(double)posarr[i]*M_PI/(double)L;
    phasearr.push_back(exp(arg));
  }

  return phasearr;
}


cmeas_t measutls_szmag
( 
  const Configuration& config,
  const posarr_t& posarr, 
  const phasearr_t& phasearr
)
{
  size_t np = posarr.size();

  // ensure sufficient number of phases have been specified
  assert(np == phasearr.size());

  cmeas_t total(0.0,0.0);
  for(size_t i=0; i<np; i++) 
    total += (double)config.getspin(posarr[i])*phasearr[i];
  
  return total;
}


meas_t measutls_realszmag
( 
  const Configuration& config,
  const posarr_t& posarr, 
  const phasearr_t& phasearr
)
{
  size_t np = posarr.size();

  // ensure sufficient number of phases have been specified
  assert(np == phasearr.size());

  cmeas_t total(0.0,0.0);
  for(size_t i=0; i<np; i++) 
    total += (double)config.getspin(posarr[i])*phasearr[i];
  
  return real(total);
}


meas_t measutls_realdzmag
( 
  const Configuration& config,
  const posarr_t& posarr, 
  const phasearr_t& phasearr
)
{
  size_t np = posarr.size();

  // ensure sufficient number of phases have been specified
  assert(np == phasearr.size());

  cmeas_t total(0.0,0.0);
  for(size_t i=0; i<np; i++) 
  {
    double s = config.getspin(posarr[i]);
    total += s*s*phasearr[i];
  }
  
  return real(total);
}



cmeassets_t measutls_szmagconf
(
  Configuration& config,
  const Parameters& params,
  const posarr_t& posarr,
  const phasearrs_t& phasearrs 
)
{
  int xo=config.xorder();
  size_t nk=phasearrs.size();
  
  cmeassets_t cmeassets;
  cmeassets.resize(phasearrs.size()); 

  cmeasarr_t cmeasarr;
  cmeasarr.resize(nk,cmeas_t{0.0,0.0});
 
  for(int p=0; p<xo; p++)
  {
    for(size_t k=0; k<nk; k++)
    {
      cmeas_t szmag = measutls_szmag(config, posarr, phasearrs[k]);
      
      cmeassets[k].avemag += szmag;
      cmeassets[k].suscep += szmag*szmag;
    }
    config.propagate();  
  }

  for(size_t k=0; k<nk; k++)
  {
    cmeassets[k].avemag = fabs(real(cmeassets[k].avemag)) / (double)xo;
    cmeassets[k].varmag = cmeassets[k].suscep/(double)xo;
    cmeassets[k].suscep = params.bt*
                         (
                            cmeassets[k].suscep/(double)((xo+1)*(xo+1))
                          + cmeasarr[k]/(double)(xo*(xo+1))
                         ); 
  }

  return cmeassets;
}

  
meassets_t measutls_realszmagconf
(
  Configuration& config,
  const Parameters& params,
  const posarr_t& posarr,
  const phasearrs_t& phasearrs 
)
{
  int xo=config.xorder();
  size_t nk=phasearrs.size();
  
  meassets_t meassets;
  meassets.resize(phasearrs.size()); 

  measarr_t measarr;
  measarr.resize(nk,0.0);
 
  for(int p=0; p<xo; p++)
  {
    for(size_t k=0; k<nk; k++)
    {
      meas_t szmag = measutls_realszmag(config, posarr, phasearrs[k]);
      
      meassets[k].avemag += szmag;
      meassets[k].varmag += szmag*szmag;
      meassets[k].suscep += szmag*szmag;
    }
    config.propagate();  
  }

  for(size_t k=0; k<nk; k++)
  {
    meassets[k].suscep = (params.bt / (xo*(xo+1)))*
                         (
                            meassets[k].suscep
                          + (meassets[k].avemag*meassets[k].avemag)
                         ); 
    meassets[k].avemag = fabs(meassets[k].avemag) / (double)xo;
    meassets[k].varmag = meassets[k].varmag / (double)xo;
    meassets[k].suscep = meassets[k].suscep;
  }

  return meassets;
}


meassets_t measutls_realdzmagconf
(
  Configuration& config,
  const Parameters& params,
  const posarr_t& posarr,
  const phasearrs_t& phasearrs 
)
{
  int xo=config.xorder();
  size_t nk=phasearrs.size();
  
  meassets_t meassets;
  meassets.resize(phasearrs.size()); 

  measarr_t measarr;
  measarr.resize(nk,0.0);
 
  for(int p=0; p<xo; p++)
  {
    for(size_t k=0; k<nk; k++)
    {
      meas_t dzmag = measutls_realdzmag(config, posarr, phasearrs[k]);
      
      meassets[k].avemag += dzmag;
      meassets[k].varmag += dzmag*dzmag;
      meassets[k].suscep += dzmag*dzmag;
    }
    config.propagate();  
  }

  for(size_t k=0; k<nk; k++)
  {
    meassets[k].suscep = (params.bt / (xo*(xo+1)))*
                         (
                            meassets[k].suscep
                          + (meassets[k].avemag*meassets[k].avemag)
                         ); 
    meassets[k].avemag = fabs(meassets[k].avemag) / (double)xo;
    meassets[k].varmag = meassets[k].varmag / (double)xo;
    meassets[k].suscep = meassets[k].suscep;
  }

  return meassets;
}


meassets_t measutls_realszmagconf
(
  Configuration& config,
  const Parameters& params,
  const posarrs_t& posarrs,
  const phasearrs_t& phasearrs 
)
{
  int xo=config.xorder();
  size_t nr=posarrs.size();
  size_t nk=phasearrs.size();
  
  meassets_t meassets;
  meassets.resize(nk*nr); 

  measarr_t measarr;
  measarr.resize(nk,0.0);
 
  for(int p=0; p<xo; p++)
  {
    for(size_t k=0; k<nk; k++)
    for(size_t r=0; r<nr; r++)
    {
      size_t index = k*nr + r;
      meas_t szmag = measutls_realszmag(config, posarrs[r], phasearrs[k]);
      
      meassets[index].avemag += szmag;
      meassets[index].varmag += szmag*szmag;
      meassets[index].suscep += szmag*szmag;
    }
    config.propagate();  
  }

  for(size_t k=0; k<nk; k++)
  for(size_t r=0; r<nr; r++)
  {
    size_t index = k*nr + r;
    meassets[index].suscep = (params.bt / (xo*(xo+1)))*
                         (
                            meassets[index].suscep
                          + (meassets[index].avemag*meassets[index].avemag)
                         ); 
    meassets[index].avemag = fabs(meassets[index].avemag) / (double)xo;
    meassets[index].varmag = meassets[index].varmag / (double)xo;
    meassets[index].suscep = meassets[index].suscep;
  }

  return meassets;
}


meassets_t measutls_realdzmagconf
(
  Configuration& config,
  const Parameters& params,
  const posarrs_t& posarrs,
  const phasearrs_t& phasearrs 
)
{
  int xo=config.xorder();
  size_t nr=posarrs.size();
  size_t nk=phasearrs.size();
  
  meassets_t meassets;
  meassets.resize(phasearrs.size()); 

  measarr_t measarr;
  measarr.resize(nk,0.0);
 
  for(int p=0; p<xo; p++)
  {
    for(size_t k=0; k<nk; k++)
    for(size_t r=0; r<nr; r++)
    {
      size_t index = k*nr + r;
      meas_t dzmag = measutls_realdzmag(config, posarrs[r], phasearrs[k]);
      
      meassets[index].avemag += dzmag;
      meassets[index].varmag += dzmag*dzmag;
      meassets[index].suscep += dzmag*dzmag;
    }
    config.propagate();  
  }

  for(size_t k=0; k<nk; k++)
  for(size_t r=0; r<nr; r++)
  {
    size_t index = k*nr + r;
    meassets[index].suscep = (params.bt / (xo*(xo+1)))*
                         (
                            meassets[index].suscep
                          + (meassets[index].avemag*meassets[index].avemag)
                         ); 
    meassets[index].avemag = fabs(meassets[index].avemag) / (double)xo;
    meassets[index].varmag = meassets[index].varmag / (double)xo;
    meassets[index].suscep = meassets[index].suscep;
  }

  return meassets;
}


}
