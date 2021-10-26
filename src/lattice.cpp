#include "lattice.h"

namespace LatticeUtls
{


bond_t Lattice::getbond(const int& bond)
{
  int total = 0;

  for(int bdlst=0; bdlst<_nbdlsts; bdlst++)
  {
    int bdlstsize = _graph[bdlst].size();
    if(bond < (total+bdlstsize)) 
      return _graph[bdlst][bond-total];
    total += bdlstsize;
  }

  return bond_t{}; 
}


bondlst_t Lattice::getbondlst(const int& i) { return _graph[i]; }

enumbond_t Lattice::getrandombond()
{
  int randbond = (*_rbond)(_mteng);
  
  int total = 0.0;

  for(int bdlst=0; bdlst<_nbdlsts; bdlst++)
  {
    int bdlstsize = _graph[bdlst].size();
    
    if(randbond < (total+bdlstsize))
      return enumbond_t{randbond-total, bdlst, _graph[bdlst][randbond-total]};
    
    total += bdlstsize;
  } 
  
  return enumbond_t{};
}


} // namespace spin1sse
