#include <random>
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>

#include "iofunctions.h"
#include "spinHalfVertices.h"
#include "configuration.h"

using namespace std;

namespace spinhalfsse
{


Configuration::Configuration(const int& nsites, const string initType) : _nsites(nsites)
{
  // start at a random site
  _rsite = new std::uniform_int_distribution<>{0, _nsites-1}; 
  if(initType=="RD")
    for(int i=0; i<nsites; i++)
    {
      if(_rspin(_mteng)==1) _state.push_back(1);
      else _state.push_back(-1);
    } 
  else if(initType=="UP")
    for(int i=0; i<nsites; i++) _state.push_back(1);
  else if(initType=="DN")
    for(int i=0; i<nsites; i++) _state.push_back(-1);
  else if(initType=="NL")
    for(int i=0; i<nsites; i++) _state.push_back(1*pow(-1,i));
  
  // Initialize vertex and bond lists
  _vtlst.resize(_xorder, -1);
  _bdlst.resize(_xorder, LatticeUtls::enumbond_t{-1,-1, LatticeUtls::bond_t{-1,-1}});

  // Initialize ledger
  for(int i=0; i<NVERTS; i++) _vertledger[i] = 0;
}


void Configuration::propagate()
{
  int vrtid = _vtlst[_pindex];
 
  if((vrtid == -1) || TYPETBL[vrtid]==0) _pindex = (_pindex+1) % _xorder;
  else 
  {
    int spin1 = VERTEXTBL[vrtid][2];
    int spin2 = VERTEXTBL[vrtid][3];
   
    int site1 = _bdlst[_pindex].bond[0];
    int site2 = _bdlst[_pindex].bond[1];

    _state[site1] = spin1;
    _state[site2] = spin2;

    _pindex = (_pindex+1) % _xorder; 
  }
}


void Configuration::resize(const int& newxo)
{
  _vtlst.resize(newxo, -1);
  _bdlst.resize(newxo, LatticeUtls::enumbond_t{-1,-1, LatticeUtls::bond_t{}});
  _xorder = newxo;
}


void Configuration::printstate() 
{
  cout << "---" << endl;
  
  for(int i=0; i<_nsites; i++) 
    cout << setfill(' ') << setw(5) << left << _state[i];
  cout << endl;
  
  cout << "---" << endl;
}


void Configuration::printverts()
{
  for(int p=0; p<_xorder; p++)
  {
    string index = "[" + to_string(p) +  "]";
    cout << setfill(' ') 
         << setw(8) << left << index 
         << setw(8) << left << _vtlst[p] 
         << endl;
  }
}


void Configuration::printbonds()
{
  for(int p=0; p<_xorder; p++)
  {
    string index = "[" + to_string(p) + "]";
    cout << setw(8) << left << index 
         << setw(4) << left << _bdlst[p].index << setw(2) << "|"
         << setw(4) << left << _bdlst[p].type  << setw(2) << "|"
         << setw(8) << left << _bdlst[p].bond[0]
         << setw(8) << left << _bdlst[p].bond[1]
         << endl;  
  }
}


void initzeroconfig(Configuration& config)
{
  const int pindex = 10;

  if(config.nsites()!=2) 
  {
    cout << "initzeroconfig(Configuraton&) requires a configuration with two sites" 
         << endl; 
  }

  config.setspin(indspn_t{0,0});
  config.setspin(indspn_t{1,0}); 

  config.setvt(pindex, 4);
  config.setbd(pindex, LatticeUtls::enumbond_t{0,0,LatticeUtls::bond_t{0,1}});
  config.setno(1);
} 


} // namespace spin1sse
