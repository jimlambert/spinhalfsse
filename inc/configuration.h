#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#define DEFAULT_XO 20

#include <vector>
#include <array>

#include "lattice.h"
#include "randomhelpers.h"
#include "spinHalfVertices.h"

namespace spinhalfsse
{   

typedef std::array<int, NVERTS> vertexledger_t; 
typedef std::vector<int> state_t;
typedef std::vector<int> vtlst_t; 
typedef std::vector<LatticeUtls::enumbond_t> bdlst_t; 


typedef 
struct IndexSpin
{
  int index;
  int spin;
} indspn_t;


/*
 *  Class containing current Monte Carlo configuration. The oplst contains a
 *  list of operators according to propagation index. 
 */
class Configuration
{
  private:

    int _pindex=0; // propagation index
    int _xorder=DEFAULT_XO; // expansion order
    int _nverts=0; // number of non-trivial operators
    int _nsites; // number of sites

    state_t _state; // spin state at current propagation index
    vtlst_t _vtlst; // list of current vertices
    bdlst_t _bdlst; // bonds for current vertices
      
    std::random_device _rd;
    std::mt19937 _mteng{_rd()};
    std::uniform_int_distribution<> _rspin{-1, 1};
    
    // Random site
    std::uniform_int_distribution<>* _rsite;

    // Configuations starts with sign 1 by default
    int _sgn=1;

    vertexledger_t _vertledger;

  public:

    Configuration(const int&);
    
    ~Configuration() { delete _rsite; }

    void propagate();
    void resize(const int&);

    int pindex() const { return _pindex; }
    int xorder() const { return _xorder; }
    int nverts() const { return _nverts; }
    int nsites() const { return _nsites; }
    int getsgn() const { return _sgn; }

    int getspin(const int& i)   const { return _state[i]; }
    int getvt(const int& i)     const { return _vtlst[i]; }
    int getvtnum(const int& v)  const { return _vertledger[v]; }
    int getbdtype(const int& i) const { return _bdlst[i].type; }
    
    LatticeUtls::enumbond_t getbd(const int& i) const { return _bdlst[i]; }

    void setspin(indspn_t indspn) { _state[indspn.index] = indspn.spin; }
    void setvt(const int& p, const int& vt) { _vtlst[p] = vt; }
    void setxo(const int& p) { _xorder = p; }
    void setno(const int& n) { _nverts = n; }
    void addvt(const int& v) { _vertledger[v] += 1; }
    void decvt(const int& v) { _vertledger[v] -= 1; }
    
    void setbd(const int& p, const LatticeUtls::enumbond_t& enmbd) { _bdlst[p] = enmbd; }   
   
    void randspin(const int& n) { _state[n] = _rspin(_mteng); }

    void setsgn(const int& f) { _sgn = f * _sgn; }

    indspn_t getrandomspin() 
    {
      int index = (*_rsite)(_mteng);
      return indspn_t{index, _state[index]};
    }
    
    void printstate();
    void printverts();
    void printbonds();
};


/*
 * Used for N=2, initialize configuration to both spins equal to zero.
 */
void initzeroconfig(Configuration&);


}


#endif // CONFIGURATION_H
