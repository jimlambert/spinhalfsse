#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <array>
#include <random>

namespace LatticeUtls
{

// == TYPEDEFS == //
typedef std::array<int, 2> bond_t;
typedef std::vector<bond_t> bondlst_t;
typedef std::vector<bondlst_t> graph_t;

typedef 
struct Enumbond
{
  int index; // index within bond list
  int type;  // type of bond 
  bond_t bond; 
} enumbond_t;


// == LATTICE CLASS == //
/*
 * Class acts as wrapper for a graph with functionality to randomly access sites
 * and bonds.
 *
 */
class Lattice
{

  private:
    
    graph_t _graph;
    
    int _nbonds=0;
    int _nbdlsts;

    std::random_device _rd;
    std::mt19937 _mteng{_rd()};
    std::uniform_int_distribution<>* _rbond;

  public:

    Lattice(const graph_t& graph) : 
      _graph(graph),
      _nbdlsts(graph.size())
    { 
      for(int i=0; i<_nbdlsts; i++) _nbonds += graph[i].size(); 
      _rbond = new std::uniform_int_distribution<>{0, _nbonds-1};
    }

    ~Lattice() {delete _rbond;}

    int nbonds() const { return _nbonds; }
    int nbdlsts() const { return _nbdlsts; }
    graph_t graph() const { return _graph; }
    
    bond_t getbond(const int&);
    bondlst_t getbondlst(const int&);

    enumbond_t getrandombond();
};


}

#endif // LATTICE_H
