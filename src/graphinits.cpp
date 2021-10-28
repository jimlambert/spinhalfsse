#include <iostream>

#include "graphinits.h"

using namespace std;

namespace LatticeUtls
{

/*
 *  Accepts number of sites and a bool which is true for periodic boundary
 *  conditions and false for open.
 *
 *  Returns graph corresponding to nearest neighbour bonds 
 */
graph_t initchain(const int& nsites, const bool& bc)
{
  bondlst_t bondlst;

  if(bc) // periodic boundary conditions
    for(int i=0; i<nsites; i++) 
      bondlst.push_back(bond_t{i, (i+1) % nsites});

  else // open boundary conditions
    for(int i=0; i<nsites-1; i++) 
      bondlst.push_back(bond_t{i, (i+1) % nsites});
  
  graph_t graph;
  graph.push_back(bondlst);

  return graph;
}


graph_t initisosqlatt
(
  const int& lx, 
  const int& ly, 
  const bool bcx, 
  const bool bcy
)
{
  graph_t graph;
  bondlst_t bondlst;

  for(int y=0; y<ly; y++)
  for(int x=0; x<lx-(bcx^1); x++)
  {
    int xmin = y*lx;
    int xmax = (y+1)*lx;  
    bondlst.push_back(bond_t{x+xmin, (x+xmin+1)%xmax + (xmin*(int)((x+1)/lx))});
  }
 
  for(int y=0; y<ly-(bcy^1); y++)
  for(int x=0; x<lx; x++)
  {
    int xmin = y*lx;
    bondlst.push_back(bond_t{x+xmin, (x+xmin+lx)%(lx*ly)}); 
  }

  graph.push_back(bondlst);
  
  return graph;
}


graph_t initsqlatt
(
  const int& lx, 
  const int& ly, 
  const bool bcx, 
  const bool bcy
)
{
  graph_t graph;
  bondlst_t xbondlst;
  bondlst_t ybondlst;
  
  for(int y=0; y<ly; y++)
  for(int x=0; x<lx-(bcx^1); x++)
  {
    int xmin = y*lx;
    int xmax = (y+1)*lx;  
    xbondlst.push_back(bond_t{x+xmin, (x+xmin+1)%xmax + (xmin*(int)((x+1)/lx))});
  }
 
  for(int y=0; y<ly-(bcy^1); y++)
  for(int x=0; x<lx; x++)
  {
    int xmin = y*lx;
    ybondlst.push_back(bond_t{x+xmin, (x+xmin+lx)%(lx*ly)}); 
  }

  graph.push_back(xbondlst);
  graph.push_back(ybondlst);

  return graph;
}


graph_t inithoneycomb
(
  const int& lx, 
  const int& ly, 
  const bool bcx,
  const bool bcy
)
{
  int chlen = 2*lx; // chain length
  int nch = ly;     // number of chains
  int nsites = 2*lx*ly;

  graph_t graph;
  bondlst_t bondlst;

  for(int chain=0; chain<nch; chain++)
  {
    int offset = chain*chlen;

    // initialize x and y bonds
    for(int r=0; r<chlen-(bcx^1); r++)
    {
      int rstrt = offset;
      int rend = offset + chlen - 1;

      int rpos2, rpos1 = r+offset;
      if(rpos1==rend) rpos2=rstrt;
      else rpos2=rpos1 + 1;

      if(rpos1%2 == 0) bondlst.push_back(bond_t{rpos1,rpos2});
      else             bondlst.push_back(bond_t{rpos2,rpos1});
    }
  
  }
    
  // initialize z bonds
  for(int chain=0; chain<nch-(bcy^1); chain++)
  {
    int offset = chain*chlen;

    for(int r=1; r<chlen; r+=2)
    {
      int rpos1 = r + offset;
      int rpos2 = (r + offset + chlen - 1) % nsites;

      bondlst.push_back(bond_t{rpos2, rpos1});
    }
  }

  graph.push_back(bondlst);

  return graph;
}


graph_t initIsoSqKondoNecklace
(
  const int& lx,
  const int& ly,
  const bool bcx,
  const bool bcy 
)
{
  graph_t graph;
  bondlst_t inPlaneBondLst;
  bondlst_t danglingBondLst;

  int N = lx*ly;

  for(int y=0; y<ly; y++)
  for(int x=0; x<lx-(bcx^1); x++)
  {
    int xmin = y*lx;
    int xmax = (y+1)*lx;  
    inPlaneBondLst.push_back(bond_t{x+xmin, (x+xmin+1)%xmax + (xmin*(int)((x+1)/lx))});
  }
 
  for(int y=0; y<ly-(bcy^1); y++)
  for(int x=0; x<lx; x++)
  {
    int xmin = y*lx;
    inPlaneBondLst.push_back(bond_t{x+xmin, (x+xmin+lx)%(lx*ly)}); 
  }

  for(int y=0; y<ly; y++)
  for(int x=0; x<lx; x++)
  {
    int xmin = y*lx;
    danglingBondLst.push_back(bond_t{x+xmin, x+xmin+N});
  }

  graph.push_back(inPlaneBondLst);
  graph.push_back(danglingBondLst);
 

  return graph;
}


graph_t initSqKondoNecklace
(
  const int& lx, 
  const int& ly, 
  const bool bcx, 
  const bool bcy
)
{
  int N = lx*ly;
  graph_t graph;
  bondlst_t xbondlst;
  bondlst_t ybondlst;
  bondlst_t danglingBondLst;
  
  for(int y=0; y<ly; y++)
  for(int x=0; x<lx-(bcx^1); x++)
  {
    int xmin = y*lx;
    int xmax = (y+1)*lx;  
    xbondlst.push_back(bond_t{x+xmin, (x+xmin+1)%xmax + (xmin*(int)((x+1)/lx))});
  }
 
  for(int y=0; y<ly-(bcy^1); y++)
  for(int x=0; x<lx; x++)
  {
    int xmin = y*lx;
    ybondlst.push_back(bond_t{x+xmin, (x+xmin+lx)%(lx*ly)}); 
  }

  for(int y=0; y<ly; y++)
  for(int x=0; x<lx; x++)
  {
    int xmin = y*lx;
    danglingBondLst.push_back(bond_t{x+xmin, x+xmin+N});
  }

  graph.push_back(xbondlst);
  graph.push_back(ybondlst);
  graph.push_back(danglingBondLst); 

  return graph;
}



}
