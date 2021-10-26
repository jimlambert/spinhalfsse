#ifndef SPINHALF_VERTICES_H
#define SPINHALF_VERTICES_H

#include <array>
#include <vector>

#define DIAGVERT_SUP 4

namespace spinhalfsse
{

const int NVERTS = 6;

typedef std::array<int,4> vert_t;

const int VERTEXTBL[NVERTS][4]
{
  { 1,  1,  1,  1}, // 0
  {-1, -1, -1, -1}, // 1
  {-1,  1, -1,  1}, // 2
  { 1, -1,  1, -1}, // 3
  {-1,  1,  1, -1}, // 4
  { 1, -1, -1,  1}, // 5
};


const int TYPETBL[NVERTS]
{
  0,
  0,
  0,
  0,
  1,
  1
};


/*
 *  Accepts four spin values.
 *
 *  If valid vertex returns index in vertex list. Else returns -1. 
 */
inline int findvrt(const int& a, const int& b, const int& c, const int& d)
{
  for(int v=0; v<NVERTS; v++)  
  {
    if
    (
      VERTEXTBL[v][0]==a && 
      VERTEXTBL[v][1]==b && 
      VERTEXTBL[v][2]==c && 
      VERTEXTBL[v][3]==d
    ) 
    return v; 
  }
  return -1;
}


}


#endif // SPIN_HALF_VERTICES_H
