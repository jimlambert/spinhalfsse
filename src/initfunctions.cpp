#include <iostream>

#include "spinHalfVertices.h"
#include "initfunctions.h"
#include "parameters.h"
#include "utls.h"

using namespace std;

namespace spinhalfsse
{


/*
 *  Returns a probabiities matrix with all zero entries. 
 */
prbmat_t zeroprbmat()
{
  prbmat_t prbmat;
  for(int x=0; x<5; x++)
  for(int e=0; e<4; e++)
    prbmat[e][x] = 0;
  
  return prbmat;
}


/*
 *  Accepts model parameters.
 *
 *  Returns wgtlst
 */
wgtlst_t initwgts(const Parameters& params)
{
  wgtlst_t wgtlst;
  for(int i=0; i<NVERTS; i++)
  {
    if(TYPETBL[i]==0)
    {
      int s1=VERTEXTBL[i][0];
      int s2=VERTEXTBL[i][1];
      
      wgtlst[i]=params.C
               +params.G*(params.ds*params.Js*0.25*s1*s2
                         +params.Jp*0.25*s1*s2
                         +0.5*params.H*(s1+s2)
                         );
    }
    else if(TYPETBL[i]==1)
    {
      wgtlst[i]=params.G*params.R*(0.5*params.Js
                                  +0.5*params.Jp   
                                  );
    }
    else cout << "UNKNOWN TYPE IN INIT FUNCTION" << endl;
  }
  return wgtlst;
}


sgnlst_t initsgns(const wgtlst_t& wgtlst)
{
  sgnlst_t sgnlst;
  for(int i=0; i<NVERTS; i++)
    sgnlst[i] = sgn<double>(wgtlst[i]);

  return sgnlst;
}


/*
 *  Accepts magnitude of spin flip at loop entrance.
 *
 *  Returns list of transformed vertices
 */
trnlst_t inittrnlst(const int& sflip)
{
  trnlst_t trnlst;
  for(int v=0; v<NVERTS; v++)
  {
    trnmat_t trnmat;

    for(int e=0; e<4; e++)
    for(int x=0; x<4; x++)
    {
      int vert[4] = 
      {
        VERTEXTBL[v][0], 
        VERTEXTBL[v][1], 
        VERTEXTBL[v][2], 
        VERTEXTBL[v][3]
      };

      vert[e]=vert[e] + sflip;
      
      // Find terminal vertex
      trnmat[e][4] = findvrt(vert[0],vert[1],vert[2],vert[3]);

      // Find pass-through vertex
      if((e<2 && x<2) || (e>1 && x>1)) vert[x]=vert[x] - sflip;
      else vert[x]=vert[x]+sflip;

      trnmat[e][x] = findvrt(vert[0], vert[1], vert[2], vert[3]);  
    }
    
    trnlst[v]=trnmat;
  }

  return trnlst;
}


/*
 *  Accepts a transfer list and a list of weights.
 *
 *  Returns probability for each exit leg given an entrace leg. 
 */
prblst_t initprblst(const trnlst_t& trnlst, const wgtlst_t& wgtlst)
{
  prblst_t prblst;
  for(int v=0; v<NVERTS; v++)
  {
    prbmat_t prbmat=zeroprbmat();

    for(int e=0; e<4; e++)
    {
      double total=0.0;

      for(int x=0; x<5; x++) 
        if(trnlst[v][e][x] > -1) total += fabs(wgtlst[trnlst[v][e][x]]);
   
      if(fabs(wgtlst[v]) < 1e-7) prbmat[e][e] = 1; // if vertex has no weight bounceback
      else
      {
        for(int x=0; x<5; x++) 
          if(trnlst[v][e][x] > -1) prbmat[e][x]=fabs(wgtlst[trnlst[v][e][x]]/total);
      }
    }

    prblst[v]=prbmat;
  }

  return prblst;
}


/*
 *  Accepts a list of probability matrices.
 *
 *  Returns a list of cumulative probabilities.
 */
prblst_t initcpdlst(const prblst_t& prblst)
{
  prblst_t cpdlst;
  for(int v=0; v<NVERTS; v++)
  {
    prbmat_t prbmat = prblst[v];
    prbmat_t cpdmat = zeroprbmat(); 
    
    for(int e=0; e<4; e++)
    {
      double total = 0.0;
      for(int x=0; x<5; x++)
      {
        total += prbmat[e][x];
        if(prbmat[e][x]<1e-15) continue;
        else cpdmat[e][x] = total;  
      }
    } 
    cpdlst[v]=cpdmat;
  }

  return cpdlst;
}


/*
 *  Accepts a configuration.
 *
 *  Returns the linked list corresponding to the operators in the configuration.
 */
lnkfrst_t getlinklst(Configuration& config)
{
  int xo = config.xorder();
  int ns = config.nsites();
  linklst_t linklst; 
  std::vector<int> vfrst;
  std::vector<int> vlast;
  linklst.resize(4*xo, 0);
  vfrst.resize(ns, -1);
  vlast.resize(ns, -1);

  int v0, v1, v2, s1, s2;

  for(int p=0; p<xo; p++)
  {
    v0 = 4 * p;
    if(config.getvt(p) == -1)
      for(int i=0; i<4; i++) linklst[v0 + i] = -1;
    else
    { 
      s1 = config.getbd(p).bond[0];
      s2 = config.getbd(p).bond[1];
      v1 = vlast[s1];
      v2 = vlast[s2];
      
      if(v1 != -1)
      {
        linklst[v1] = v0; 
        linklst[v0] = v1;
      }
      else vfrst[s1] = v0;
      
      if(v2 != -1)
      {
        linklst[v2]   = v0 + 1;
        linklst[v0+1] = v2;
      }
      else vfrst[s2] = v0 + 1;
      
      vlast[s1] = v0 + 2;
      vlast[s2] = v0 + 3;
    }
  }

  // close links across imaginary time boundary
  for(int i=0; i<ns; i++){
    if(vfrst[i] != -1){
      linklst[vlast[i]] = vfrst[i];
      linklst[vfrst[i]] = vlast[i];
    }
  }
    
  return lnkfrst_t{linklst,vfrst};
}
}
