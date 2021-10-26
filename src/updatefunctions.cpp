#include <iostream>
#include <iomanip>
#include <assert.h>

#include "spinHalfVertices.h"
#include "updatefunctions.h"
#include "iofunctions.h"

using namespace std;

namespace spinhalfsse
{

void spinupdt
(
  Configuration& config, 
  const freesplst_t& freesplst
)
{
  for(size_t i=0; i<freesplst.size(); i++) 
    if(freesplst[i] == -1) config.randspin(i);
}


void xordupdt(Configuration& config)
{
  int newxo = 1.33 * config.nverts();
  if(newxo > config.xorder()) config.resize(newxo);
}


void diagupdt
(
  Configuration& config, 
  LatticeUtls::Lattice& latt, 
  Parameters& params, 
  wgtlsts_t& wgts,
  RandomHelpers& randhelp
)
{
  int xo = config.xorder();
  int nb = latt.nbonds();

  for(int p=0; p<xo; p++)
  {
    int no = config.nverts();
    double r = randhelp.random();
    if(config.getvt(p) == -1) // attempt insertion
    {  
      LatticeUtls::enumbond_t enumbond = latt.getrandombond();  
      
      int type = enumbond.type; // get bond type
      int spin1 = config.getspin(enumbond.bond[0]);
      int spin2 = config.getspin(enumbond.bond[1]);
      int vrtid = findvrt(spin1, spin2, spin1, spin2);

      if(r*(xo - no) < nb*params.bt*fabs(wgts[type][vrtid]))
      {
        config.setvt(p,vrtid);
        config.addvt(vrtid);
        config.setbd(p,enumbond);
        config.setno(no+1);
      }
    }
    else if(config.getvt(p) < DIAGVERT_SUP) // attempt removal
    {
      int vrtid = config.getvt(p);   
      int type = config.getbdtype(p);
    
      if((r*nb*params.bt*fabs(wgts[type][vrtid])<(xo-no+1)))
      {
        config.setvt(p,-1);
        config.decvt(vrtid);
        config.setbd(p,LatticeUtls::enumbond_t{-1,-1, LatticeUtls::bond_t{}});
        config.setno(no-1);
      }
    }
    config.propagate(); // cycle state
  }  
}


void oneloopupdt
(
  Configuration& config,
  sgnlst_t& sgnlst,
  trnlstset_t& trnlstset, // spin down followed by spin up
  prblstset_t& cpdlstset,
  RandomHelpers& randhelp
)
{
  lnkfrst_t lnkfrst = getlinklst(config);
  
  // attempt to start update at each operator index
  for(int p=0; p<config.xorder(); p++)
  {
    if(config.getvt(p) == -1) continue;
    else
    {
      int sfi, sfc; // initial and current spin flip
      
      int e = randhelp.randleg();
        
      int vt = config.getvt(p);

      // determine direction for initial spin flip      
      if(VERTEXTBL[vt][e] == 1) sfi = 0;
      else if(VERTEXTBL[vt][e] == -1) sfi = 1;
      else
      {
        if(randhelp.random() < 0.5) sfi = 0;
        else sfi = 1;
      }
      
      sfc = sfi; 
      
      int v0 = 4*p + e;
      int vc = v0;
     
      bool endterm = false; // has the head and tail hit a termination vertex? 

      do
      {
        vt = config.getvt((int)vc/4);

        double r = randhelp.random();
        int x=e;

        int oldsgn = sgnlst[vt];

        for(int i=0; i<5; i++)
        {
          double cpd = cpdlstset[sfc][vt][e][i];
          if(cpd < 1e-7) continue; // if probability is zero, skip entry
          if(fabs(cpd - 1) < 1e-7) {x=i; break;} // if probability is 1 accept
          if(r < cpd) { x=i; break; } 
        }
        
        // update sign of configuration
        int nvt = trnlstset[sfc][vt][e][x];
        int newsgn = sgnlst[nvt];
        config.setsgn(newsgn*oldsgn);
        config.decvt(vt);
        config.addvt(nvt);

        // update vertex configuration
        config.setvt((int)(vc/4), nvt); 
        
        if(x==4) // end termination condition
        {
          if(endterm) break; // end loop update
          else // first end terminates
          {
            endterm = true;
            vc = lnkfrst.linklst[v0]; // proceed out of starting vertex
            sfc = sfi; // revert to original spin flip sign
          } 
        }
        else
        { 
          // determine new spin flip direction
          if((e<2 && x<2) || (e>1 && x>1)) sfc = sfc^1; // if loop reverses, flip spin   
          
          // set vertex coordinate and entrance leg for next step
          vc = lnkfrst.linklst[vc - e + x];
        }
        
        // set entrance leg
        e = vc % 4;

      } while((vc != v0) && !endterm);
    }   
  }
  
  // update the spins in the spin configuration
  for(int i=0; i<config.nsites(); i++) 
  {
    if(lnkfrst.vfrst[i] != -1)
    {
      int p = (lnkfrst.vfrst[i]) / 4; // vertex list index
      int l = (lnkfrst.vfrst[i]) % 4; // leg
      config.setspin(indspn_t{i, VERTEXTBL[config.getvt(p)][l]});
    }
    else config.randspin(i);
  }
}


void oneloopupdt
(
  Configuration& config,
  trnlstset_t& trnlstset, // spin down followed by spin up
  prblstsets_t& cpdlstsets,
  RandomHelpers& randhelp
)
{
  lnkfrst_t lnkfrst = getlinklst(config);
  
  // attempt to start update at each operator index
  for(int p=0; p<config.xorder(); p++)
  {
    if(config.getvt(p) == -1) continue;
    else
    {
      int sfi, sfc; // initial and current spin flip
      
      int e = randhelp.randleg();
        
      int vt = config.getvt(p);

      // determine direction for initial spin flip      
      if(VERTEXTBL[vt][e] == 1) sfi = 0;
      else if(VERTEXTBL[vt][e] == -1) sfi = 1;
      else
      {
        if(randhelp.random() < 0.5) sfi = 0;
        else sfi = 1;
      }
      
      sfc = sfi; 
      
      int v0 = 4*p + e;
      int vc = v0;
     
      bool endterm = false; // has the head or tail hit a termination vertex? 

      do
      {
        vt = config.getvt((int)vc/4);
        int tp = config.getbdtype((int)vc/4); // get bond type

        double r = randhelp.random();
        int x=e;

        for(int i=0; i<5; i++)
        {
          double cpd = cpdlstsets[tp][sfc][vt][e][i];
          if(cpd < 1e-7) continue; // if probability is zero, skip entry
          if(fabs(cpd - 1) < 1e-7) {x=i; break;} // if probability is 1 accept
          if(r < cpd) { x=i; break; } 
        }
        
        // update sign of configuration
        int nvt = trnlstset[sfc][vt][e][x];
        config.decvt(vt);
        config.addvt(nvt);

        // update vertex configuration
        config.setvt((int)(vc/4), nvt); 
        
        if(x==4) // end termination condition
        {
          if(endterm) break; // end loop update
          else // first end terminates
          {
            endterm = true;
            vc = lnkfrst.linklst[v0]; // proceed out of starting vertex
            sfc = sfi; // revert to original spin flip sign
          } 
        }
        else
        { 
          // determine new spin flip direction
          if((e<2 && x<2) || (e>1 && x>1)) sfc = sfc^1; // if loop reverses, flip spin   
          
          // set vertex coordinate and entrance leg for next step
          vc = lnkfrst.linklst[vc - e + x];
        }
        
        // set entrance leg
        e = vc % 4;

      } while((vc != v0) && !endterm);
    }   
  }
  
  // update the spins in the spin configuration
  for(int i=0; i<config.nsites(); i++) 
  {
    if(lnkfrst.vfrst[i] != -1)
    {
      int p = (lnkfrst.vfrst[i]) / 4; // vertex list index
      int l = (lnkfrst.vfrst[i]) % 4; // leg
      config.setspin(indspn_t{i, VERTEXTBL[config.getvt(p)][l]});
    }
    else config.randspin(i);
  }
}

}
