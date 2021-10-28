#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <string>
#include <fstream>

#include "initfunctions.h"
#include "spinHalfVertices.h"
#include "iofunctions.h"
#include "updatefunctions.h"
#include "graphinits.h"
#include "Accumulator.h"

using namespace std;
using namespace spinhalfsse;
using namespace LatticeUtls;


int main(int argc, char *argv[])
{
  Parameters params, inPlaneParams, onSiteParams;
  paramlst_t paramlst;

  if(argc==1) // initializeparameters explicitly
  {
    double bt = 30;  // inverse temperature
    double ds = 1.0;   // easy axis anisotropy
    double Js = 1.0;   // in plane heisenberg exchange
    double Jp = 0.0;   // on site heisenberg exchange
    double H  = 0.0; // applied magnetic field along z

    double C  =  0.25; // offset term
    double R  = -1;  // sublattice rotation
    double G  = -1;  // add or subtract Hamiltonian
   
    size_t nequil = 10000;
    size_t nsimul = 100000;
    size_t mspace = 1;
   
    size_t nx = 1;
    size_t ny = 10;
    bool bcx = false;
    bool bcy = true;

    string outname = "test.out";

    // initialize paramters
      params = Parameters{bt,ds,Js,Jp,H,C,R,G,nequil,nsimul,mspace,nx,ny,bcx,bcy, "./", outname};
  }
  else 
  {
    paramfile_t paramfile = extractparams(argv[1]);
    params = Parameters{paramfile};
  }
  params.write();

  inPlaneParams = Parameters
  {
   params.bt,
   params.ds,
   params.Js,
   0.0,
   0.0,
   params.C,
   params.R,
   params.G,
   params.nequil,
   params.nsimul,
   params.mspace,
   params.nx,
   params.ny,
   params.bcx,
   params.bcy,
   params.outdir,
   params.outfile
  };

  onSiteParams = Parameters
  {
   params.bt,
   0.0,
   0.0,
   params.Jp,
   params.H,
   params.C,
   params.R,
   params.G,
   params.nequil,
   params.nsimul,
   params.mspace,
   params.nx,
   params.ny,
   params.bcx,
   params.bcy,
   params.outdir,
   params.outfile
  };

  paramlst.push_back(inPlaneParams); // x bond parameters
  paramlst.push_back(inPlaneParams); // y bond parameters
  paramlst.push_back(onSiteParams);


// LATTICE AND CONFIGURATION
// =========================

  int nx = params.nx;
  int ny = params.ny;
  int nsites = 2*nx*ny;
  bool bcx = params.bcx;
  bool bcy = params.bcy;

  graph_t sqKondoNecklace = LatticeUtls::initSqKondoNecklace(nx,ny,bcx,bcy);
  Lattice sqKondoLattice(sqKondoNecklace);

  Configuration config(nsites, "RD");

  RandomHelpers randhelp;

// INITIALIZE PROBABILITY TABLES
// =============================

  // weights
  wgtlsts_t wgtlsts;
  
  wgtlsts.push_back(initwgts(paramlst[0]));
  wgtlsts.push_back(initwgts(paramlst[1]));
  wgtlsts.push_back(initwgts(paramlst[2]));

  // Transfer lists
  trnlst_t trnlstup = inittrnlst(2);
  trnlst_t trnlstdn = inittrnlst(-2);
  trnlstset_t trnLstSet{trnlstup,trnlstdn};

  // Probability distributions for in plane x bonds
  prblst_t prbLstUpInPlaneX = initprblst(trnlstup,wgtlsts[0]);
  prblst_t prbLstDnInPlaneX = initprblst(trnlstdn,wgtlsts[0]);
  prblst_t cpdLstUpInPlaneX = initcpdlst(prbLstUpInPlaneX);
  prblst_t cpdLstDnInPlaneX = initcpdlst(prbLstDnInPlaneX); 
  prblstset_t cpdLstSetInPlaneX{cpdLstUpInPlaneX,cpdLstDnInPlaneX};
  
  // Probability distributions for in plane y bonds
  prblst_t prbLstUpInPlaneY = initprblst(trnlstup,wgtlsts[1]);
  prblst_t prbLstDnInPlaneY = initprblst(trnlstdn,wgtlsts[1]);
  prblst_t cpdLstUpInPlaneY = initcpdlst(prbLstUpInPlaneY);
  prblst_t cpdLstDnInPlaneY = initcpdlst(prbLstDnInPlaneY); 
  prblstset_t cpdLstSetInPlaneY{cpdLstUpInPlaneY,cpdLstDnInPlaneY};
  
  // Probability distributions for on site bonds 
  prblst_t prbLstUpOnSite = initprblst(trnlstup,wgtlsts[2]);
  prblst_t prbLstDnOnSite = initprblst(trnlstdn,wgtlsts[2]);
  prblst_t cpdLstUpOnSite = initcpdlst(prbLstUpOnSite);
  prblst_t cpdLstDnOnSite = initcpdlst(prbLstDnOnSite); 
  prblstset_t cpdLstSetOnSite{cpdLstUpOnSite,cpdLstDnOnSite};

  prblstsets_t cpdLstSets;

  cpdLstSets.push_back(cpdLstSetInPlaneX);
  cpdLstSets.push_back(cpdLstSetInPlaneY);
  cpdLstSets.push_back(cpdLstSetOnSite);

// SIMULATION
// ==========

  for(size_t i=0; i<params.nequil; i++)
  {
    progress("Equilibration", params.nequil-1, i);

    diagupdt
    (
      config,
      sqKondoLattice,
      params,
      wgtlsts,
      randhelp 
    );

    oneloopupdt
    (
      config, 
      trnLstSet,
      cpdLstSets,
      randhelp
    );

    xordupdt(config);
  }  

  int binsize = 100;

  Accumulator EnergyAccumulator(binsize);
  Accumulator OpNumAccumulator(binsize);
  Accumulator OpNumSqAccumulator(binsize);
  
  for(size_t i=0; i<params.nsimul; i++)
  {
    progress("Simulation", params.nsimul-1, i); 
    diagupdt
    (
      config, 
      sqKondoLattice, 
      params, 
      wgtlsts, 
      randhelp
    );

    oneloopupdt
    (
      config,
      trnLstSet,
      cpdLstSets,
      randhelp   
    );
    
    int nvt = config.nverts();

    OpNumAccumulator.push_back(nvt);
  }

  Measurement opnum = OpNumAccumulator.getMeasurement();
  double opnumerr = opnum.err / opnum.ave;
  double energy = params.G*((1.0/params.bt)*(opnum.ave) 
                       -sqKondoLattice.nbonds()*params.C);
  
  double energyerr = abs(energy)*opnumerr;
  cout << "Energy: " << '\t' << setprecision(10)
       << energy << '\t' << "+/-" << energyerr << endl;
  
  return 0;
}
