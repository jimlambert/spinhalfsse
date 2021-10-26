#include <iostream>
#include <iomanip>

#include "iofunctions.h"
#include "parameters.h"
#include "utls.h"

#define COLW 10

using namespace std;

namespace spinhalfsse
{

Parameters::Parameters
(
  double bt, // inverse temperature
  double ds, // In lattice easy axis anisotropy
  double Js, // In lattice Heisenberg exchange
  double Jp, // On site exchange
  double H,  // Onsite magnetic field
  double C,  // Global diagonal shift
  double R,  // Sublattice rotation
  double G,  // Add to/Subtract from C
  size_t nequil,
  size_t nsimul,
  size_t mspace,
  size_t nx,
  size_t ny,
  bool   bcx,
  bool   bcy,
  std::string outdir,
  std::string outfile
) 
: bt(bt),
  ds(ds),
  Js(Js),
  Jp(Jp),
  H(H),
  C(C),
  R(R),
  G(G),
  nequil(nequil),
  nsimul(nsimul),
  mspace(mspace),
  nx(nx),
  ny(ny),
  bcx(bcx),
  bcy(bcy),
  outdir(outdir),
  outfile(outfile)
{}


Parameters::Parameters(const paramfile_t& paramfile)
{
  for(auto v : paramfile)
  {
    if(v[0]=="bt")           bt = convert_to<double>(v[1]);
    else if(v[0]=="T")       bt = 1.0/convert_to<double>(v[1]);
    else if(v[0]=="ds")      ds = convert_to<double>(v[1]);
    else if(v[0]=="Js")      Js = convert_to<double>(v[1]);
    else if(v[0]=="Jp")      Jp = convert_to<double>(v[1]);
    else if(v[0]=="H")       H = convert_to<double>(v[1]);
    else if(v[0]=="C")       C = convert_to<double>(v[1]);
    else if(v[0]=="R")       R = convert_to<double>(v[1]);
    else if(v[0]=="G")       G = convert_to<double>(v[1]);
    else if(v[0]=="nequil")  nequil = convert_to<int>(v[1]);
    else if(v[0]=="nsimul")  nsimul = convert_to<int>(v[1]);
    else if(v[0]=="mspace")  mspace = convert_to<int>(v[1]);
    else if(v[0]=="nx")      nx = convert_to<int>(v[1]);
    else if(v[0]=="ny")      ny = convert_to<int>(v[1]);
    else if(v[0]=="bcx")     bcx = convert_to<bool>(v[1]);
    else if(v[0]=="bcy")     bcy = convert_to<bool>(v[1]);
    else if(v[0]=="outfile") outfile = v[1];
    else 
    {
      cout << "UNRECOGNIZED PARAMETER NAME:" << '\t' << v[0] << endl;
      cout << "ASSOCIATED VALUE:" << '\t' << v[1] << endl;
      exit(1);
    }
  }
}


void Parameters::write()
{
  cout << left << setw(COLW) << "Simulation Parameters" << endl;
  cout << left << setw(COLW*2+1) << setfill('-') << '-' << endl;
  cout << setfill(' ') << left << setw(COLW) << "bt:" 
       << left << setw(COLW) << bt 
       << endl;
  cout << left << setw(COLW) << "ds:" 
       << left << setw(COLW) << ds
       << endl;
  cout << left << setw(COLW) << "Js:" 
       << left << setw(COLW) << Js
       << endl;
  cout << left << setw(COLW) << "Jp:" 
       << left << setw(COLW) << Jp
       << endl;
  cout << left << setw(COLW) << "H:" 
       << left << setw(COLW) << H 
       << endl;
  cout << left << setw(COLW) << "C:" 
       << left << setw(COLW) << C 
       << endl;
  cout << left << setw(COLW) << "R:" 
       << left << setw(COLW) << R 
       << endl;
  cout << left << setw(COLW) << "G:" 
       << left << setw(COLW) << G 
       << endl;
  cout << left << setw(COLW) << "nequil:" 
       << left << setw(COLW) << nequil 
       << endl;
  cout << left << setw(COLW) << "nsimul:" 
       << left << setw(COLW) << nsimul 
       << endl;
  cout << left << setw(COLW) << "mspace:" 
       << left << setw(COLW) << mspace 
       << endl;
  cout << left << setw(COLW) << "nx:" 
       << left << setw(COLW) << nx 
       << endl;
  cout << left << setw(COLW) << "ny:" 
       << left << setw(COLW) << ny 
       << endl;
  cout << left << setw(COLW) << "bcx:" 
       << left << setw(COLW) << bcx
       << endl;
  cout << left << setw(COLW) << "bcy:" 
       << left << setw(COLW) << bcy
       << endl;
  cout << left << setw(COLW) << "outfile:" 
       << left << setw(COLW) << outfile
       << endl;
  cout << left << setw(COLW*2+1) << setfill('-') << '-' << endl;
}


}

