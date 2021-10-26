#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <array>

namespace spinhalfsse
{


typedef std::array<std::string,2> paramline_t; 
typedef std::vector<paramline_t> paramfile_t;


struct Parameters
{
  double bt; // inverse temperature
  double ds; // In lattice easy axis anisotropy
  double Js; // In lattice Heisenberg exchange
  double Jp; // On site exchange
  double H;  // Onsite magnetic field
  double C;  // Global diagonal shift
  double R;  // Sublattice rotation
  double G;  // Add to/Subtract from C
  
  size_t nequil; // number of equilibration steps
  size_t nsimul; // number of simulation steps
  size_t mspace; // spacing between measurements

  size_t nx;  // number of sites in x
  size_t ny;  // number of sites in y
  bool   bcx; // boundary conditions in x
  bool   bcy; // boundary conditions in y

  std::string outdir; // name output directory
  std::string outfile; // name of output file

  Parameters() {} 

  Parameters
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
  );


  /*
   *  Accepts paramfile_t produced by extractparams function.
   */
  Parameters(const paramfile_t&);

  void write();
};

typedef std::vector<Parameters> paramlst_t;

}

#endif // PARAMETERS_H
