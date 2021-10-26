#ifndef ACCUMULATOR_H
#define ACCUMULATOR_H

#include <vector>

struct Measurement
{
  double ave;
  double err;
};


class Accumulator
{
  typedef std::vector<double> binvals_t;

  private:
    
    size_t _binsize; 
    size_t _nmeas=0; // number of measurements in current bin
    double _total=0.0;

    binvals_t _binvals;


  public:

    Accumulator(int bs) : _binsize(bs) {}

    void push_back(double val)
    {
      _total += val;
      _nmeas += 1;

      if(_nmeas == _binsize)
      {
        _binvals.push_back(_total/(double)_nmeas);
        _total=0.0;
        _nmeas=0;
      }    
    }

    int getNumBins() {return _binvals.size();}

    Measurement getMeasurement()
    {
      double tot = 0.0;
      double sqtot = 0.0;
      double nbins = _binvals.size();
      for(auto iter : _binvals)
      {
        tot += iter;
        sqtot += iter*iter;
      }
      tot = tot / nbins;
      sqtot = sqtot / nbins;

      return Measurement{tot, (sqtot - (tot*tot))/nbins};
    }

    double operator[] (const size_t index)
    { return _binvals[index]; }
};

#endif // ACCUMULATOR_H
