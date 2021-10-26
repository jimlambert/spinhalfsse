#ifndef UTLS_H
#define UTLS_H

#include <string>
#include <sstream>

template <typename T>
int sgn(const T& val)
{ return (T(0) < val) - (val < T(0)); }


template <typename T> 
T convert_to(const std::string& str)
{
  std::istringstream ss(str);
  
  T val;
  ss >> val;
  
  return val;
}

#endif // UTLS_H
