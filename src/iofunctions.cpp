#include <iostream>
#include <sstream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

#include "iofunctions.h"

#define COLW 5
#define LLW 5 // spacing for linked list ouptut

using namespace std;

namespace spinhalfsse
{

// TRANSFER MATRIX
void write(const trnmat_t& trnmat)
{
  for(int e=0; e<4; e++)
  {
    for(int x=0; x<5; x++)
      cout << setw(10) << setprecision(3) << left << trnmat[e][x];
    cout << endl;
  }
}


// LIST OF TRANSFER MATRICES
void write(const trnlst_t& trnlst)
{
  for(int v=0; v<NVERTS; v++)
  {
    cout << "vertex: [" << v << "]" << endl;  
    for(int e=0; e<4; e++)
    {
      for(int x=0; x<5; x++)
        cout << setw(COLW) << left << trnlst[v][e][x];     
      cout << endl;
    }
  }
}


// PROBABILITY MATRIX
void write(const prbmat_t& prbmat)
{
  for(int e=0; e<4; e++)
  {
    for(int x=0; x<5; x++)
      cout << setw(10) << setprecision(3) << left << prbmat[e][x]; 
    
    cout << endl;
  }
}


// LIST OF PROBABILITY MATRICES
void write(const prblst_t& prblst)
{
  for(int v=0; v<NVERTS; v++)
  {
    cout << "index:" << '\t' << v << endl;
    write(prblst[v]);  
    cout << "----" << endl;
  }
}


// WEIGHT LIST
void write(const wgtlst_t& wgtlst)
{
  for(int v=0; v<NVERTS; v++)
  {
    string index = "[" + to_string(v) + "]";
    cout << setw(COLW) 
         << left << index 
         << setw(COLW) << left << wgtlst[v]
         << endl; 
  }
}


// SIGN LIST
void write(const sgnlst_t& sgnlst)
{
  for(int v=0; v<NVERTS; v++)
  {
    string index = "[" + to_string(v) + "]";
    cout << setw(COLW) 
         << left << index 
         << setw(COLW) << left << sgnlst[v]
         << endl; 
  }
}


// BOND
void write(const LatticeUtls::bond_t& bond)
{
  cout << setw(5) << left << bond[0]
       << setw(5) << left << bond[1] 
       << endl;
}


// BOND LIST
void write(const LatticeUtls::bondlst_t& bondlst)
{
  for(size_t i=0; i<bondlst.size(); i++)
  {
    string index = "[" + to_string(i) + "]";
    cout << setw(8) << left << index; 
    write(bondlst[i]);
  }
  cout << "---" << endl;
}


// LIST OF BOND LISTS
void write(const LatticeUtls::graph_t& graph)
{
  for(auto git=graph.begin(); git!=graph.end(); git++)
    write((*git));   
}


// LINKED LIST
void write(const linklst_t& linklst)
{
  for(size_t i=0; i<linklst.size(); i++)
  {
    string index = "[" + to_string(i) + "]";
    cout << setw(LLW) << left << index << setw(LLW) << left << linklst[i];
    if(i%4 == 3) cout << endl;
  }
}


void write(const phasearr_t& phasearr)
{
  for(size_t i=0; i<phasearr.size(); i++)
  {
    string index = "[" + to_string(i) + "]";
    cout << setw(LLW) << left << index << setw(LLW) << left 
         << phasearr[i] << endl;
  }
}


// WVERTEX
void writevertex(const int& vtindex)
{
  for(int i=0; i<4; i++) cout << VERTEXTBL[vtindex][i] << '\t';
  cout << endl;
}


void progress(const string& procname, const int& maxstp, const int& curstp)
{
  int lenprnm = procname.length();
  string comper = to_string((int)(100*((double)curstp/(double)maxstp))) + "%";
  cout << setw(lenprnm) << left << procname 
       << setw(2) << setfill(' ') << ' '
       << setw(8) << left << comper << '\r';
  if(curstp==maxstp)
    cout << setw(lenprnm+10) << setfill(' ') << procname << " complete" << endl;
  
}


// READ PARAMETER FILE
paramfile_t extractparams(char* fname)
{
  ifstream infile;
  infile.open(fname);
  
  string line;

  paramfile_t paramfile;
 
  while(getline(infile,line,'\n'))
  {
    stringstream ss(line);
    string ele;

    paramline_t paramline;

    int index = 0;
    while(getline(ss,ele,'='))
    {
      paramline[index] = ele;
      index = index^1;
    }

    paramfile.push_back(paramline);
  }

  paramfile.push_back(paramline_t{"outfile",convertfilename(fname)});

  return paramfile;
}


string convertfilename(char* fname)
{
  string line;
  stringstream ss(fname);

  getline(ss, line, '.');

  return line+".out";
}

}
