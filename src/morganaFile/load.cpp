/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "load.h"

load::
load()
{
  strProjects   = "projects/";
  strGeometries = "geometries/";
  strResumePath = "resume/";
}

void
load::
controllo(ifstream &file, string &file_name)
{
  if(!file.good())
  {
    cout << "ERRORE: File " << file_name << " inesistente"<< endl;
    file.close();
  }
}

UInt
load::
leggiriga(ifstream &file)
{
  //Legge una riga
  string strRead;
  getline(file, strRead, '\n');
  
  //Split della riga
  istringstream ISS(strRead);
  VS = vector<string>(istream_iterator<string>(ISS), istream_iterator<string>());
  
  //Rimozione dei commenti
  for(unsigned int i = 0; i<VS.size(); i++)
  {
    if (VS[i] == "//")
    {
      VS.erase (VS.begin()+i,VS.end());
      break;
    }
  }
  
  return(VS.size());
}

void
load::
assegnavalore(const UInt &i, UInt &val)
{
  stringstream SS;
  
  SS << VS[i]; 
  SS >> val;
}

void
load::
assegnavalore(const UInt &i, Real &val)
{
  stringstream SS;
  
  SS << VS[i];
  SS >> val;
}

void
load::
assegnavalore(const UInt &i, Real *val)
{
  stringstream SS;
  
  SS << VS[i];
  SS >> *val;
}

void
load::
assegnavalore(const UInt &i, string &val)
{
  val = VS[i];
}

void
load::
assegnavalore(const UInt &i, char &val)
{
  stringstream SS;
  
  SS << VS[i];
  SS >> val;
}

void
load::
stoupper(string &s)
{
  for (string::iterator IT = s.begin() ; IT != s.end() ; IT++)
    *IT = toupper((unsigned char)*IT);
}

void
load::
stolower(string &s)
{
  for (string::iterator IT = s.begin() ; IT != s.end() ; IT++)
    *IT = tolower((unsigned char)*IT);
}

UInt
load::
segmentationUpper(const UInt & pid, const UInt & npid, const UInt & n) const
{
  return( UInt(( (Real(pid + 1)) / Real(npid) ) * Real(n)) );
}

UInt
load::
segmentationLower(const UInt & pid, const UInt & npid, const UInt & n) const
{
  return( UInt(( Real(pid) / Real(npid) ) * Real(n)) + 1 );
}

string
load::
readString(ifstream & file)
{
  string outString, item;
  UInt size = leggiriga(file);
  
  for(UInt i=0; i<size; ++i)
  {
    assegnavalore(i,item);    
    outString += item;
  }
  
  return(outString);
}

string
load::
readFirstString(ifstream & file)
{
  UInt num = leggiriga(file);
  assert(num >= 1);
  
  string val;
  assegnavalore(0,val);
  
  return(val);
}

UInt
load::
readFirstUInt(ifstream & file)
{
  UInt num = leggiriga(file);
  assert(num >= 1);
  
  UInt val;
  assegnavalore(0,val);
  
  return(val);
}

int
load::
readFirstInt(ifstream & file)
{
  UInt num = leggiriga(file);
  assert(num >= 1);
  
  UInt val;
  assegnavalore(0,val);
  
  return(int(val));
}

Real
load::
readFirstReal(ifstream & file)
{
  UInt num = leggiriga(file);
  assert(num >= 1);
  
  Real val;
  assegnavalore(0,val);
  
  return(val);
}
