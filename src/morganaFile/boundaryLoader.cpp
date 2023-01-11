/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "boundaryLoader.h"

boundaryLoader::
boundaryLoader()
{
}

void
boundaryLoader::
loadCards(const string & fileName)
{
  //File open
  ifstream file(fileName.c_str());
  
  
  //File check
  if(!file.good())
  { cout << "ERROR - Boundary Loader: File not found" << endl; }
  
  assert(file.good());
  
  
  //Loading number of conditions
  UInt numData = load::leggiriga(file);
  UInt numBoundary;
  
  if(numData == 1)
  { load::assegnavalore(0,numBoundary); }
  else
  { cout << "ERROR - Boundary Loader: the first row must conaint only the number of boundary surfaces" << endl; }
  
  feCards.resize(numBoundary);
  
  
  //Loading the boundary conditions
  UInt geoId, bcFlag;
  
  for(UInt i=1; i <= numBoundary; ++i)
  {
    numData = load::leggiriga(file);
    
    if(numData == 2)
    {
      load::assegnavalore(0,geoId);
      load::assegnavalore(1,bcFlag);
      
      feCards(i).setGeoId(geoId);
      feCards(i).setBcFlag(bcFlag);
    }
    else
    { cout << "ERROR - Boundary Loader: each boundary condition must contain only the geoId and the bcFlag" << endl; } 
  }
  
  file.close();
}

sVect<bcCardFE>
boundaryLoader::
getCards()
{
  return(feCards);
}

sVect<UInt>
boundaryLoader::
getGeoIds(const UInt & bcFlag)
{
  sVect<UInt> outVect;
  
  for(UInt i=1; i <= feCards.size(); ++i)
  {
    if(feCards(i).getBcFlag() == bcFlag)
    { outVect.push_back(feCards(i).getGeoId()); }
  }
  
  return(outVect);
}