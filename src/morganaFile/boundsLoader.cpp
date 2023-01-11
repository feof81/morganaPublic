/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "boundsLoader.h"

boundsLoader::
boundsLoader()
{
}

void
boundsLoader::
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
  
  bCards.resize(numBoundary);
  
  
  //Loading the boundary conditions
  UInt bcFlag, numGeoIds;
  Real realVal;
  
  for(UInt i=1; i <= numBoundary; ++i)
  {
    numData = load::leggiriga(file);
    
    if(numData >= 3)
    {
      load::assegnavalore(0,numGeoIds);
      bCards(i).resetGeoIds();
      
      assert(numGeoIds >= 1);
      
      for(UInt j=1; j <= numGeoIds; ++j)
      {
	load::assegnavalore(j,realVal);
	bCards(i).addGeoId(realVal);
      }
      
      load::assegnavalore(numData-2,bcFlag);
      bCards(i).setFlag(bcFlag);
      
      load::assegnavalore(numData-1,realVal);
      bCards(i).setValue(realVal);
    }
    else
    { cout << "ERROR - Boundary Loader: each boundary condition must contain only the geoId and the bcFlag" << endl; } 
  }
  
  file.close();
}

sVect<boundCard>
boundsLoader::
getCards()
{
  return(bCards);
}
