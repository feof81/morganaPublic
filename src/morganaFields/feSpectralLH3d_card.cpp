/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "feSpectralLH3d_card.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
feSpectralLH3d_card::
feSpectralLH3d_card()
{
  isActive = true;
  rx = 0;
  ry = 0;
  rz = 0;
}

feSpectralLH3d_card::
feSpectralLH3d_card(const feSpectralLH3d_card & C)
{
  isActive = C.isActive;
  rx       = C.rx;
  ry       = C.ry;
  rz       = C.rz;
}

feSpectralLH3d_card
feSpectralLH3d_card::
operator=(const feSpectralLH3d_card & C)
{
  isActive = C.isActive;
  rx       = C.rx;
  ry       = C.ry;
  rz       = C.rz;
  
  return(*this);
}

bool
feSpectralLH3d_card::
operator!=(const feSpectralLH3d_card & C) const
{
  bool flag = true;
  
  flag = flag & (isActive == C.isActive);
  flag = flag & (rx == C.rx);
  flag = flag & (ry == C.ry);
  flag = flag & (rz == C.rz);
  
  return(!flag);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
feSpectralLH3d_card::
setIsActive(const bool & IsActive)
{
  isActive = IsActive;
}

void
feSpectralLH3d_card::
setR(const UInt & RX, const UInt & RY, const UInt & RZ)
{
  rx = RX;
  ry = RY;
  rz = RZ;
}


//_________________________________________________________________________________________________
// GET CONST REFERENCE
//-------------------------------------------------------------------------------------------------
const bool &
feSpectralLH3d_card::
getIsActive() const
{
  return(isActive);
}

const UInt &
feSpectralLH3d_card::
getRx() const
{
  return(rx);
}

const UInt &
feSpectralLH3d_card::
getRy() const
{
  return(ry);
}

const UInt &
feSpectralLH3d_card::
getRz() const
{
  return(rz);
}



//_________________________________________________________________________________________________
// GET REFERENCE
//-------------------------------------------------------------------------------------------------
bool &
feSpectralLH3d_card::
getIsActive()
{
  return(isActive);
}

UInt &
feSpectralLH3d_card::
getRx()
{
  return(rx);
}

UInt &
feSpectralLH3d_card::
getRy()
{
  return(ry);
}

UInt &
feSpectralLH3d_card::
getRz()
{
  return(rz);
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feSpectralLH3d_card & V)
{
  f << "isActive : " << V.isActive << endl;
  f << "rx       : " << V.rx << endl;
  f << "ry       : " << V.ry << endl;
  f << "rz       : " << V.rz << endl;
  
  return(f);
}

