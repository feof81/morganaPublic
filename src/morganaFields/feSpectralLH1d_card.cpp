/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "feSpectralLH1d_card.h"


//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
feSpectralLH1d_card::
feSpectralLH1d_card()
{
  isActive = true;
  rx = 0;
}

feSpectralLH1d_card::
feSpectralLH1d_card(const feSpectralLH1d_card & C)
{
  isActive = C.isActive;
  rx       = C.rx;
}

feSpectralLH1d_card
feSpectralLH1d_card::
operator=(const feSpectralLH1d_card & C)
{
  isActive = C.isActive;
  rx       = C.rx;
  
  return(*this);
}

bool
feSpectralLH1d_card::
operator!=(const feSpectralLH1d_card & C) const
{
  bool flag = true;
  
  flag = flag & (isActive == C.isActive);
  flag = flag & (rx == C.rx);
  
  return(!flag);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
feSpectralLH1d_card::
setIsActive(const bool & IsActive)
{
  isActive = IsActive;
}

void
feSpectralLH1d_card::
setR(const UInt & RX)
{
  rx = RX;
}


//_________________________________________________________________________________________________
// GET CONST REFERENCE
//-------------------------------------------------------------------------------------------------
const bool &
feSpectralLH1d_card::
getIsActive() const
{
  return(isActive);
}

const UInt &
feSpectralLH1d_card::
getRx() const
{
  return(rx);
}



//_________________________________________________________________________________________________
// GET REFERENCE
//-------------------------------------------------------------------------------------------------
bool &
feSpectralLH1d_card::
getIsActive()
{
  return(isActive);
}

UInt &
feSpectralLH1d_card::
getRx()
{
  return(rx);
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feSpectralLH1d_card & V)
{
  f << "isActive : " << V.isActive << endl;
  f << "rx       : " << V.rx << endl;
  
  return(f);
}


