/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "intPolicySTC.h"

intPolicySTC::
intPolicySTC()
{
  isActive = true;
}

intPolicySTC::
intPolicySTC(const bool & IsActive)
{
  isActive = IsActive;
}

intPolicySTC::
intPolicySTC(const intPolicySTC & Policy)
{
  isActive = Policy.isActive;
}

intPolicySTC
intPolicySTC::
operator=(const intPolicySTC & Policy)
{
  isActive = Policy.isActive;
  return(*this);
}

void
intPolicySTC::
setIsActive(const bool & IsActive)
{
  isActive = IsActive;
}

const bool &
intPolicySTC::
getIsActive() const
{
  return(isActive);
}

bool &
intPolicySTC::
getIsActive()
{
  return(isActive);
}

ostream & operator<<(ostream & f, const intPolicySTC & P)
{
  f << P.isActive << endl;
  return(f);
}

