/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEXFEM_LT3D_CARD_HPP
#define FEXFEM_LT3D_CARD_HPP

#include "morganaTypes.hpp"

#include "geoMapInterface.hpp"
#include "geoMapSupport3d.hpp"

#include "morganaFields.hpp"
#include "fePr3d.hpp"
#include "elCard3d.hpp"
#include "feDynamicDofCard3d.h"


/*! The information carrying class for the XFEM-FE. Contains the local ids of the enrichment-dofs,
 the four dofs that represent the level set function, a logic that states whether the element is 
 active and the variable xtoll: represent the minimum volume fraction of wich the element is 
 considered splitted. */
class feXFEM_LT3d_card
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
  
    /*! @name Internal data */ //@{
  public:
    static const FECardLabel cardLabel = CR_XF_3d;  
    bool isActive;
    sVect<UInt> enrichedDof;
    sVect<Real> phi;
    Real xToll;
    //@}
  
    /*! @name Constructor and functions */ //@{
  public:
    feXFEM_LT3d_card();
    feXFEM_LT3d_card(const feXFEM_LT3d_card & C);
    feXFEM_LT3d_card operator=(const feXFEM_LT3d_card & C);
    bool operator!=(const feXFEM_LT3d_card & C) const;
    //@}
   
    /*! @name Set functions */ //@{
  public:
    void setActive(const bool & IsActive);
    void setXtoll(const Real & XToll);
    void setPhi(const sVect<Real> & Phi);
    void setEnrichedDof(const sVect<UInt> & EnrichedDof);
    void resizeEnrichedDof(const UInt & n);
    //@}
    
    /*! @name Get Phi */ //@{
  public:
    sVect<Real>       & getPhi();
    const sVect<Real> & getPhi() const;
    Real              & getPhi(const UInt & i);
    const Real        & getPhi(const UInt & i) const;
    //@}
    
    /*! @name Get Active dof */ //@{
  public:
    sVect<UInt>       & getEnrichedDof();
    const sVect<UInt> & getEnrichedDof() const;
    UInt              & getEnrichedDof(const UInt & i);
    const UInt        & getEnrichedDof(const UInt & i) const;
    UInt                getNumEnrichedBasis() const;
    //@}
    
    /*! @name Get Xtoll */ //@{
  public:
    Real       & getXtoll();
    const Real & getXtoll() const;
    //@}
    
    /*! @name Get Active */ //@{
  public:
    bool       & getIsActive();
    const bool & getIsActive() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    friend ostream & operator<<(ostream & f, const feXFEM_LT3d_card & V);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
feXFEM_LT3d_card::
feXFEM_LT3d_card()
{
  xToll    = geoToll;
  isActive = true;
  
  phi.resize(4);
}

feXFEM_LT3d_card::
feXFEM_LT3d_card(const feXFEM_LT3d_card & C)
{
  xToll     = C.xToll;
  isActive  = C.isActive;
  phi       = C.phi;
  enrichedDof = C.enrichedDof;
}

feXFEM_LT3d_card
feXFEM_LT3d_card::
operator=(const feXFEM_LT3d_card & C)
{
  xToll     = C.xToll;
  isActive  = C.isActive;
  phi       = C.phi;
  enrichedDof = C.enrichedDof;
  
  return(*this);
}

bool
feXFEM_LT3d_card::
operator!=(const feXFEM_LT3d_card & C) const
{
  bool flag = true;
  
  flag = flag & (isActive == C.isActive);
  
  for(UInt i=1; i <= enrichedDof.size(); ++i)
  { flag = flag & (enrichedDof(i) == C.enrichedDof(i)); }
  
  for(UInt i=1; i <= phi.size(); ++i)
  { flag = flag & (phi(i) == C.phi(i)); }
  
  flag = flag & (xToll == C.xToll);
  
  return(!flag);
}


//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
feXFEM_LT3d_card::
setPhi(const sVect<Real> & Phi)
{
  phi = Phi;
}

void
feXFEM_LT3d_card::
setActive(const bool & IsActive)
{
  isActive = IsActive;
}

void
feXFEM_LT3d_card::
setEnrichedDof(const sVect<UInt> & EnrichedDof)
{
  enrichedDof = EnrichedDof;
}

void
feXFEM_LT3d_card::
resizeEnrichedDof(const UInt & n)
{
  enrichedDof.resize(n);
}

void
feXFEM_LT3d_card::
setXtoll(const Real & XToll)
{
  xToll = XToll;
}



//_________________________________________________________________________________________________
// GET PHI
//-------------------------------------------------------------------------------------------------
sVect<Real> &
feXFEM_LT3d_card::
getPhi()
{
  return(phi);
}

const sVect<Real> &
feXFEM_LT3d_card::
getPhi() const
{
  return(phi);
}

Real &
feXFEM_LT3d_card::
getPhi(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  return(phi(i));
}

const Real &
feXFEM_LT3d_card::
getPhi(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= 4);
  return(phi(i));
}


//_________________________________________________________________________________________________
// GET ENRICHED
//-------------------------------------------------------------------------------------------------
sVect<UInt> &
feXFEM_LT3d_card::
getEnrichedDof()
{
  return(enrichedDof);
}

const sVect<UInt> &
feXFEM_LT3d_card::
getEnrichedDof() const
{
  return(enrichedDof);
}

UInt &
feXFEM_LT3d_card::
getEnrichedDof(const UInt & i)
{
  assert(i >= 1);
  assert(i <= enrichedDof.size());
  
  return(enrichedDof(i));
}

const UInt &
feXFEM_LT3d_card::
getEnrichedDof(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= enrichedDof.size());
  
  return(enrichedDof(i));
}

UInt
feXFEM_LT3d_card::
getNumEnrichedBasis() const
{
  return(enrichedDof.size());
}



//_________________________________________________________________________________________________
// GET XTOLL
//-------------------------------------------------------------------------------------------------
Real &
feXFEM_LT3d_card::
getXtoll()
{
  return(xToll);
}

const Real &
feXFEM_LT3d_card::
getXtoll() const
{
  return(xToll);
}


//_________________________________________________________________________________________________
// GET ISACTIVE
//-------------------------------------------------------------------------------------------------
bool &
feXFEM_LT3d_card::
getIsActive()
{
  return(isActive);
}

const bool &
feXFEM_LT3d_card::
getIsActive() const
{
  return(isActive);
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
ostream & operator<<(ostream & f, const feXFEM_LT3d_card & V)
{
  f << "isActive : " << V.isActive << endl;
  f << "xToll    : " << V.xToll << endl;
  f << "enriched : ";
  
  for(UInt i=1; i <= V.enrichedDof.size(); ++i)
  {
    f << " " << V.enrichedDof(i);
  }
  f << endl;
  
  f << "phi      : ";
  
  for(UInt i=1; i <= V.phi.size(); ++i)
  {
    f << " " << V.phi(i);
  }
  f << endl;
  
  return(f);
}



//_________________________________________________________________________________________________
// PARALLEL SUPPORT
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
feXFEM_LT3d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & isActive;
  ar & enrichedDof;
  ar & phi;
  ar & xToll;
}



#endif
