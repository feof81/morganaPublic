/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTEGRATORFHY3D_STD_HPP
#define INTEGRATORFHY3D_STD_HPP

#include "intPolicySTD.hpp"
#include "elCard3d.hpp"
#include "intStaticCards.h"
#include "geoMapInterface.hpp"

#include "morganaIntegrator.hpp"
#include "morganaFiniteElements.hpp"


/*! Standard 3d integrator (functional) */
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
class integratorFHY3d_STD
{
  /*! @name Typedefs */ //@{
  public:        
    typedef typename FUNCTIONAL::TEST_FECARD   TEST_FECARD;
    typedef typename FUNCTIONAL::TEST_GEOSHAPE TEST_GEOSHAPE;
    typedef typename FUNCTIONAL::TEST_PMAPTYPE TEST_PMAPTYPE;
    
    typedef typename FUNCTIONAL::INTGRID  INTGRID;
    typedef typename INTGRID::GRID_ELMAP  PMAPTYPE;
    
    typedef typename INTGRID::GEOSHAPE3D   GEOSHAPE3D;
    typedef typename GEOSHAPE3D::GEOBSHAPE GEOSHAPE2D;
    
    typedef elCard3d<GEOSHAPE3D,PMAPTYPE> ELCARD;
    typedef intPolicySTD                  INTCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const intClass intFlag = FHY3d_STD;
    static const OPType   opFlag  = fnHY3d;
    
    geoMapInterface<GEOSHAPE2D> geoInterface;
    //@}
    
    /*! @name Functions */ //@{
  public:
    /*! Constructor */
    integratorFHY3d_STD();
    
    /*! Constructor */
    integratorFHY3d_STD(const INTCARD & IntCard);
    
    /*! Set integration card */
    void setIntCard(const INTCARD & IntCard);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(FUNCTIONAL & op);
    
    /*! The variable \c mat should be zero */
    sVect<Real> integration(const Teuchos::RCP<FUNCTIONAL> & op);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFHY3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integratorFHY3d_STD()
{
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
integratorFHY3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integratorFHY3d_STD(const INTCARD & IntCard)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim  == 3>::returnValue);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
void
integratorFHY3d_STD<FUNCTIONAL,TYPE,PRECISION>::
setIntCard(const INTCARD & IntCard)
{
  assert( (&IntCard) == (&IntCard) );
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFHY3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integration(FUNCTIONAL & op)
{
  typedef intStaticCard<GEOSHAPE2D,TYPE,PRECISION> INTSTATIC;
  
  //Alloc
  UInt matSize = op.numIndex_test();
  sVect<Real> bufMat(matSize);
  sVect<Real> mat(matSize);
  point3d Vx, Vy, N;
  
  //Integration loop
  for(UInt j=1; j <= GEOSHAPE3D::numFaces; ++j)
  {
    op.setLocFace(j);
    
    for(UInt i=1; i <= INTSTATIC::N; ++i)
    {
      //Evaluate field
      op.eval(INTSTATIC::getYn(i), bufMat);
      
      //Evaluate volume element
      Vx = geoInterface.getDerX(op.getGlobFacePoints(), INTSTATIC::getYn(i));
      Vy = geoInterface.getDerY(op.getGlobFacePoints(), INTSTATIC::getYn(i));
      N  = Vx ^ Vy; 
  
      //Sum
      for(UInt k=1; k <= matSize; ++k)
      {
        mat(k) += bufMat(k) * N.norm2() * INTSTATIC::getWn(i);
      }
    }
  }
  
  return(mat);
}

template<typename FUNCTIONAL, intTypes TYPE, UInt PRECISION>
sVect<Real>
integratorFHY3d_STD<FUNCTIONAL,TYPE,PRECISION>::
integration(const Teuchos::RCP<FUNCTIONAL> & op)
{
  return(integration(*op));
}

#endif
