/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTERPOLATORFVFLUX3D_HPP
#define INTERPOLATORFVFLUX3D_HPP

#include "geoMapInterface.hpp"

#include "feStaticField3d.hpp"
#include "traitsInterpolator.hpp"

#include "traitsFvEl3d.hpp"
#include "fvEL_dofMapAdapter.h"


/*! Interpolation class: specialized form, SOURCEFIELD -> flux Finite Volume field */
template<typename SOURCEFIELD>
class interpolatorFvFlux3d
{
    /*! @name Typedefs - Source */ //@{
  public:
    typedef typename SOURCEFIELD::GEOSHAPE     SOURCE_GEOSHAPE;
    typedef typename SOURCEFIELD::OUTTYPE      SOURCE_OUTTYPE;
    typedef typename SOURCEFIELD::FIELD_FETYPE SOURCE_FETYPE;
    typedef typename SOURCEFIELD::PMAPTYPE     SOURCE_PMAPTYPE;
    typedef typename SOURCEFIELD::FEINTERFACE  SOURCE_FEINTERFACE;
    //@}
    
    /*! @name Typedefs - Geometry */ //@{
  public:
    typedef typename interpolatorTrait<SOURCEFIELD>::MESH    MESH;
    typedef typename interpolatorTrait<SOURCEFIELD>::CONNECT CONNECT;
    //@}
    
    /*! @name Typedefs - Target */ //@{
  public:
    typedef typename traitsFvEl3d<SOURCE_GEOSHAPE,point3d>::FLUX_FIELD TARGETFIELD;
    typedef typename TARGETFIELD::GEOSHAPE     TARGET_GEOSHAPE;
    typedef typename TARGETFIELD::OUTTYPE      TARGET_OUTTYPE;
    typedef typename TARGETFIELD::FIELD_FETYPE TARGET_FETYPE;
    typedef typename TARGETFIELD::PMAPTYPE     TARGET_PMAPTYPE;
    typedef typename TARGETFIELD::FEINTERFACE  TARGET_FEINTERFACE;
    //@}
    
    /*! @name Internal links */ //@{
  public:
    Teuchos::RCP<communicator> commDev;
    bool commLoaded;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    interpolatorFvFlux3d();
    interpolatorFvFlux3d(const Teuchos::RCP<communicator> & CommDev);
    interpolatorFvFlux3d(communicator & CommDev);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Computing functions */ //@{
  public:
    void interpolate(const Teuchos::RCP<SOURCEFIELD> & SourceField, const Teuchos::RCP<TARGETFIELD> & TargetField);
    void interpolate(const SOURCEFIELD & SourceField, const TARGETFIELD & TargetField);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename SOURCEFIELD>
interpolatorFvFlux3d<SOURCEFIELD>::
interpolatorFvFlux3d()
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = false;
}

template<typename SOURCEFIELD>
interpolatorFvFlux3d<SOURCEFIELD>::
interpolatorFvFlux3d(const Teuchos::RCP<communicator> & CommDev)
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = true;
}

template<typename SOURCEFIELD>
interpolatorFvFlux3d<SOURCEFIELD>::
interpolatorFvFlux3d(communicator & CommDev)
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = true;
}

template<typename SOURCEFIELD>
void
interpolatorFvFlux3d<SOURCEFIELD>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commLoaded = true;
  commDev    = CommDev;
}

template<typename SOURCEFIELD>
void
interpolatorFvFlux3d<SOURCEFIELD>::
setCommDev(communicator & CommDev)
{
  commLoaded = true;
  commDev    = Teuchos::rcpFromRef(CommDev);
}

template<typename SOURCEFIELD>
void
interpolatorFvFlux3d<SOURCEFIELD>::
interpolate(const Teuchos::RCP<SOURCEFIELD> & SourceField, const Teuchos::RCP<TARGETFIELD> & TargetField)
{
  assert(commLoaded);
  interpolate(*SourceField,*TargetField);
}

template<typename SOURCEFIELD>
void
interpolatorFvFlux3d<SOURCEFIELD>::
interpolate(const SOURCEFIELD & SourceField, const TARGETFIELD & TargetField)
{
  typedef SOURCE_GEOSHAPE                     GEOSHAPE3D;
  typedef typename SOURCE_GEOSHAPE::GEOBSHAPE GEOSHAPE2D;
  
  //Assert-----------------------------------------------------------
  assert(commLoaded);
  assert(SourceField.getStartupOk());
  assert(TargetField.getStartupOk());
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<MESH>           grid = SourceField.getMesh();
  Teuchos::RCP<CONNECT> gridConnect = SourceField.getMeshConnect();
  
  //Adapter buildup--------------------------------------------------
  fvEL_dofMapAdapter adapter(SourceField.getDofMapper().getFaceIsActive(), SourceField.getDofMapper().getNewFaceLid());
  adapter.startup();
  
  //DofLoop----------------------------------------------------------
  UInt fc, newFc, pt;
  point3d Y, V, dof;
  
  for(UInt i=1; i <= grid->getNumElements(); ++i)
  {
    for(UInt j=1; j <= geoMapInterface<GEOSHAPE3D>::getNumFaces(); ++j)
    {
      fc = gridConnect->getElementToFace(i,j);
      
      if(adapter.getOldIsActive(fc))
      {
	//Compute the mean flux
	dof.set(0.0, 0.0, 0.0);
	
	for(UInt k=1; k <= geoMapInterface<GEOSHAPE2D>::getNumPoints(); ++k)
	{  
	  pt = geoMapInterface<GEOSHAPE3D>::faceToPoint(j,k);
	  Y  = geoMapInterface<GEOSHAPE3D>::getRefNodes(pt);
	  
	  SourceField.evalL(i,Y,V);
	  dof += V;
	}
	
	dof /= Real(geoMapInterface<GEOSHAPE2D>::getNumPoints());
	
	//Set dof
	newFc = adapter.getOldToNew(fc);
	TargetField.setDofL(newFc,dof);
      }
    }
  }
}



#endif
