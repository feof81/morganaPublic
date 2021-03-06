/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef INTERPOLATORNODALSURFACE3D_HPP
#define INTERPOLATORNODALSURFACE3D_HPP

#include "traitsInterpolator.hpp"
#include "supportLagrange3d.hpp"


/*! Interpolator for the nodal static fields. The \c sourceField should be \c static and \c nodal since the finite element should
 implement the reference nodes. The \c targetField should only implement the \c eval method and the \c targetField should be a d-1 manyfold
 where d is the dymension of the \c sourceField.*/
template<typename SOURCEFIELD, typename TARGETFIELD>
class interpolatorNodalSurface3d
{
  /*! @name Typedefs */ //@{
  public:
    typedef typename SOURCEFIELD::GEOSHAPE     SOURCE_GEOSHAPE;
    typedef typename SOURCEFIELD::OUTTYPE      SOURCE_OUTTYPE;
    typedef typename SOURCEFIELD::FIELD_FETYPE SOURCE_FETYPE;
    typedef typename SOURCEFIELD::PMAPTYPE     SOURCE_PMAPTYPE;
    typedef typename SOURCEFIELD::FEINTERFACE  SOURCE_FEINTERFACE;
    typedef typename interpolatorTrait<SOURCEFIELD>::MESH    SOURCE_MESH;
    typedef typename interpolatorTrait<SOURCEFIELD>::CONNECT SOURCE_CONNECT;
    
    typedef typename TARGETFIELD::GEOSHAPE     TARGET_GEOSHAPE;
    typedef typename TARGETFIELD::OUTTYPE      TARGET_OUTTYPE;
    typedef typename TARGETFIELD::FIELD_FETYPE TARGET_FETYPE;
    typedef typename TARGETFIELD::PMAPTYPE     TARGET_PMAPTYPE;
    typedef typename TARGETFIELD::FEINTERFACE  TARGET_FEINTERFACE;
    typedef typename interpolatorTrait<TARGETFIELD>::MESH    TARGET_MESH;
    typedef typename interpolatorTrait<TARGETFIELD>::CONNECT TARGET_CONNECT;
    //@}
   
    /*! @name Internal links */ //@{
  public:
    Teuchos::RCP<communicator> commDev;
    bool commLoaded;
    //@}
  
    /*! @name Constructors and set functions */ //@{
  public:
    interpolatorNodalSurface3d();
    interpolatorNodalSurface3d(const Teuchos::RCP<communicator> & CommDev);
    interpolatorNodalSurface3d(communicator & CommDev);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    //@}
    
    /*! @name Computing functions */ //@{
  public:
    //WARNING the data should be parallel updated
    void interpolate(const Teuchos::RCP<SOURCEFIELD> & SourceField,
                     const Teuchos::RCP<TARGETFIELD> & TargetField);
    
    //WARNING the data should be parallel updated
    void interpolate(const SOURCEFIELD & SourceField,
                     const TARGETFIELD & TargetField);
    
    //WARNING the data should be parallel updated
    void interpolate(const Teuchos::RCP<SOURCEFIELD> & SourceField,
                     const Teuchos::RCP<TARGETFIELD> & TargetField,
                     std::set<UInt> activeGeoIds3d);
    
    //WARNING the data should be parallel updated
    void interpolate(const SOURCEFIELD & SourceField,
                     const TARGETFIELD & TargetField,
                     std::set<UInt> activeGeoIds3d);
    //@}
};


template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
interpolatorNodalSurface3d()
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  assert(staticAssert< TARGET_FETYPE::dim        == 2        >::returnValue);
  
  assert(staticAssert< SOURCE_FETYPE::dim == 3 >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = false;
}

template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
interpolatorNodalSurface3d(const Teuchos::RCP<communicator> & CommDev)
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  assert(staticAssert< TARGET_FETYPE::dim        == 2        >::returnValue);
  
  assert(staticAssert< SOURCE_FETYPE::dim == 3 >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = true;
  commDev    = CommDev;
}

template<typename SOURCEFIELD, typename TARGETFIELD>
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
interpolatorNodalSurface3d(communicator & CommDev)
{
  assert(staticAssert< TARGET_FETYPE::isNodal                >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feBaseType == scalar   >::returnValue);
  assert(staticAssert< TARGET_FETYPE::feClass    == feStatic >::returnValue);
  assert(staticAssert< TARGET_FETYPE::dim        == 2        >::returnValue);
  
  assert(staticAssert< SOURCE_FETYPE::dim == 3 >::returnValue);
  
  assert(staticAssert< traitsBasic<SOURCE_OUTTYPE>::myType  ==  traitsBasic<TARGET_OUTTYPE>::myType >::returnValue);
  assert(staticAssert< SOURCE_PMAPTYPE::parallelType        == TARGET_PMAPTYPE::parallelType>::returnValue);
  
  commLoaded = true;
  commDev    = Teuchos::rcpFromRef(CommDev);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commLoaded = true;
  commDev    = CommDev;
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
setCommDev(communicator & CommDev)
{
  commLoaded = true;
  commDev    = Teuchos::rcpFromRef(CommDev);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
interpolate(const Teuchos::RCP<SOURCEFIELD> & SourceField,
            const Teuchos::RCP<TARGETFIELD> & TargetField)
{
  assert(commLoaded);
  interpolate(*SourceField,*TargetField);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
interpolate(const SOURCEFIELD & SourceField,
            const TARGETFIELD & TargetField)
{
  //Assert-----------------------------------------------------------
  assert(commLoaded);
  assert(SourceField.getStartupOk());
  assert(TargetField.getStartupOk());
  
  //Typedefs---------------------------------------------------------
  typedef geoMapInterface<TARGET_GEOSHAPE> TARGET_GEOMAPINTERFACE;
  typedef geoMapSupport3d<SOURCE_GEOSHAPE> SOURCE_GEOMAPSUPPORT;
  typedef typename TARGETFIELD::DOFCARD    TARGET_DOFCARD;
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<TARGET_MESH>       targetGrid2d = TargetField.getMesh();
  Teuchos::RCP<SOURCE_MESH>       sourceGrid3d = SourceField.getMesh();
  Teuchos::RCP<TARGET_CONNECT> targetConnect2d = TargetField.getMeshConnect();
  Teuchos::RCP<SOURCE_CONNECT> sourceConnect3d = SourceField.getMeshConnect();
  
  supportLagrange3d<SOURCE_GEOSHAPE,SOURCE_PMAPTYPE> support(sourceGrid3d,sourceConnect3d,targetGrid2d,targetConnect2d);
  TARGET_GEOMAPINTERFACE targetMapInterface;
  SOURCE_GEOMAPSUPPORT   sourceMapSupport;
  TARGET_DOFCARD         targetDofCard;
  
  //Loop-------------------------------------------------------------
  UInt el2d, lid;
  point3d Y2d, Y3d, P;
  sVect<point3d> points3d, points2d;
  SOURCE_OUTTYPE val;
  
  
  for(UInt el3d=1; el3d <= sourceGrid3d->getNumElements(); ++el3d)
  {
    support.setElement3d(el3d);
    points3d = sourceGrid3d->getElementNodesL(el3d);
    sourceMapSupport.setPoints(points3d);
    
    for(UInt f=1; f <= SOURCE_GEOSHAPE::numFaces; ++f)
    {
      support.setLocalFace(f);
      
      if(support.isBoundary())
      {
	el2d     = support.getElement2d();
	points2d = targetGrid2d->getElementNodesL(el2d);
	
	for(UInt j=1; j <= TARGET_FETYPE::numBasis; ++j)
	{
	  //Extract the solution
	  Y2d = TARGET_FETYPE::getRefNode(j);
          P   = targetMapInterface.getPosition(points2d,Y2d);
	  
	  sourceMapSupport.mapGlobalToLocal(P,Y3d);
	  SourceField.evalL(el3d,Y3d,val);
	  
	  //Indentify the lid
	  targetDofCard = TARGET_FETYPE::getDofCard(j);
          targetDofCard.setLocalElId(el2d);
	  
	  //Insert
	  if(TargetField.getDofMapper().isActive(targetDofCard))
	  {
	    lid = TargetField.getDofMapper().mapDofL(targetDofCard);
	    TargetField.setDofL(lid,val);
	  }
	}
      }
    }
  }
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
interpolate(const Teuchos::RCP<SOURCEFIELD> & SourceField,
            const Teuchos::RCP<TARGETFIELD> & TargetField,
            std::set<UInt> activeGeoIds3d)
{
  assert(commLoaded);
  interpolate(*SourceField,*TargetField,activeGeoIds3d);
}

template<typename SOURCEFIELD, typename TARGETFIELD>
void
interpolatorNodalSurface3d<SOURCEFIELD,TARGETFIELD>::
interpolate(const SOURCEFIELD & SourceField,
            const TARGETFIELD & TargetField,
            std::set<UInt> activeGeoIds3d)
{
  //Assert-----------------------------------------------------------
  assert(commLoaded);
  assert(SourceField.getStartupOk());
  assert(TargetField.getStartupOk());
  
  //Typedefs---------------------------------------------------------
  typedef geoMapInterface<TARGET_GEOSHAPE> TARGET_GEOMAPINTERFACE;
  typedef geoMapSupport3d<SOURCE_GEOSHAPE> SOURCE_GEOMAPSUPPORT;
  typedef typename TARGETFIELD::DOFCARD    TARGET_DOFCARD;
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<TARGET_MESH>       targetGrid2d = TargetField.getMesh();
  Teuchos::RCP<SOURCE_MESH>       sourceGrid3d = SourceField.getMesh();
  Teuchos::RCP<TARGET_CONNECT> targetConnect2d = TargetField.getMeshConnect();
  Teuchos::RCP<SOURCE_CONNECT> sourceConnect3d = SourceField.getMeshConnect();
  
  supportLagrange3d<SOURCE_GEOSHAPE,SOURCE_PMAPTYPE> support(sourceGrid3d,sourceConnect3d,targetGrid2d,targetConnect2d);
  TARGET_GEOMAPINTERFACE targetMapInterface;
  SOURCE_GEOMAPSUPPORT   sourceMapSupport;
  TARGET_DOFCARD         targetDofCard;
  
  //Loop-------------------------------------------------------------
  UInt el2d, lid, geoId;
  point3d Y2d, Y3d, P;
  sVect<point3d> points3d, points2d;
  SOURCE_OUTTYPE val;
  
  
  for(UInt el3d=1; el3d <= sourceGrid3d->getNumElements(); ++el3d)
  {
    geoId = sourceGrid3d->getElementL(el3d).getGeoId();
    
    if(activeGeoIds3d.count(geoId) >= 1)
    {
      support.setElement3d(el3d);
      points3d = sourceGrid3d->getElementNodesL(el3d);
      sourceMapSupport.setPoints(points3d);
      
      //Face loop
      for(UInt f=1; f <= SOURCE_GEOSHAPE::numFaces; ++f)
      {
	support.setLocalFace(f);
	
	if(support.isBoundary())
	{
	  el2d     = support.getElement2d();
	  points2d = targetGrid2d->getElementNodesL(el2d);
	  
	  for(UInt j=1; j <= TARGET_FETYPE::numBasis; ++j)
	  {
	    //Extract the solution
	    Y2d = TARGET_FETYPE::getRefNode(j);
	    P   = targetMapInterface.getPosition(points2d,Y2d);
	    
	    sourceMapSupport.mapGlobalToLocal(P,Y3d);
	    SourceField.evalL(el3d,Y3d,val);
	    
	    //Indentify the lid
	    targetDofCard = TARGET_FETYPE::getDofCard(j);
	    targetDofCard.setLocalElId(el2d);

	    //Insert
	    if(TargetField.getDofMapper().isActive(targetDofCard))
	    {
	      lid = TargetField.getDofMapper().mapDofL(targetDofCard);
	      TargetField.setDofL(lid,val);
	    }
	  }
	}
      } //End face loop
      
    }
  }
  
  //Parallel update--------------------------------------------------
}


#endif
