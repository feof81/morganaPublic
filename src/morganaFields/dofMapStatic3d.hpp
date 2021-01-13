/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef DOFMAPSTATIC3D_HPP
#define DOFMAPSTATIC3D_HPP

#include "typesInterface.hpp"
#include "traitsBasic.h"

#include "pVectManip.hpp"

#include "morganaGeometry.hpp"
#include "geoMapInterface.hpp"
#include "connect3d.hpp"
#include "mesh3dGlobalManip.hpp"

#include "dofMapStatic3d_options.h"
#include "feStaticDofCard3d.h"



enum dms3d_order {dms3d_vectMajor, dms3d_componentMajor};
enum dms3d_mode  {dms3d_allMode, dms3d_geoIdMode, dms3d_disjointMode};




//! Forward declaration
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER = dms3d_vectMajor, dms3d_mode MODE = dms3d_allMode> class dofMapStatic3d;


//_________________________________________________________________________________________________
// ALL MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 3d static finite elements -> \c allMode */
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
class dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic3d_options                OPTIONS;
    
    typedef pMap<PMAPTYPE> DOFMAP;
    typedef pMap<PMAPTYPE> LISTMAP;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    FETYPE          refFe;
    GEOMAPINTERFACE refShape;
    FECARDS         feCards;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator>  commDev;
    Teuchos::RCP<MESH3D>        grid3d;
    Teuchos::RCP<CONNECT3D>     connectGrid3d;
    Teuchos::RCP<OPTIONS>       options;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool commDevLoaded, geometryLoaded, optionsLoaded, startupOk;
    //@}
    
     /*! @name Static data */ //@{
  public:
    static const UInt block = traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
    //@}
    
    /*! @name Internal Numbers */ //@{
  public:
    UInt numVerticesL, numVerticesG;
    UInt numEdgesL,    numEdgesG;
    UInt numFacesL,    numFacesG;
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
    UInt numEdgesDofsL,  numEdgesDofsG;
    UInt numFacesDofsL,  numFacesDofsG;
    UInt numVolumeDofsL, numVolumeDofsG;
    
    UInt numDofsL, numDofsG;
    UInt sizeListL, sizeListG;
    //@}
    
    /*! @name Internal Maps */ //@{
  public:
    DOFMAP  dofMap;
    LISTMAP listMap;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    dofMapStatic3d();
    dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    dofMapStatic3d(const dofMapStatic3d & DofMap);
    dofMapStatic3d operator=(const dofMapStatic3d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    //@}
    
    /*! @name Options loading and startup */ //@{
  public:
    void setFeCardL(const UInt & lid, const FECARD & FeCards);
    void setFeCardG(const UInt & gid, const FECARD & FeCards);
    void setFeCards(const FECARDS & FECards);
    void setOptions(const Teuchos::RCP<OPTIONS> & Options);
    void setOptions(OPTIONS & Options);
    void startup();
    //@}
    
    /*! @name Get functions */ //@{
  public:
    bool isActive(const DOFCARD & dofCard) const;
    UInt dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const;
    UInt dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const;
    UInt mapDofL(const DOFCARD & dofCard) const;
    UInt mapDofG(const DOFCARD & dofCard) const;    
    UInt mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    UInt mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    const UInt & getNumDofsL() const;
    const UInt & getNumDofsG() const;
    const UInt & getSizeListL() const;
    const UInt & getSizeListG() const;
    //@}
    
     /*! @name Dump functions */ //@{
  public:
    const DOFMAP                     & getDofMap()  const;
    const LISTMAP                    & getListMap() const;
    const FECARDS                    & getFeCards() const;
          FECARDS                    & getFeCards();
    const Teuchos::RCP<communicator> & getCommDev() const;
    const Teuchos::RCP<MESH3D>       & getGrid3d()  const;
    const Teuchos::RCP<CONNECT3D>    & getConnectGrid3d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms3d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_allMode> & V);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid3d        = Grid3d;
  connectGrid3d = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofMapStatic3d(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
operator=(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;

  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev        = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev        = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  geometryLoaded = true;
  grid3d         = Grid3d;
  connectGrid3d  = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  geometryLoaded = true; 
  grid3d         = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}



//_________________________________________________________________________________________________
// OPTIONS LOADING AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
startup()
{
  //Assert and flags_______________________________________________________________________________
  assert(commDevLoaded);
  assert(geometryLoaded);
  assert(optionsLoaded);
  
  assert(grid3d->getElements().colIsLocal());
  assert(grid3d->getFaces().colIsLocal());
  assert(grid3d->getEdges().colIsLocal());
  
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid3d));
  
  startupOk = true;
  
  
  //Count the geotypes_____________________________________________________________________________ 
  numVerticesL = grid3d->getNumVertices();
  numEdgesL    = grid3d->getNumEdges();
  numFacesL    = grid3d->getNumFaces();
  numElementsL = grid3d->getNumElements();
  
  numVerticesG = manipulator.getNumGlobalVertices(grid3d);
  numEdgesG    = manipulator.getNumGlobalEdges(grid3d);
  numFacesG    = manipulator.getNumGlobalFaces(grid3d);
  numElementsG = manipulator.getNumGlobalElements(grid3d);
  
  //Count the dofs
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numEdgesDofsL  = FETYPE::dofPerEdge   * numEdgesL;
  numFacesDofsL  = FETYPE::dofPerFace   * numFacesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numEdgesDofsG  = FETYPE::dofPerEdge   * numEdgesG;
  numFacesDofsG  = FETYPE::dofPerFace   * numFacesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numEdgesDofsL + numFacesDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numEdgesDofsG + numFacesDofsG + numVolumeDofsG;
  
  sizeListL = numDofsL * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  sizeListG = numDofsG * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  
  
  //Build dof map__________________________________________________________________________________
  PMAPTYPE mapItem;
  UInt k=1;
  
  dofMap.resize(numDofsL);
  
  for(UInt lev=1; lev <= FETYPE::dofPerVertex; ++lev)  //Vertices
  {
    for(UInt i=1; i <= numVerticesL; ++i)
    {
      mapItem = grid3d->getNodes().getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numVerticesG) );
      
      dofMap(k) = mapItem;
      ++k;
    } 
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerEdge; ++lev)  //Edges
  {
    for(UInt i=1; i <= numEdgesL; ++i)
    {
      mapItem = grid3d->getEdges().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numEdgesG) + numVertexDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerFace; ++lev)  //Faces
  {
    for(UInt i=1; i <= numFacesL; ++i)
    {
      mapItem = grid3d->getFaces().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numFacesG) + numVertexDofsG + numEdgesDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL; ++i)
    {
      mapItem = grid3d->getElements().getRowMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG + numEdgesDofsG + numFacesDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  
  //Build list map_________________________________________________________________________________
  listMap.resize(sizeListL);
  
  UInt lid, gid;
  
  for(UInt i=1; i <= dofMap.size(); ++i)
  {   
    for(UInt K=1; K <= block; ++K)
    {
      mapItem = dofMap(i);
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor);
      
      mapItem.setLid(lid);
      mapItem.setGid(gid);
      
      listMap(lid) = mapItem;
    }
  }
  
  
  //Testing________________________________________________________________________________________
  pMapGlobalManip<PMAPTYPE> tester(commDev); 
  assert(tester.check(dofMap));
  assert(tester.check(listMap));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
bool
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
isActive(const DOFCARD & dofCard) const
{
  assert( (&dofCard) == (&dofCard) );
  return(true);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms3d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms3d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());     
      return( grid3d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      return( connectGrid3d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numEdgesL  + numVertexDofsL );
    
    case FACE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerFace);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numFaces);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      return( connectGrid3d->getElementToFace(dofCard.getLocalElId(), dofCard.getLocalId()) + (dofCard.getLevel() - 1) * numFacesL  + numVertexDofsL + numEdgesDofsL );
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      return(dofCard.getLocalElId() + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL + numEdgesDofsL + numFacesDofsL );
      
    default :
      assert(1==2);
  }
  
  return(0);             
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::DOFMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getDofMap() const
{
  assert(startupOk);
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::LISTMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getListMap() const
{
  assert(startupOk);
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getNumDofsL() const
{
  assert(startupOk);
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getNumDofsG() const
{
  assert(startupOk);
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getSizeListL() const
{
  assert(startupOk);
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getSizeListG() const
{
  assert(startupOk);
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::MESH3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getGrid3d() const
{
  assert(geometryLoaded);
  return(grid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::CONNECT3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getConnectGrid3d() const
{
  assert(geometryLoaded);
  return(connectGrid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::OPTIONS> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const bool &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_allMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms3d_order O>
ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_allMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numEdges      : " << V.numEdgesL      << " " << V.numEdgesG << endl;
  f << "numFaces      : " << V.numFacesL      << " " << V.numFacesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numFacesDofs  : " << V.numFacesDofsL  << " " << V.numFacesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}



//_________________________________________________________________________________________________
// GEO ID MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 3d static finite elements -> \c geoIdMode */
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
class dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH3D;
    typedef connect3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT3D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapStatic3d_options                OPTIONS;
    
    typedef pMap<PMAPTYPE> DOFMAP;
    typedef pMap<PMAPTYPE> LISTMAP;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    FETYPE          refFe;
    GEOMAPINTERFACE refShape;
    FECARDS         feCards;
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<communicator>  commDev;
    Teuchos::RCP<MESH3D>        grid3d;
    Teuchos::RCP<CONNECT3D>     connectGrid3d;
    Teuchos::RCP<OPTIONS>       options;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool commDevLoaded, geometryLoaded, optionsLoaded, startupOk;
    //@}
    
     /*! @name Static data */ //@{
  public:
    static const UInt block = traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    sVect<bool> vertexIsActive;
    sVect<bool> edgeIsActive;
    sVect<bool> faceIsActive;
    sVect<bool> elementIsActive;
    
    sVect<UInt> newVertexLid;
    sVect<UInt> newEdgeLid;
    sVect<UInt> newFaceLid;
    sVect<UInt> newElementLid;
    
    pVect<UInt,PMAPTYPE>  globVertices;
    pVect<UInt,PMAPTYPE>  globEdges;
    pVect<UInt,PMAPTYPE>  globFaces;
    pVect<UInt,PMAPTYPE>  globElements;
    
    UInt numVerticesL, numVerticesG;
    UInt numEdgesL,    numEdgesG;
    UInt numFacesL,    numFacesG;
    UInt numElementsL, numElementsG;
    
    UInt numVertexDofsL, numVertexDofsG;
    UInt numEdgesDofsL,  numEdgesDofsG;
    UInt numFacesDofsL,  numFacesDofsG;
    UInt numVolumeDofsL, numVolumeDofsG;
    
    UInt numDofsL, numDofsG;
    UInt sizeListL, sizeListG;
    //@}
    
    /*! @name Internal Maps */ //@{
  public:
    DOFMAP  dofMap;
    LISTMAP listMap;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    dofMapStatic3d();
    dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    dofMapStatic3d(const dofMapStatic3d & DofMap);
    dofMapStatic3d operator=(const dofMapStatic3d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    void setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    //@}
    
    /*! @name Options loading and startup */ //@{
  public:
    void setFeCardL(const UInt & lid, const FECARD & FeCards);
    void setFeCardG(const UInt & gid, const FECARD & FeCards);
    void setFeCards(const FECARDS & FECards);
    void setOptions(const Teuchos::RCP<OPTIONS> & Options);
    void setOptions(OPTIONS & Options);
    void startup();
    //@}
    
    /*! @name Get functions */ //@{
  public:
    bool isActive(const DOFCARD & dofCard) const;
    UInt dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const;
    UInt dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const;
    UInt mapDofL(const DOFCARD & dofCard) const;
    UInt mapDofG(const DOFCARD & dofCard) const;    
    UInt mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    UInt mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const;
    const UInt & getNumDofsL() const;
    const UInt & getNumDofsG() const;
    const UInt & getSizeListL() const;
    const UInt & getSizeListG() const;
    //@}
    
    /*! @name Dump functions */ //@{
  public:
    const DOFMAP                     & getDofMap() const;
    const LISTMAP                    & getListMap() const;
    const FECARDS                    & getFeCards() const;
          FECARDS                    & getFeCards();
    const Teuchos::RCP<communicator> & getCommDev() const;
    const Teuchos::RCP<MESH3D>       & getGrid3d() const;
    const Teuchos::RCP<CONNECT3D>    & getConnectGrid3d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Extra dump functions */ //@{
  public:
    const sVect<bool> & getVertexIsActive() const;
    const sVect<bool> & getEdgeIsActive() const;
    const sVect<bool> & getFaceIsActive() const;
    const sVect<bool> & getElementIsActive() const;
    
    const sVect<UInt> & getNewVertexLid() const;
    const sVect<UInt> & getNewEdgeLid() const;
    const sVect<UInt> & getNewFaceLid() const;
    const sVect<UInt> & getNewElementLid() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dms3d_order O>
    friend ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_allMode> & V);
    //@}
};



template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid3d        = Grid3d;
  connectGrid3d = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid3d        = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofMapStatic3d(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid; 
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
operator=(const dofMapStatic3d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid3d        = DofMap.grid3d;
  connectGrid3d = DofMap.connectGrid3d;
  options       = DofMap.options;

  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
     
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
   
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid; 
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
    
  numVerticesL = DofMap.numVerticesL;
  numVerticesG = DofMap.numVerticesG;
  numEdgesL    = DofMap.numEdgesL;
  numEdgesG    = DofMap.numEdgesG;
  numFacesL    = DofMap.numFacesL;
  numFacesG    = DofMap.numFacesG;
  numElementsL = DofMap.numElementsL;
  numElementsG = DofMap.numElementsG;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numEdgesDofsL  = DofMap.numEdgesDofsL;
  numEdgesDofsG  = DofMap.numEdgesDofsG;
  numFacesDofsL  = DofMap.numFacesDofsL;
  numFacesDofsG  = DofMap.numFacesDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;
    
  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  geometryLoaded = true;
  grid3d         = Grid3d;
  connectGrid3d  = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setGeometry(MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
{
  geometryLoaded = true; 
  grid3d         = Teuchos::rcp(new MESH3D(Grid3d));
  connectGrid3d  = Teuchos::rcp(new CONNECT3D(ConnedGrid3d));
  
  feCards.resize(grid3d->getNumElements());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
void
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
startup()
{
  assert(commDevLoaded);
  assert(optionsLoaded);
  
  assert(grid3d->getElements().colIsLocal());
  assert(grid3d->getFaces().colIsLocal());
  assert(grid3d->getEdges().colIsLocal());
  
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid3d));
  
  startupOk = true;
  
  //Alloc__________________________________________________________________________________________
  vertexIsActive.resize(  grid3d->getNumVertices() );
  edgeIsActive.resize(    grid3d->getNumEdges()    );
  faceIsActive.resize(    grid3d->getNumFaces()    );
  elementIsActive.resize( grid3d->getNumElements() );
   
  newVertexLid.resize(  grid3d->getNumVertices() );
  newEdgeLid.resize(    grid3d->getNumEdges()    );
  newFaceLid.resize(    grid3d->getNumFaces()    );
  newElementLid.resize( grid3d->getNumElements() );
  
  for(UInt i=1; i <= grid3d->getNumVertices(); ++i)
  { vertexIsActive(i) = false; }
  
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
  { edgeIsActive(i) = false; }
  
  for(UInt i=1; i <= grid3d->getNumFaces(); ++i)
  { faceIsActive(i) = false;}
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  { elementIsActive(i) = false; }

  
  //Identify the active items______________________________________________________________________
  UInt geoId;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    geoId = grid3d->getElementL(i).getGeoId();
    
    if(options->isGeoId(geoId))
    {
      for(UInt j=1; j <= GEOSHAPE::numVertices; ++j)
      {
	vertexIsActive(grid3d->getElementL(i).getCid(j)) = true;
      }
      
      for(UInt j=1; j <= GEOSHAPE::numEdges; ++j)
      {
	edgeIsActive(connectGrid3d->getElementToEdge(i,j)) = true;
      }
      
      for(UInt j=1; j <= GEOSHAPE::numFaces; ++j)
      {
	faceIsActive(connectGrid3d->getElementToFace(i,j)) = true;
      }
      
      elementIsActive(i) = true;
    }
  }
  
  //Assign new local ids___________________________________________________________________________
  numVerticesL = 0;
  numEdgesL    = 0;
  numFacesL    = 0;
  numElementsL = 0;
  
  globVertices.clear();
  globEdges.clear();
  globFaces.clear();
  globElements.clear();
  
  PMAPTYPE nodeMapItem;
  PMAPTYPE elMapItem;
  
  for(UInt i=1; i <= grid3d->getNumVertices(); ++i)  //Vertices
  {
    if(vertexIsActive(i))
    {
      ++numVerticesL;
      newVertexLid(i) = numVerticesL;
      
      nodeMapItem = grid3d->getNodes().getMapL(i);
      nodeMapItem.setLid(numVerticesL);
      
      globVertices.push_back(nodeMapItem, nodeMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid3d->getNumEdges(); ++i)     //Edges
  {
    if(edgeIsActive(i))
    {
      ++numEdgesL;
      newEdgeLid(i) = numEdgesL;
      
      elMapItem = grid3d->getEdges().getRowMapL(i);
      elMapItem.setLid(numEdgesL);
      
      globEdges.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid3d->getNumFaces(); ++i)     //Faces
  {
    if(faceIsActive(i))
    {
      ++numFacesL;
      newFaceLid(i) = numFacesL;
      
      elMapItem = grid3d->getFaces().getRowMapL(i);
      elMapItem.setLid(numFacesL);
      
      globFaces.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)  //Elements
  {
    if(elementIsActive(i))
    {
      ++numElementsL;
      newElementLid(i) = numElementsL;
      
      elMapItem = grid3d->getElements().getRowMapL(i);
      elMapItem.setLid(numElementsL);
      
      globElements.push_back(elMapItem, elMapItem.getGid());
    }
  }
  
  
  //Assign new global ids__________________________________________________________________________
  pVectGlobalManip<UInt,PMAPTYPE> nodeManip(commDev);
  pVectGlobalManip<UInt,PMAPTYPE>   elManip(commDev);
   
  nodeManip.buildGlobalNumbering(globVertices);
  elManip.buildGlobalNumbering(globEdges);
  elManip.buildGlobalNumbering(globFaces);
  elManip.buildGlobalNumbering(globElements);
  
  numVerticesG = nodeManip.sizeG(globVertices);
  numEdgesG    = elManip.sizeG(globEdges);
  numFacesG    = elManip.sizeG(globFaces);
  numElementsG = elManip.sizeG(globElements);
   
  //Count the dofs_________________________________________________________________________________
  numVertexDofsL = FETYPE::dofPerVertex * numVerticesL;
  numEdgesDofsL  = FETYPE::dofPerEdge   * numEdgesL;
  numFacesDofsL  = FETYPE::dofPerFace   * numFacesL;
  numVolumeDofsL = FETYPE::dofPerVolume * numElementsL;
  
  numVertexDofsG = FETYPE::dofPerVertex * numVerticesG;
  numEdgesDofsG  = FETYPE::dofPerEdge   * numEdgesG;
  numFacesDofsG  = FETYPE::dofPerFace   * numFacesG;
  numVolumeDofsG = FETYPE::dofPerVolume * numElementsG;
  
  //Total number of dofs and lists_________________________________________________________________
  numDofsL = numVertexDofsL + numEdgesDofsL + numFacesDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numEdgesDofsG + numFacesDofsG + numVolumeDofsG;
  
  sizeListL = numDofsL * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  sizeListG = numDofsG * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  
  
  //Build dof map__________________________________________________________________________________
  PMAPTYPE mapItem;
  UInt k=1;
  
  dofMap.resize(numDofsL);
  
  for(UInt lev=1; lev <= FETYPE::dofPerVertex; ++lev)  //Vertices
  {
    for(UInt i=1; i <= numVerticesL; ++i)
    {
      mapItem = globVertices.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numVerticesG) );
      
      dofMap(k) = mapItem;
      ++k;
    } 
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerEdge; ++lev)  //Edges
  {
    for(UInt i=1; i <= numEdgesL; ++i)
    {
      mapItem = globEdges.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numEdgesG) + numVertexDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerFace; ++lev)  //Faces
  {
    for(UInt i=1; i <= numFacesL; ++i)
    {
      mapItem = globFaces.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numFacesG) + numVertexDofsG + numEdgesDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  for(UInt lev=1; lev <= FETYPE::dofPerVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL; ++i)
    {
      mapItem = globElements.getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + ((lev-1) * numElementsG) + numVertexDofsG + numEdgesDofsG + numFacesDofsG );
      
      dofMap(k) = mapItem;
      ++k;
    }
  }
  
  
  //Build list map_________________________________________________________________________________
  listMap.resize(sizeListL);
  
  UInt lid, gid;
  
  for(UInt i=1; i <= dofMap.size(); ++i)
  {   
    for(UInt K=1; K <= block; ++K)
    {
      mapItem = dofMap(i);
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dms3d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor);
      
      mapItem.setLid(lid);
      mapItem.setGid(gid);
      
      listMap(lid) = mapItem;
    }
  }
  
  //Testing________________________________________________________________________________________
  pMapGlobalManip<PMAPTYPE> tester(commDev); 
  assert(tester.check(dofMap));
  assert(tester.check(listMap));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
bool
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
isActive(const DOFCARD & dofCard) const
{
  assert(dofCard.getLocalElId() <= elementIsActive.size());
  return(elementIsActive(dofCard.getLocalElId()));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)  * (ORDER == dms3d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dms3d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dms3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(elementIsActive(dofCard.getLocalElId()));
  
  UInt oldId = 0;
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = grid3d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newVertexLid(oldId) + (dofCard.getLevel() - 1) * numVerticesL);
    
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newEdgeLid(oldId) + (dofCard.getLevel() - 1) * numEdgesL + numVertexDofsL );
    
    case FACE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerFace);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numFaces);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToFace(dofCard.getLocalElId(), dofCard.getLocalId());
      return( newFaceLid(oldId) + (dofCard.getLevel() - 1) * numFacesL + numVertexDofsL + numEdgesDofsL );
    
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= FETYPE::dofPerVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      return( newElementLid(oldId) + (dofCard.getLevel() - 1) * numElementsL + numVertexDofsL + numEdgesDofsL + numFacesDofsL );
      
    default :
      assert(1==2);
   }
  
  return(0); 
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(elementIsActive(dofCard.getLocalElId()));
  
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
UInt
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::DOFMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getDofMap() const
{
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::LISTMAP &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getListMap() const
{
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNumDofsL() const
{
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNumDofsG() const
{
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getSizeListL() const
{
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const UInt &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getSizeListG() const
{
  return(sizeListG);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::FECARDS &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::MESH3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getGrid3d() const
{
  assert(geometryLoaded);
  return(grid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::CONNECT3D> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getConnectGrid3d() const
{
  assert(geometryLoaded);
  return(connectGrid3d);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const Teuchos::RCP<typename dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::OPTIONS> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getVertexIsActive() const
{
  assert(startupOk);
  return(vertexIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getEdgeIsActive() const
{
  assert(startupOk);
  return(edgeIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getFaceIsActive() const
{
  assert(startupOk);
  return(faceIsActive);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<bool> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getElementIsActive() const
{
  assert(startupOk);
  return(elementIsActive);
}
    
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewVertexLid() const
{
  assert(startupOk);
  return(newVertexLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewEdgeLid() const
{
  assert(startupOk);
  return(newEdgeLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewFaceLid() const
{
  assert(startupOk);
  return(newFaceLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const sVect<UInt> &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getNewElementLid() const
{
  assert(startupOk);
  return(newElementLid);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER>
const bool &
dofMapStatic3d<FETYPE,DOFTYPE,ORDER,dms3d_geoIdMode>::
getStartupOk() const
{
  return(startupOk);
}

template<typename F, typename D, dms3d_order O>
ostream & operator<<(ostream & f, const dofMapStatic3d<F,D, O,dms3d_geoIdMode> & V)
{
  f << "numVertices   : " << V.numVerticesL   << " " << V.numVerticesG << endl;
  f << "numEdges      : " << V.numEdgesL      << " " << V.numEdgesG << endl;
  f << "numFaces      : " << V.numFacesL      << " " << V.numFacesG << endl;
  f << "numElements   : " << V.numElementsL   << " " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numFacesDofs  : " << V.numFacesDofsL  << " " << V.numFacesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}



#endif
