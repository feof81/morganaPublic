/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef DOFMAPDYNAMIC3D_HPP
#define DOFMAPDYNAMIC3D_HPP

#include "typesInterface.hpp"
#include "traitsBasic.h"

#include "pVectManip.hpp"

#include "morganaGeometry.hpp"
#include "geoMapInterface.hpp"
#include "connect3d.hpp"
#include "mesh3dGlobalManip.hpp"

#include "elCardFeeder3d.hpp"
#include "dofMapDynamic3d_options.h"
#include "feDynamicDofCard3d.h"
#include "morganaGeometry.hpp"


enum dmd3d_order {dmd3d_vectMajor, dmd3d_componentMajor};
enum dmd3d_mode  {dmd3d_standard};


//! Forward declaration
template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER = dmd3d_vectMajor, dmd3d_mode MODE = dmd3d_standard> class dofMapDynamic3d;


//_________________________________________________________________________________________________
// STANDARD MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 3d dynamic finite elements -> standard mode */
template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
class dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>
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
    typedef dofMapDynamic3d_options               OPTIONS;
    
    typedef pMap<PMAPTYPE>       DOFMAP;
    typedef pMap<PMAPTYPE>       LISTMAP;
    typedef pVect<UInt,PMAPTYPE> GEOVECT;
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
    sVect<sVect<UInt> > newVertexLid;
    sVect<sVect<UInt> > newEdgeLid;
    sVect<sVect<UInt> > newFaceLid;
    sVect<sVect<UInt> > newElementLid;
    
    sVect<GEOVECT>  globVertices;
    sVect<GEOVECT>  globEdges;
    sVect<GEOVECT>  globFaces;
    sVect<GEOVECT>  globElements;
    
    sVect<sVect<bool> > vertexIsActive;
    sVect<sVect<bool> > edgeIsActive;
    sVect<sVect<bool> > faceIsActive;
    sVect<sVect<bool> > elementIsActive;
    
    sVect<UInt> numVerticesL, numVerticesG, offsetVertices;
    sVect<UInt> numEdgesL,    numEdgesG,    offsetEdges;
    sVect<UInt> numFacesL,    numFacesG,    offsetFaces;
    sVect<UInt> numElementsL, numElementsG, offsetElements;
    
    UInt maxLevVertex, maxLevEdge, maxLevFace, maxLevVolume;
    
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
    dofMapDynamic3d();
    dofMapDynamic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d);
    dofMapDynamic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d);
    dofMapDynamic3d(const dofMapDynamic3d & DofMap);
    dofMapDynamic3d operator=(const dofMapDynamic3d & DofMap);
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
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dmd3d_order O>
    friend ostream & operator<<(ostream & f, const dofMapDynamic3d<F,D, O,dmd3d_standard> & V);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTOR AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
dofMapDynamic3d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
dofMapDynamic3d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
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

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
dofMapDynamic3d(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnedGrid3d)
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

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
dofMapDynamic3d(const dofMapDynamic3d & DofMap)
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
 
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid;
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
    
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
    
  numVerticesL   = DofMap.numVerticesL;
  numVerticesG   = DofMap.numVerticesG;
  offsetVertices = DofMap.offsetVertices;
  
  numEdgesL   = DofMap.numEdgesL;
  numEdgesG   = DofMap.numEdgesG;
  offsetEdges = DofMap.offsetEdges;
  
  numFacesL   = DofMap.numFacesL;
  numFacesG   = DofMap.numFacesG;
  offsetFaces = DofMap.offsetFaces;
  
  numElementsL  = DofMap.numElementsL;
  numElementsG  = DofMap.numElementsG;
  offsetElements = DofMap.offsetElements;
    
  maxLevVertex = DofMap.maxLevVertex;
  maxLevEdge   = DofMap.maxLevEdge;
  maxLevFace   = DofMap.maxLevFace;
  maxLevVolume = DofMap.maxLevVolume;
    
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

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
operator=(const dofMapDynamic3d & DofMap)
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
 
  newVertexLid  = DofMap.newVertexLid;
  newEdgeLid    = DofMap.newEdgeLid;
  newFaceLid    = DofMap.newFaceLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globEdges    = DofMap.globEdges;
  globFaces    = DofMap.globFaces;
  globElements = DofMap.globElements;
    
  vertexIsActive  = DofMap.vertexIsActive;
  edgeIsActive    = DofMap.edgeIsActive;
  faceIsActive    = DofMap.faceIsActive;
  elementIsActive = DofMap.elementIsActive;
    
  numVerticesL   = DofMap.numVerticesL;
  numVerticesG   = DofMap.numVerticesG;
  offsetVertices = DofMap.offsetVertices;
  
  numEdgesL   = DofMap.numEdgesL;
  numEdgesG   = DofMap.numEdgesG;
  offsetEdges = DofMap.offsetEdges;
  
  numFacesL   = DofMap.numFacesL;
  numFacesG   = DofMap.numFacesG;
  offsetFaces = DofMap.offsetFaces;
  
  numElementsL  = DofMap.numElementsL;
  numElementsG  = DofMap.numElementsG;
  offsetElements = DofMap.offsetElements;
    
  maxLevVertex = DofMap.maxLevVertex;
  maxLevEdge   = DofMap.maxLevEdge;
  maxLevFace   = DofMap.maxLevFace;
  maxLevVolume = DofMap.maxLevVolume;
    
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

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev        = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev        = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setGeometry(const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnedGrid3d)
{
  geometryLoaded = true;
  grid3d         = Grid3d;
  connectGrid3d  = ConnedGrid3d;
  
  feCards.resize(grid3d->getNumElements());
  feCards.setMap(grid3d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
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
template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
void
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
startup()
{ 
  //Asserts
  assert(commDevLoaded);
  assert(geometryLoaded);
  assert(optionsLoaded);
  
  assert(grid3d->getElements().colIsLocal());
  assert(grid3d->getFaces().colIsLocal());
  assert(grid3d->getEdges().colIsLocal());
  
  //Checking
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid3d));
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsCheker(commDev);
  assert(feCardsCheker.check(feCards));
  
  //Flag update
  startupOk = true;
  
  //Allocate structures
  UInt numProc = commDev->size();
  elCardFeeder3d<GEOSHAPE,PMAPTYPE> feeder(grid3d,connectGrid3d);
  
  
  //Maximum levels---------------------------------------------------------------------------------
  maxLevVertex = 0;
  maxLevEdge   = 0;
  maxLevFace   = 0;
  maxLevVolume = 0;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    refFe.setCards(feCards(i), feeder.getCardLocal(i));
    
    maxLevVertex = std::max(maxLevVertex, refFe.getMaxLevVertex());
    maxLevEdge   = std::max(maxLevEdge,   refFe.getMaxLevEdge());
    maxLevFace   = std::max(maxLevFace,   refFe.getMaxLevFace());
    maxLevVolume = std::max(maxLevVolume, refFe.getMaxLevVolume());
  }
  
  sVect<UInt> levsVertex(numProc);
  sVect<UInt> levsEdge(numProc);
  sVect<UInt> levsFace(numProc);
  sVect<UInt> levsVolume(numProc);
  
  
  all_gather(*commDev, maxLevVertex, levsVertex);
  all_gather(*commDev, maxLevEdge,   levsEdge);
  all_gather(*commDev, maxLevFace,   levsFace);
  all_gather(*commDev, maxLevVolume, levsVolume);
  
  for(UInt i=1; i <= numProc; ++i)
  {
    maxLevVertex = std::max(maxLevVertex, levsVertex(i));
    maxLevEdge   = std::max(maxLevEdge,   levsEdge(i));
    maxLevFace   = std::max(maxLevFace,   levsFace(i));
    maxLevVolume = std::max(maxLevVolume, levsVolume(i));
  }
  
  //GeoItems are active----------------------------------------------------------------------------
  vertexIsActive.resize(maxLevVertex);
  edgeIsActive.resize(maxLevEdge);
  faceIsActive.resize(maxLevFace);
  elementIsActive.resize(maxLevVolume);
  
  for(UInt i=1; i<=maxLevVertex; ++i)
  { vertexIsActive(i).resize( grid3d->getNumVertices() ); }
  
  for(UInt i=1; i<=maxLevEdge; ++i)
  { edgeIsActive(i).resize( grid3d->getNumEdges() ); }
  
  for(UInt i=1; i<=maxLevFace; ++i)
  { faceIsActive(i).resize( grid3d->getNumFaces() ); }
  
  for(UInt i=1; i<=maxLevVolume; ++i)
  { elementIsActive(i).resize( grid3d->getNumElements() ); }
  
  //Default to false
  for(UInt i=1; i <= maxLevVertex; ++i)
  {
    for(UInt j=1; j <= vertexIsActive(i).size(); ++j)
    { vertexIsActive(i)(j) = false; }
  }
  
  for(UInt i=1; i <= maxLevEdge; ++i)
  {
    for(UInt j=1; j<=edgeIsActive(i).size(); ++j)
    { edgeIsActive(i)(j) = false; }
  }
  
  for(UInt i=1; i <= maxLevFace; ++i)
  {
    for(UInt j=1; j <= faceIsActive(i).size(); ++j)
    { faceIsActive(i)(j) = false; }
  }
  
  for(UInt i=1; i <= maxLevVolume; ++i)
  {
    for(UInt j=1; j <= elementIsActive(i).size(); ++j)
    { elementIsActive(i)(j) = false; }
  }
  
  //Identify the active items
  UInt lev, locId;
  sVect<DOFCARD> dofCards;
  
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    refFe.setCards(feCards(i), feeder.getCardLocal(i));
    dofCards = refFe.getDofCards();
    assert(dofCards.size() == refFe.getNumBasis());
        
    for(UInt j=1; j<= refFe.getNumBasis(); ++j)
    {
      lev   = dofCards(j).getLevel();
      locId = dofCards(j).getLocalId();
      
      switch(dofCards(j).getGeoType())
      {
	case VERTEX :
	  vertexIsActive(lev)(grid3d->getElementL(i).getCid(locId)) = true;
	break;
	
	case EDGE :
	  edgeIsActive(lev)(connectGrid3d->getElementToEdge(i,locId)) = true;
	  break;
	  
	case FACE :
	  faceIsActive(lev)(connectGrid3d->getElementToFace(i,locId)) = true;
	  break;
	  
	case VOLUME :
	  elementIsActive(lev)(i) = true;
	  break;
	  
	default : 
	  assert(1==2);
      }
    }
  }
  

  //Assign new lids--------------------------------------------------------------------------------
  PMAPTYPE mapItem;
  
  newVertexLid.resize(maxLevVertex);
  for(UInt i=1; i<=maxLevVertex; ++i)
  { newVertexLid(i).resize( grid3d->getNumVertices() ); }
  
  newEdgeLid.resize(maxLevEdge);
  for(UInt i=1; i<=maxLevEdge; ++i)
  { newEdgeLid(i).resize( grid3d->getNumEdges() ); }
  
  newFaceLid.resize(maxLevFace);
  for(UInt i=1; i<=maxLevFace; ++i)
  { newFaceLid(i).resize( grid3d->getNumFaces() ); }
  
  newElementLid.resize(maxLevVolume);
  for(UInt i=1; i<=maxLevVolume; ++i)
  { newElementLid(i).resize( grid3d->getNumElements() ); }
  
  //Reset the numbers
  numVerticesL.resize(maxLevVertex);
  numVerticesG.resize(maxLevVertex);
  
  for(UInt i=1; i<=maxLevVertex; ++i)
  {
    numVerticesL(i) = 0;
    numVerticesG(i) = 0;
  }
  
  numEdgesL.resize(maxLevEdge);
  numEdgesG.resize(maxLevEdge);
  
  for(UInt i=1; i<=maxLevEdge; ++i)
  {
    numEdgesL(i) = 0;
    numEdgesG(i) = 0;
  }
  
  numFacesL.resize(maxLevFace);
  numFacesG.resize(maxLevFace);
  
  for(UInt i=1; i<=maxLevFace; ++i)
  {
    numFacesL(i) = 0;
    numFacesG(i) = 0;
  }
  
  numElementsL.resize(maxLevVolume);
  numElementsG.resize(maxLevVolume);
  
  for(UInt i=1; i<=maxLevVolume; ++i)
  {
    numElementsL(i) = 0;
    numElementsG(i) = 0;
  }
  
  //Vertex new lids
  globVertices.resize(maxLevVertex);
  
  for(UInt lev=1; lev <= maxLevVertex; ++lev)
  { globVertices(lev).clear(); }
  
  for(UInt lev=1; lev <= maxLevVertex; ++lev)
  {
    for(UInt i=1; i <= grid3d->getNumVertices(); ++i)
    {
      if(vertexIsActive(lev)(i))
      {
	numVerticesL(lev)++;
	newVertexLid(lev)(i) = numVerticesL(lev);
	
	mapItem = grid3d->getNodes().getMapL(i);
        mapItem.setLid(numVerticesL(lev));
	
	globVertices(lev).push_back(mapItem, mapItem.getGid());
      }
    }
  }
  
  //Edges new lids
  globEdges.resize(maxLevEdge);
  
  for(UInt lev=1; lev <= maxLevEdge; ++lev)
  { globEdges(lev).clear(); }
  
  for(UInt lev=1; lev <= maxLevEdge; ++lev)
  {
    for(UInt i=1; i <= grid3d->getNumEdges(); ++i)
    {
      if(edgeIsActive(lev)(i))
      {
	numEdgesL(lev)++;
	newEdgeLid(lev)(i) = numEdgesL(lev);
	
	mapItem = grid3d->getEdges().getMapL(i);
        mapItem.setLid(numEdgesL(lev));
	
	globEdges(lev).push_back(mapItem, mapItem.getGid());
      }
    }
  }
  
  //Faces new lids
  globFaces.resize(maxLevFace);
  
  for(UInt lev=1; lev <= maxLevFace; ++lev)
  { globFaces(lev).clear(); }
  
  for(UInt lev=1; lev <= maxLevFace; ++lev)
  {
    for(UInt i=1; i <= grid3d->getNumFaces(); ++i)
    {
      if(faceIsActive(lev)(i))
      {
	numFacesL(lev)++;
	newFaceLid(lev)(i) = numFacesL(lev);
	
	mapItem = grid3d->getFaces().getMapL(i);
        mapItem.setLid(numFacesL(lev));
	
	globFaces(lev).push_back(mapItem, mapItem.getGid());
      }
    }
  }
  
  //Elements new lids
  globElements.resize(maxLevVolume);
  
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  { globElements(lev).clear(); }
  
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  {
    for(UInt i=1; i <= grid3d->getNumElements(); ++i)  //Elements
    {
      if(elementIsActive(lev)(i))
      {
	numElementsL(lev)++;
	newElementLid(lev)(i) = numElementsL(lev);
	
	mapItem = grid3d->getElements().getRowMapL(i);
	mapItem.setLid(numElementsL(lev));
	
	globElements(lev).push_back(mapItem, mapItem.getGid());
      }
    }
  }
  
  
  //Global mapping---------------------------------------------------------------------------------
  pVectGlobalManip<UInt,PMAPTYPE> geoItemsManip(commDev);
  
  for(UInt lev=1; lev <= maxLevVertex; ++lev)
  { 
    geoItemsManip.buildGlobalNumbering(globVertices(lev));
    numVerticesG(lev) = geoItemsManip.sizeG(globVertices(lev));
  }  
  
  for(UInt lev=1; lev <= maxLevEdge; ++lev)
  {
    geoItemsManip.buildGlobalNumbering(globEdges(lev));
    numEdgesG(lev) = geoItemsManip.sizeG(globEdges(lev));
  }
    
  for(UInt lev=1; lev <= maxLevFace; ++lev)
  {
    geoItemsManip.buildGlobalNumbering(globFaces(lev));
    numFacesG(lev) = geoItemsManip.sizeG(globFaces(lev));
  }  
    
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  { 
    geoItemsManip.buildGlobalNumbering(globElements(lev));
    numElementsG(lev) = geoItemsManip.sizeG(globElements(lev));
  }
  
  
  //Count the dofs---------------------------------------------------------------------------------
  numVertexDofsL = 0;
  numEdgesDofsL  = 0;
  numFacesDofsL  = 0;
  numVolumeDofsL = 0;
  
  numVertexDofsG = 0;
  numEdgesDofsG  = 0;
  numFacesDofsG  = 0;
  numVolumeDofsG = 0;
  
  for(UInt lev=1; lev <= maxLevVertex; ++lev)
  {
    numVertexDofsL += numVerticesL(lev);
    numVertexDofsG += numVerticesG(lev);
  }
  
  for(UInt lev=1; lev <= maxLevEdge; ++lev)
  {
    numEdgesDofsL += numEdgesL(lev);
    numEdgesDofsG += numEdgesG(lev);
  }
  
  for(UInt lev=1; lev <= maxLevFace; ++lev)
  {
    numFacesDofsL += numFacesL(lev);
    numFacesDofsG += numFacesG(lev);
  }
  
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  {
    numVolumeDofsL += numElementsL(lev);
    numVolumeDofsG += numElementsG(lev);
  }
  
  
  //Total number of dofs and lists-----------------------------------------------------------------
  numDofsL = numVertexDofsL + numEdgesDofsL + numFacesDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numEdgesDofsG + numFacesDofsG + numVolumeDofsG;  
  
  sizeListL = numDofsL * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  sizeListG = numDofsG * traitsBasic<DOFTYPE>::numI * traitsBasic<DOFTYPE>::numJ;
  
  //Build dof map__________________________________________________________________________________
  UInt offsetG = 0, k=1;
  dofMap.resize(numDofsL);
  
  //Vertices
  for(UInt lev=1; lev <= maxLevVertex; ++lev)
  {
    for(UInt i=1; i <= numVerticesL(lev); ++i)
    {
      mapItem = globVertices(lev).getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + offsetG );
      
      dofMap(k) = mapItem;
      ++k;
    }
    
    offsetG += numVerticesG(lev);
  }
  
  //Edges
  for(UInt lev=1; lev <= maxLevEdge; ++lev)  //Edges
  {
    for(UInt i=1; i <= numEdgesL(lev); ++i)
    {
      mapItem = globEdges(lev).getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + offsetG );
      
      dofMap(k) = mapItem;
      ++k;
    }
    
    offsetG += numEdgesG(lev);
  }
  
  //Faces
  for(UInt lev=1; lev <= maxLevFace; ++lev)  //Faces
  {
    for(UInt i=1; i <= numFacesL(lev); ++i)
    {
      mapItem = globFaces(lev).getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + offsetG );
      
      dofMap(k) = mapItem;
      ++k;
    }
    
    offsetG += numFacesG(lev);
  }
  
  //Elements
  for(UInt lev=1; lev <= maxLevVolume; ++lev)  //Volume
  {
    for(UInt i=1; i <= numElementsL(lev); ++i)
    {
      mapItem = globElements(lev).getMapL(i);
      mapItem.setLid( k );
      mapItem.setGid( mapItem.getGid() + offsetG );
      
      dofMap(k) = mapItem;
      ++k;
    }
    
    offsetG += numElementsG(lev);
  }
  
  
  //Build list map---------------------------------------------------------------------------------
  listMap.resize(sizeListL);
  
  UInt lid, gid;
  
  for(UInt i=1; i <= dofMap.size(); ++i)
  {   
    for(UInt K=1; K <= block; ++K)
    {
      mapItem = dofMap(i);
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dmd3d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dmd3d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dmd3d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dmd3d_componentMajor);
      
      mapItem.setLid(lid);
      mapItem.setGid(gid);
      
      listMap(lid) = mapItem;
    }
  }
  
  
  //Testing----------------------------------------------------------------------------------------
  pMapGlobalManip<PMAPTYPE> tester(commDev); 
  assert(tester.check(dofMap));
  assert(tester.check(listMap));
  
  
  //Build offsets----------------------------------------------------------------------------------
  offsetVertices.resize(maxLevVertex);
  offsetEdges.resize(maxLevEdge);
  offsetFaces.resize(maxLevFace);
  offsetElements.resize(maxLevVolume);
  
  if(maxLevVertex != 0) {offsetVertices(1) = 0;}
  if(maxLevEdge   != 0) {offsetEdges(1)    = 0;}
  if(maxLevFace   != 0) {offsetFaces(1)    = 0;}
  if(maxLevVolume != 0) {offsetElements(1) = 0;}
  
  for(UInt lev=2; lev <= maxLevVertex; ++lev)
  { offsetVertices(lev) = offsetVertices(lev-1) + numVerticesL(lev-1); }
  
  for(UInt lev=2; lev <= maxLevEdge; ++lev)
  { offsetEdges(lev) = offsetEdges(lev-1) + numEdgesL(lev-1); }
  
  for(UInt lev=2; lev <= maxLevFace; ++lev)
  { offsetFaces(lev) = offsetFaces(lev-1) + numFacesL(lev-1); }
  
  for(UInt lev=2; lev <= maxLevVolume; ++lev)
  { offsetElements(lev) = offsetElements(lev-1) + numElementsL(lev-1); }
}


//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
bool
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
isActive(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  UInt oldId;
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = grid3d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      return(vertexIsActive(dofCard.getLevel())(oldId));
      
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId());
      return(edgeIsActive(dofCard.getLevel())(oldId));
      
    case FACE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevFace);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numFaces);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToFace(dofCard.getLocalElId(), dofCard.getLocalId());
      return(faceIsActive(dofCard.getLevel())(oldId));
      
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      return(elementIsActive(dofCard.getLevel())(oldId));
      
    default :
      assert(1==2);
  }
  
  return(false);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
UInt
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)   * (ORDER == dmd3d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dmd3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
UInt
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dmd3d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dmd3d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
UInt
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  UInt oldId;
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = grid3d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      assert(vertexIsActive(dofCard.getLevel())(oldId));
      
      return( newVertexLid(dofCard.getLevel())(oldId) + offsetVertices(dofCard.getLevel()));
      
    case EDGE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevEdge);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numEdges);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToEdge(dofCard.getLocalElId(), dofCard.getLocalId());
      assert(edgeIsActive(dofCard.getLevel())(oldId));
      
      return( newEdgeLid(dofCard.getLevel())(oldId) + offsetEdges(dofCard.getLevel()) + numVertexDofsL );
      
    case FACE   :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevFace);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numFaces);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = connectGrid3d->getElementToFace(dofCard.getLocalElId(), dofCard.getLocalId());
      assert(faceIsActive(dofCard.getLevel())(oldId));
      
      return( newFaceLid(dofCard.getLevel())(oldId) + offsetFaces(dofCard.getLevel()) + numVertexDofsL + numEdgesDofsL );
      
      
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid3d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      assert(elementIsActive(dofCard.getLevel())(oldId));
      
      return( newElementLid(dofCard.getLevel())(oldId) + offsetElements(dofCard.getLevel()) + numVertexDofsL + numEdgesDofsL + numFacesDofsL );
      
    default :
      assert(1==2);
  }
  
  return(0);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
UInt
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
UInt
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
UInt
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const UInt &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getNumDofsL() const
{
  assert(startupOk);
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const UInt &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getNumDofsG() const
{
  assert(startupOk);
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const UInt &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getSizeListL() const
{
  assert(startupOk);
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const UInt &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getSizeListG() const
{
  assert(startupOk);
  return(sizeListG);
}



//_________________________________________________________________________________________________
// DUMP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const typename dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::DOFMAP  &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getDofMap() const
{
  assert(startupOk);
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const typename dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::LISTMAP &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getListMap() const
{
  assert(startupOk);
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const typename dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::FECARDS &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
typename dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::FECARDS &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const Teuchos::RCP<typename dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::MESH3D> &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getGrid3d() const
{
  assert(geometryLoaded);
  return(grid3d);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const Teuchos::RCP<typename dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::CONNECT3D> &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getConnectGrid3d() const
{
  assert(geometryLoaded);
  return(connectGrid3d);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const Teuchos::RCP<typename dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::OPTIONS> &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dmd3d_order ORDER>
const bool &
dofMapDynamic3d<FETYPE,DOFTYPE,ORDER,dmd3d_standard>::
getStartupOk() const
{
  return(startupOk);
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename F, typename D, dmd3d_order O>
ostream & operator<<(ostream & f, const dofMapDynamic3d<F,D, O, dmd3d_standard> & V)
{
  f << "numVerticesL   : " << V.numVerticesL   << endl;
  f << "numEdgesL      : " << V.numEdgesL      << endl;
  f << "numFacesL      : " << V.numFacesL      << endl;
  f << "numElementsL   : " << V.numElementsL   << endl << endl;
  
  f << "numVerticesG   : " << V.numVerticesG << endl;
  f << "numEdgesG      : " << V.numEdgesG << endl;
  f << "numFacesG      : " << V.numFacesG << endl;
  f << "numElementsG   : " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numEdgesDofs  : " << V.numEdgesDofsL  << " " << V.numEdgesDofsG << endl;
  f << "numFacesDofs  : " << V.numFacesDofsL  << " " << V.numFacesDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}


#endif
