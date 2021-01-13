/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef DOFMAPDYNAMIC1D_HPP
#define DOFMAPDYNAMIC1D_HPP

#include "typesInterface.hpp"
#include "traitsBasic.h"

#include "pVectManip.hpp"

#include "morganaGeometry.hpp"
#include "geoMapInterface.hpp"
#include "connect1d.hpp"
#include "mesh1dGlobalManip.hpp"

#include "elCardFeeder1d.hpp"
#include "dofMapDynamic1d_options.h"
#include "feDynamicDofCard1d.h"
#include "morganaGeometry.hpp"


enum dmd1d_order {dmd1d_vectMajor, dmd1d_componentMajor};
enum dmd1d_mode  {dmd1d_standard};


//! Forward declaration
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER = dmd1d_vectMajor, dmd1d_mode MODE = dmd1d_standard> class dofMapDynamic1d;


//_________________________________________________________________________________________________
// STANDARD MODE
//-------------------------------------------------------------------------------------------------

/*! Support for the mapping of 1d dynamic finite elements -> standard mode */
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
class dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>
{
    /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    
    typedef mesh1d<GEOSHAPE,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef connect1d<GEOSHAPE,PMAPTYPE,PMAPTYPE> CONNECT1D;
    typedef geoMapInterface<GEOSHAPE>             GEOMAPINTERFACE;
    typedef pVect<FECARD,PMAPTYPE>                FECARDS;
    typedef dofMapDynamic1d_options               OPTIONS;
    
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
    Teuchos::RCP<MESH1D>        grid1d;
    Teuchos::RCP<CONNECT1D>     connectGrid1d;
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
    sVect<sVect<UInt> > newElementLid;
    
    sVect<GEOVECT>  globVertices;
    sVect<GEOVECT>  globElements;
    
    sVect<sVect<bool> > vertexIsActive;
    sVect<sVect<bool> > elementIsActive;
    
    sVect<UInt> numVerticesL, numVerticesG, offsetVertices;
    sVect<UInt> numElementsL, numElementsG, offsetElements;
    
    UInt maxLevVertex, maxLevVolume;
    
    UInt numVertexDofsL, numVertexDofsG;
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
    dofMapDynamic1d();
    dofMapDynamic1d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    dofMapDynamic1d(communicator & CommDev, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
    dofMapDynamic1d(const dofMapDynamic1d & DofMap);
    dofMapDynamic1d operator=(const dofMapDynamic1d & DofMap);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d);
    void setGeometry(MESH1D & Grid1d, CONNECT1D & ConnedGrid1d);
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
    const Teuchos::RCP<MESH1D>       & getGrid1d() const;
    const Teuchos::RCP<CONNECT1D>    & getConnectGrid1d() const;
    const Teuchos::RCP<OPTIONS>      & getOptions() const;
    const bool                       & getStartupOk() const;
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename F, typename D, dmd1d_order O>
    friend ostream & operator<<(ostream & f, const dofMapDynamic1d<F,D, O, dmd1d_standard> & V);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTOR AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
dofMapDynamic1d()
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
dofMapDynamic1d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = CommDev;
  grid1d        = Grid1d;
  connectGrid1d = ConnedGrid1d;
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
dofMapDynamic1d(communicator & CommDev, MESH1D & Grid1d, CONNECT1D & ConnedGrid1d)
{
  typedef typename FETYPE::GEOSHAPE FE_GEOSHAPE;
  assert(GEOSHAPE::geoName == FE_GEOSHAPE::geoName);
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  commDevLoaded  = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  commDev       = Teuchos::rcpFromRef(CommDev);
  grid1d        = Teuchos::rcp(new MESH1D(Grid1d));
  connectGrid1d = Teuchos::rcp(new CONNECT1D(ConnedGrid1d));
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
dofMapDynamic1d(const dofMapDynamic1d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid1d        = DofMap.grid1d;
  connectGrid1d = DofMap.connectGrid1d;
  options       = DofMap.options;
  
  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
 
  newVertexLid  = DofMap.newVertexLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globElements = DofMap.globElements;
    
  vertexIsActive  = DofMap.vertexIsActive;
  elementIsActive = DofMap.elementIsActive;
    
  numVerticesL   = DofMap.numVerticesL;
  numVerticesG   = DofMap.numVerticesG;
  offsetVertices = DofMap.offsetVertices;
  
  numElementsL  = DofMap.numElementsL;
  numElementsG  = DofMap.numElementsG;
  offsetElements = DofMap.offsetElements;
    
  maxLevVertex = DofMap.maxLevVertex;
  maxLevVolume = DofMap.maxLevVolume;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
  numVolumeDofsL = DofMap.numVolumeDofsL;
  numVolumeDofsG = DofMap.numVolumeDofsG;
    
  numDofsL  = DofMap.numDofsL;
  numDofsG  = DofMap.numDofsG;
  sizeListL = DofMap.sizeListL;
  sizeListG = DofMap.sizeListG;

  dofMap  = DofMap.dofMap;
  listMap = DofMap.listMap;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
operator=(const dofMapDynamic1d & DofMap)
{
  feCards = DofMap.feCards;
  
  commDev       = DofMap.commDev;
  grid1d        = DofMap.grid1d;
  connectGrid1d = DofMap.connectGrid1d;
  options       = DofMap.options;
  
  commDevLoaded  = DofMap.commDevLoaded;
  geometryLoaded = DofMap.geometryLoaded;
  optionsLoaded  = DofMap.optionsLoaded;
  startupOk      = DofMap.startupOk;
 
  newVertexLid  = DofMap.newVertexLid;
  newElementLid = DofMap.newElementLid;
    
  globVertices = DofMap.globVertices;
  globElements = DofMap.globElements;
    
  vertexIsActive  = DofMap.vertexIsActive;
  elementIsActive = DofMap.elementIsActive;
    
  numVerticesL   = DofMap.numVerticesL;
  numVerticesG   = DofMap.numVerticesG;
  offsetVertices = DofMap.offsetVertices;
  
  numElementsL  = DofMap.numElementsL;
  numElementsG  = DofMap.numElementsG;
  offsetElements = DofMap.offsetElements;
    
  maxLevVertex = DofMap.maxLevVertex;
  maxLevVolume = DofMap.maxLevVolume;
    
  numVertexDofsL = DofMap.numVertexDofsL;
  numVertexDofsG = DofMap.numVertexDofsG;
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

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded  = true;
  commDev        = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded  = true;
  commDev        = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setGeometry(const Teuchos::RCP<MESH1D> & Grid1d, const Teuchos::RCP<CONNECT1D> & ConnedGrid1d)
{
  geometryLoaded = true;
  grid1d         = Grid1d;
  connectGrid1d  = ConnedGrid1d;
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setGeometry(MESH1D & Grid1d, CONNECT1D & ConnedGrid1d)
{
  geometryLoaded = true; 
  grid1d         = Teuchos::rcp(new MESH1D(Grid1d));
  connectGrid1d  = Teuchos::rcp(new CONNECT1D(ConnedGrid1d));
  
  feCards.resize(grid1d->getNumElements());
  feCards.setMap(grid1d->getElements().getRowMap());
  feCards.updateFinder();
}



//_________________________________________________________________________________________________
// OPTIONS LOADING AND STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(lid <= feCards.sizeL());
  feCards.getL(lid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  assert(feCards.isG(gid));
  feCards.getG(gid) = FeCards;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  assert(feCards.size() == FECards.size());
  feCards = FECards;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  options = Options;
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  options = Teuchos::rcp(new OPTIONS(Options));
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
void
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
startup()
{ 
  //Asserts
  assert(commDevLoaded);
  assert(geometryLoaded);
  assert(optionsLoaded);
  
  assert(grid1d->getElements().colIsLocal());
  
  //Checking
  mesh1dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> manipulator(commDev);
  assert(manipulator.check(grid1d));
  
  pVectGlobalManip<FECARD,PMAPTYPE> feCardsCheker(commDev);
  assert(feCardsCheker.check(feCards));
  
  //Flag update
  startupOk = true;
  
  //Allocate structures
  UInt numProc = commDev->size();
  elCardFeeder1d<GEOSHAPE,PMAPTYPE> feeder(grid1d,connectGrid1d);
  
  
  //Maximum levels---------------------------------------------------------------------------------
  maxLevVertex = 0;
  maxLevVolume = 0;
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    refFe.setCards(feCards(i), feeder.getCardLocal(i));
    
    maxLevVertex = std::max(maxLevVertex, refFe.getMaxLevVertex());
    maxLevVolume = std::max(maxLevVolume, refFe.getMaxLevVolume());
  }
  
  sVect<UInt> levsVertex(numProc);
  sVect<UInt> levsVolume(numProc);
  
  
  all_gather(*commDev, maxLevVertex, levsVertex);
  all_gather(*commDev, maxLevVolume, levsVolume);
  
  for(UInt i=1; i <= numProc; ++i)
  {
    maxLevVertex = std::max(maxLevVertex, levsVertex(i));
    maxLevVolume = std::max(maxLevVolume, levsVolume(i));
  }
  
  //GeoItems are active----------------------------------------------------------------------------
  vertexIsActive.resize(maxLevVertex);
  elementIsActive.resize(maxLevVolume);
  
  for(UInt i=1; i<=maxLevVertex; ++i)
  { vertexIsActive(i).resize( grid1d->getNumVertices() ); }
  
  for(UInt i=1; i<=maxLevVolume; ++i)
  { elementIsActive(i).resize( grid1d->getNumElements() ); }
  
  //Default to false
  for(UInt i=1; i <= maxLevVertex; ++i)
  {
    for(UInt j=1; j <= vertexIsActive(i).size(); ++j)
    { vertexIsActive(i)(j) = false; }
  }
  
  for(UInt i=1; i <= maxLevVolume; ++i)
  {
    for(UInt j=1; j <= elementIsActive(i).size(); ++j)
    { elementIsActive(i)(j) = false; }
  }
  
  //Identify the active items
  UInt lev, locId;
  sVect<DOFCARD> dofCards;
  
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
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
	  vertexIsActive(lev)(grid1d->getElementL(i).getCid(locId)) = true;
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
  { newVertexLid(i).resize( grid1d->getNumVertices() ); }
  
  newElementLid.resize(maxLevVolume);
  for(UInt i=1; i<=maxLevVolume; ++i)
  { newElementLid(i).resize( grid1d->getNumElements() ); }
  
  //Reset the numbers
  numVerticesL.resize(maxLevVertex);
  numVerticesG.resize(maxLevVertex);
  
  for(UInt i=1; i<=maxLevVertex; ++i)
  {
    numVerticesL(i) = 0;
    numVerticesG(i) = 0;
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
    for(UInt i=1; i <= grid1d->getNumVertices(); ++i)
    {
      if(vertexIsActive(lev)(i))
      {
	numVerticesL(lev)++;
	newVertexLid(lev)(i) = numVerticesL(lev);
	
	mapItem = grid1d->getNodes().getMapL(i);
        mapItem.setLid(numVerticesL(lev));
	
	globVertices(lev).push_back(mapItem, mapItem.getGid());
      }
    }
  }
  
  //Elements new lids
  globElements.resize(maxLevVolume);
  
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  { globElements(lev).clear(); }
  
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  {
    for(UInt i=1; i <= grid1d->getNumElements(); ++i)  //Elements
    {
      if(elementIsActive(lev)(i))
      {
	numElementsL(lev)++;
	newElementLid(lev)(i) = numElementsL(lev);
	
	mapItem = grid1d->getElements().getRowMapL(i);
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
    
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  { 
    geoItemsManip.buildGlobalNumbering(globElements(lev));
    numElementsG(lev) = geoItemsManip.sizeG(globElements(lev));
  }
  
  
  //Count the dofs---------------------------------------------------------------------------------
  numVertexDofsL = 0;
  numVolumeDofsL = 0;
  
  numVertexDofsG = 0;
  numVolumeDofsG = 0;
  
  for(UInt lev=1; lev <= maxLevVertex; ++lev)
  {
    numVertexDofsL += numVerticesL(lev);
    numVertexDofsG += numVerticesG(lev);
  }
  
  for(UInt lev=1; lev <= maxLevVolume; ++lev)
  {
    numVolumeDofsL += numElementsL(lev);
    numVolumeDofsG += numElementsG(lev);
  }
  
  
  //Total number of dofs and lists-----------------------------------------------------------------
  numDofsL = numVertexDofsL + numVolumeDofsL;
  numDofsG = numVertexDofsG + numVolumeDofsG;  
  
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
      
      lid = (block * (mapItem.getLid() - 1) + K) * (ORDER == dmd1d_vectMajor)  + (mapItem.getLid() + (K-1) * numDofsL) * (ORDER == dmd1d_componentMajor);
      gid = (block * (mapItem.getGid() - 1) + K) * (ORDER == dmd1d_vectMajor)  + (mapItem.getGid() + (K-1) * numDofsG) * (ORDER == dmd1d_componentMajor);
      
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
  offsetElements.resize(maxLevVolume);
  
  if(maxLevVertex != 0) {offsetVertices(1) = 0;}
  if(maxLevVolume != 0) {offsetElements(1) = 0;}
  
  for(UInt lev=2; lev <= maxLevVertex; ++lev)
  { offsetVertices(lev) = offsetVertices(lev-1) + numVerticesL(lev-1); }
  
  for(UInt lev=2; lev <= maxLevVolume; ++lev)
  { offsetElements(lev) = offsetElements(lev-1) + numElementsL(lev-1); }
}


//_________________________________________________________________________________________________
// GET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
bool
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
isActive(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  UInt oldId;
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());
      
      oldId = grid1d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      return(vertexIsActive(dofCard.getLevel())(oldId));    
      
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      return(elementIsActive(dofCard.getLevel())(oldId));
      
    default :
      assert(1==2);
  }
  
  return(false);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
UInt
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
dofToListL(const UInt & dofL, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofL - 1) + K)   * (ORDER == dmd1d_vectMajor) +
          (dofL  + (K-1) * numDofsL) * (ORDER == dmd1d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
UInt
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
dofToListG(const UInt & dofG, const UInt & I, const UInt & J) const
{
  UInt K  = I + (J-1) * traitsBasic<DOFTYPE>::numI;
  
  return( (block * (dofG - 1) + K)   * (ORDER == dmd1d_vectMajor) +
          (dofG  + (K-1) * numDofsG) * (ORDER == dmd1d_componentMajor) );
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
UInt
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
mapDofL(const DOFCARD & dofCard) const
{
  assert(startupOk);
  
  UInt oldId;
  
  switch(dofCard.getGeoType())
  {
    case VERTEX :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVertex);
      assert(dofCard.getLocalId()   >= 1);  assert(dofCard.getLocalId()   <= GEOSHAPE::numVertices);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());
      
      oldId = grid1d->getElements().getCid_LL(dofCard.getLocalElId(), dofCard.getLocalId());
      assert(vertexIsActive(dofCard.getLevel())(oldId));
      
      return( newVertexLid(dofCard.getLevel())(oldId) + offsetVertices(dofCard.getLevel()));      
      
    case VOLUME :
      assert(dofCard.getLevel()     >= 1);  assert(dofCard.getLevel()     <= maxLevVolume);
      assert(dofCard.getLocalId()   == 1);
      assert(dofCard.getLocalElId() >= 1);  assert(dofCard.getLocalElId() <= grid1d->getNumElements());
      
      oldId = dofCard.getLocalElId();
      assert(elementIsActive(dofCard.getLevel())(oldId));
      
      return( newElementLid(dofCard.getLevel())(oldId) + offsetElements(dofCard.getLevel()) + numVertexDofsL );
      
    default :
      assert(1==2);
  }
  
  return(0);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
UInt
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
mapDofG(const DOFCARD & dofCard) const
{
  assert(startupOk);
  return(dofMap(mapDofL(dofCard)).getGid());
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
UInt
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
mapListL(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofL = mapDofL(dofCard);
  return(dofToListL(dofL,I,J));
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
UInt
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
mapListG(const UInt & I, const UInt & J, const DOFCARD & dofCard) const
{
  assert(startupOk);
  assert(I >= 1); assert(I <= traitsBasic<DOFTYPE>::numI);
  assert(J >= 1); assert(J <= traitsBasic<DOFTYPE>::numJ);
  
  UInt dofG = mapDofG(dofCard);
  return(dofToListG(dofG,I,J));
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const UInt &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getNumDofsL() const
{
  assert(startupOk);
  return(numDofsL);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const UInt &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getNumDofsG() const
{
  assert(startupOk);
  return(numDofsG);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const UInt &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getSizeListL() const
{
  assert(startupOk);
  return(sizeListL);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const UInt &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getSizeListG() const
{
  assert(startupOk);
  return(sizeListG);
}



//_________________________________________________________________________________________________
// DUMP FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const typename dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::DOFMAP  &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getDofMap() const
{
  assert(startupOk);
  return(dofMap);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const typename dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::LISTMAP &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getListMap() const
{
  assert(startupOk);
  return(listMap);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const typename dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::FECARDS &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getFeCards() const
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
typename dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::FECARDS &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getFeCards()
{
  return(feCards);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const Teuchos::RCP<communicator> &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getCommDev() const
{
  assert(commDevLoaded);
  return(commDev);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const Teuchos::RCP<typename dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::MESH1D> &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getGrid1d() const
{
  assert(geometryLoaded);
  return(grid1d);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const Teuchos::RCP<typename dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::CONNECT1D> &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getConnectGrid1d() const
{
  assert(geometryLoaded);
  return(connectGrid1d);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const Teuchos::RCP<typename dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::OPTIONS> &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getOptions() const
{
  assert(optionsLoaded);
  return(options);
}

template<typename FETYPE, typename DOFTYPE, dmd1d_order ORDER>
const bool &
dofMapDynamic1d<FETYPE,DOFTYPE,ORDER,dmd1d_standard>::
getStartupOk() const
{
  return(startupOk);
}



//_________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------
template<typename F, typename D, dmd1d_order O>
ostream & operator<<(ostream & f, const dofMapDynamic1d<F,D, O, dmd1d_standard> & V)
{
  f << "numVerticesL   : " << V.numVerticesL   << endl;
  f << "numElementsL   : " << V.numElementsL   << endl << endl;
  
  f << "numVerticesG   : " << V.numVerticesG << endl;
  f << "numElementsG   : " << V.numElementsG << endl << endl;
  
  f << "numVertexDofs : " << V.numVertexDofsL << " " << V.numVertexDofsG << endl;
  f << "numVolumeDofs : " << V.numVolumeDofsL << " " << V.numVolumeDofsG << endl << endl;
  
  f << "numTotalDofs  : " << V.numDofsL       << " " << V.numDofsG << endl;
  f << "sizeList      : " << V.sizeListL      << " " << V.sizeListG << endl;
  
  return(f);
}


#endif
