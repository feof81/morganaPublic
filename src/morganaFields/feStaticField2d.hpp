/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FESTATICFIELD2D_HPP
#define FESTATICFIELD2D_HPP

#include "dofMapStatic2d.hpp"
#include "feStaticInterface2d.hpp"
#include "elCardFeeder2d.hpp"


/*! Field 2d for the static finite elements */
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER = dms2d_vectMajor, dms2d_mode MODE = dms2d_allMode>
class feStaticField2d
{
     /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef typename FETYPE::GEOSHAPE          GEOSHAPE;
    typedef typename FETYPE::DOFCARD           DOFCARD;
    typedef typename FETYPE::FECARD            FECARD;
    typedef typename FETYPE::ELCARD            ELCARD;
    typedef typename FETYPE::BASETYPE          BASETYPE;
    
    typedef dofMapStatic2d<FETYPE,DOFTYPE,ORDER,MODE>       DOFMAPPER;
    typedef feStaticInterface2d<FETYPE,DOFTYPE,ORDER,MODE>  FEINTERFACE;
    typedef elCardFeeder2d<GEOSHAPE,PMAPTYPE>               ELCARDFEEDER;
    typedef pVect<DOFTYPE,PMAPTYPE>                         DOFVECT;
    typedef pVect<Real,PMAPTYPE>                            LISTVECT;
    typedef pMap<PMAPTYPE>                                  LISTMAP;
    
    typedef typename                          DOFMAPPER::OPTIONS   OPTIONS;
    typedef typename traitsMultiply<BASETYPE,DOFTYPE>::DIRECTTYPE  OUTTYPE;
    
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>     MESH2D;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>  CONNECT2D;
    
    typedef mesh2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>     FIELD_MESH;
    typedef connect2d<GEOSHAPE,PMAPTYPE,PMAPTYPE>  FIELD_CONNECT;
    typedef FETYPE  FIELD_FETYPE;
    typedef DOFTYPE FIELD_DOFTYPE;
    
    typedef pVect<FECARD,PMAPTYPE>  FECARDS;
    //@}
  
    /*! @name Internal links */ //@{
  public:
    Teuchos::RCP<MESH2D>    grid2d;
    Teuchos::RCP<CONNECT2D> connectGrid2d;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    Teuchos::RCP<DOFMAPPER> dofMapper;
    FEINTERFACE  feInterface;
    ELCARDFEEDER feeder;
    DOFVECT      dofVect;
    //@}
    
    /*! @name Internal flags */ //@{
  public:
    bool commLoaded, geometryLoaded, optionsLoaded, startupOk;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    feStaticField2d();
    feStaticField2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    feStaticField2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    feStaticField2d(const Teuchos::RCP<DOFMAPPER> & DofMapper);
    feStaticField2d(const DOFMAPPER & DofMapper);
    feStaticField2d(const feStaticField2d & Field);
    feStaticField2d & operator=(const feStaticField2d & Field);
    void clone(const feStaticField2d & Field);
    //@}
    
    /*! @name Mandatory links */ //@{
  public:
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    void setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d);
    void setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d);
    Teuchos::RCP<MESH2D>       getMesh2d() const;
    Teuchos::RCP<MESH2D>       getMesh() const;
    Teuchos::RCP<CONNECT2D>    getMeshConnect() const;
    const DOFMAPPER          & getDofMapper() const;
    const ELCARDFEEDER       & getElementFeeder() const;
    const FECARDS            & getFeCards() const;
          FECARDS            & getFeCards();
    const OPTIONS              getOptions() const;
    //@}
    
    /*! @name Optional links */ //@{
  public:
    void setFeCardL(const UInt & lid, const FECARD & FeCards);
    void setFeCardG(const UInt & gid, const FECARD & FeCards);
    void setFeCards(const FECARDS & FECards);
    void setOptions(const Teuchos::RCP<OPTIONS> & Options);
    void setOptions(OPTIONS & Options);
    //@}
    
    /*! @name Startup */ //@{
  public:
    void startup();
    //@}
    
     /*! @name Get flags */ //@{
  public:
    const bool & getGeoemtryLoaded() const;
    const bool & getOptionsLoaded() const;
    const bool & getStartupOk() const;
    //@}
    
    /*! @name Dof functions */ //@{
  public:
    void            setDofVect(const DOFVECT & DofVect);
    void            setListVect(const LISTVECT & ListVect);
    const DOFVECT & getDofVect() const;
    LISTVECT        getList() const;
    const LISTMAP & getListMap() const;
    //@}
    
    /*! @name Dof items functions */ //@{
  public:
    DOFTYPE       & getDofL(const UInt & lid);
    DOFTYPE       & getDofG(const UInt & gid); 
    const DOFTYPE & getDofL(const UInt & lid) const;
    const DOFTYPE & getDofG(const UInt & gid) const; 
    void            setDofL(const UInt & lid, const DOFTYPE & dof);
    void            setDofG(const UInt & gid, const DOFTYPE & dof);
    DOFTYPE &       operator() (const UInt & lid);
    const DOFTYPE & operator() (const UInt & lid) const;
    UInt            size() const;
    void            scale(const Real & scale);
    
    bool mapDofL(UInt & dofL, const DOFCARD & dofCard) const;
    bool mapDofG(UInt & dofG, const DOFCARD & dofCard) const;
    bool getDof(DOFTYPE & Val, const DOFCARD & dofCard) const;
    //@}
    
    /*! @name Evaluation functions */ //@{
  public:
    void evalL(const UInt & elL, const point3d & Y, OUTTYPE & val);
    void evalG(const UInt & elG, const point3d & Y, OUTTYPE & val);
    void gradL(const UInt & elL, const point3d & Y, OUTTYPE & gradX, OUTTYPE & gradY, OUTTYPE & gradZ);
    void gradG(const UInt & elG, const point3d & Y, OUTTYPE & gradX, OUTTYPE & gradY, OUTTYPE & gradZ);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SETTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField2d()
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  assert(feInterface.feTest());
  
  commLoaded     = false;
  geometryLoaded = false;
  optionsLoaded  = false;
  startupOk      = false;
  
  dofMapper = Teuchos::rcp(new DOFMAPPER);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField2d(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  assert(feInterface.feTest());
  
  commLoaded     = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  grid2d        = Grid2d;
  connectGrid2d = ConnedGrid2d;
  
  dofMapper = Teuchos::rcp(new DOFMAPPER);  
  dofMapper->setCommunicator(CommDev);
  dofMapper->setGeometry(Grid2d,ConnedGrid2d);
  feeder.setGeometry(Grid2d,ConnedGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField2d(communicator & CommDev, MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  assert(feInterface.feTest());
  
  commLoaded     = true;
  geometryLoaded = true;
  optionsLoaded  = false;
  startupOk      = false;
  
  grid2d        = Grid2d;
  connectGrid2d = ConnedGrid2d;
  
  dofMapper = Teuchos::rcp(new DOFMAPPER);
  dofMapper->setCommunicator(CommDev);
  dofMapper->setGeometry(Grid2d,ConnedGrid2d);
  feeder.setGeometry(Grid2d,ConnedGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField2d(const Teuchos::RCP<DOFMAPPER> & DofMapper)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  assert(DofMapper->getStartupOk());
  
  commLoaded     = true;
  geometryLoaded = true;
  optionsLoaded  = true;
  startupOk      = true;
  
  dofMapper     = DofMapper;
  grid2d        = DofMapper->getGrid2d();
  connectGrid2d = DofMapper->getConnectGrid2d();
   
  //Startup
  feeder.setGeometry(grid2d,connectGrid2d);  
  dofVect.resize(dofMapper->getNumDofsL());
  dofVect.setMap(dofMapper->getDofMap());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField2d(const DOFMAPPER & DofMapper)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  assert(DofMapper.getStartupOk());
  
  commLoaded     = true;
  geometryLoaded = true;
  optionsLoaded  = true;
  startupOk      = true;
  
  dofMapper     = Teuchos::rcp(new DOFMAPPER(DofMapper));
  grid2d        = DofMapper.getGrid2d();
  connectGrid2d = DofMapper.getConnectGrid2d();
   
  //Startup
  feeder.setGeometry(grid2d,connectGrid2d);  
  dofVect.resize(dofMapper->getNumDofsL());
  dofVect.setMap(dofMapper->getDofMap());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField2d(const feStaticField2d & Field)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
  
  commLoaded     = Field.commLoaded;
  geometryLoaded = Field.geometryLoaded;
  optionsLoaded  = Field.optionsLoaded;
  startupOk      = Field.startupOk;
  
  grid2d        = Field.grid2d;
  connectGrid2d = Field.connectGrid2d;
  
  dofMapper = Field.dofMapper;
  dofVect   = Field.dofVect;
  
  feeder.setGeometry(grid2d,connectGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::feStaticField2d &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
operator=(const feStaticField2d & Field)
{
  commLoaded     = Field.commLoaded;
  geometryLoaded = Field.geometryLoaded;
  optionsLoaded  = Field.optionsLoaded;
  startupOk      = Field.startupOk;
  
  grid2d        = Field.grid2d;
  connectGrid2d = Field.connectGrid2d;
  
  dofMapper = Field.dofMapper;
  dofVect   = Field.dofVect;
  
  feeder.setGeometry(grid2d,connectGrid2d);
  
  return(*this);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
clone(const feStaticField2d & Field)
{
  commLoaded     = Field.commLoaded;
  geometryLoaded = Field.geometryLoaded;
  optionsLoaded  = Field.optionsLoaded;
  startupOk      = Field.startupOk;
  
  grid2d        = Field.grid2d;
  connectGrid2d = Field.connectGrid2d;
  
  dofMapper = Field.dofMapper;
  dofVect   = Field.dofVect;
  
  feeder.setGeometry(grid2d,connectGrid2d);
}


//_________________________________________________________________________________________________
// MANDATORY LINKS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commLoaded = true;
  
  dofMapper->setCommunicator(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commLoaded = true;
  
  dofMapper->setCommunicator(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setGeometry(const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnedGrid2d)
{
  geometryLoaded = true;
  
  grid2d        = Grid2d;
  connectGrid2d = ConnedGrid2d;
  
  dofMapper->setGeometry(Grid2d,ConnedGrid2d);
  feeder.setGeometry(Grid2d,ConnedGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setGeometry(MESH2D & Grid2d, CONNECT2D & ConnedGrid2d)
{
  geometryLoaded = true;
  
  grid2d        = Teuchos::rcp(new MESH2D(Grid2d));
  connectGrid2d = Teuchos::rcp(new CONNECT2D(ConnedGrid2d));
  
  dofMapper->setGeometry(Grid2d,ConnedGrid2d);
  feeder.setGeometry(Grid2d,ConnedGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
Teuchos::RCP<typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::MESH2D>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getMesh2d() const
{
  assert(geometryLoaded);
  return(grid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
Teuchos::RCP<typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::MESH2D>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getMesh() const
{
  assert(geometryLoaded);
  return(grid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
Teuchos::RCP<typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::CONNECT2D>
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getMeshConnect() const
{
  assert(geometryLoaded);
  return(connectGrid2d);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::DOFMAPPER &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getDofMapper() const
{
  assert(startupOk);
  return(*dofMapper);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::ELCARDFEEDER &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getElementFeeder() const
{
  assert(startupOk);
  return(feeder);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::FECARDS &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getFeCards() const
{
  return(dofMapper->getFeCards());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::FECARDS &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getFeCards()
{
  return(dofMapper->getFeCards());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::OPTIONS
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getOptions() const
{
  return(*dofMapper->getOptions());
}



//_________________________________________________________________________________________________
// OPTIONAL LINKS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setFeCardL(const UInt & lid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  dofMapper->setFeCardL(lid,FeCards);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setFeCardG(const UInt & gid, const FECARD & FeCards)
{
  assert(geometryLoaded);
  dofMapper->setFeCardG(gid,FeCards);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setFeCards(const FECARDS & FECards)
{
  assert(geometryLoaded);
  dofMapper->setFeCards(FECards);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setOptions(const Teuchos::RCP<OPTIONS> & Options)
{
  optionsLoaded = true;
  
  dofMapper->setOptions(Options);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setOptions(OPTIONS & Options)
{
  optionsLoaded = true;
  
  dofMapper->setOptions(Options);
}



//_________________________________________________________________________________________________
// STARTUP
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
startup()
{
  //Logic flags
  assert(commLoaded);
  assert(geometryLoaded);
  assert(optionsLoaded);
  
  startupOk = true;
  
  //DofMapper startup
  dofMapper->startup();
  
  //Alloc the dof vector
  dofVect.resize(dofMapper->getNumDofsL());
  dofVect.setMap(dofMapper->getDofMap());
  dofVect.updateFinder();
}



//_________________________________________________________________________________________________
// GET FLAGS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const bool &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getGeoemtryLoaded() const
{
  return(geometryLoaded);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const bool &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getOptionsLoaded() const
{
  return(optionsLoaded);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const bool &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getStartupOk() const
{
  return(startupOk);
}



//_________________________________________________________________________________________________
// DOF FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setDofVect(const DOFVECT & DofVect)
{
  assert(startupOk);
  assert(DofVect.size() == dofVect.size());
  
  dofVect = DofVect;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setListVect(const LISTVECT & ListVect)
{
  assert(startupOk);
  
  traitsBasic<DOFTYPE> dofTrait; 
  UInt lid, listLid;
  DOFTYPE dof;
  
  pMap<PMAPTYPE> checkList(dofMapper->getListMap());
  
  for(UInt i=1; i <= dofVect.size(); ++i)
  {
    lid = dofVect.getMapL(i).getLid();
    
    assert(i == lid);
    
    for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)
    {
      for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
      {
	listLid = dofMapper->dofToListL(lid,I,J);
	
	assert(  !(checkList(listLid) != ListVect.getMapL(listLid)) );
	dofTrait.setIJ(I,J, ListVect(listLid), dof);
      }
    }
    
    dofVect(i) = dof;
  }
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const
typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::DOFVECT &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getDofVect() const
{
  assert(startupOk);
  return(dofVect);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::LISTVECT
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getList() const
{
  assert(startupOk);
  
  typedef pVect<Real,PMAPTYPE>  LISTVECT; 
  
  LISTVECT listVect(dofMapper->getSizeListL());
  listVect.setMap(dofMapper->getListMap());
  
  traitsBasic<DOFTYPE> dofTrait; 
  UInt lid, gid, listLid, listGid;
  
  for(UInt i=1; i <= dofVect.size(); ++i)
  {
    lid = dofVect.getMapL(i).getLid();
    gid = dofVect.getMapL(i).getGid();
    
    assert(i == lid);
    
    for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)
    {
      for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
      {
	listLid = dofMapper->dofToListL(lid,I,J);
	listGid = dofMapper->dofToListG(gid,I,J);
	
	assert(listVect.getMapL(listLid).getGid() == listGid);
	listVect(listLid) = dofTrait.getIJ(I,J,dofVect.getDataL(i));
      }
    }
  }
  
  return(listVect);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const
typename feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::LISTMAP &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getListMap() const
{
  assert(startupOk);
  return(dofMapper->getListMap());
}



//_________________________________________________________________________________________________
// DOF ITEMS FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
DOFTYPE &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getDofL(const UInt & lid)
{
  assert(lid >= 1);
  assert(lid <= dofVect.size());
  return(dofVect.getL(lid));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
DOFTYPE &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getDofG(const UInt & gid)
{
  assert(dofVect.isG(gid));
  return(dofVect.getG(gid));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const DOFTYPE &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getDofL(const UInt & lid) const
{
  assert(lid >= 1);
  assert(lid <= dofVect.size());
  return(dofVect.getL(lid));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const DOFTYPE &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getDofG(const UInt & gid) const
{
  assert(dofVect.isG(gid));
  return(dofVect.getG(gid));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void 
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setDofL(const UInt & lid, const DOFTYPE & dof)
{
  assert(lid >= 1);
  assert(lid <= dofVect.size());
  
  dofVect.getL(lid) = dof;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
setDofG(const UInt & gid, const DOFTYPE & dof)
{
  assert(dofVect.isG(gid));
  
  dofVect.getG(gid) = dof;
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
DOFTYPE &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
operator() (const UInt & lid)
{
  assert(lid >= 1);
  assert(lid <= dofVect.size());
  
  return(dofVect(lid));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
const DOFTYPE &
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
operator() (const UInt & lid) const
{
  assert(lid >= 1);
  assert(lid <= dofVect.size());
  
  return(dofVect(lid));
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
UInt
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
size() const
{
  return(dofVect.size());
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
scale(const Real & scale)
{
  for(UInt i=1; i <= dofVect.size(); ++i)
  { dofVect(i) *= scale; }
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
bool
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
mapDofL(UInt & dofL, const DOFCARD & dofCard) const
{
  assert(startupOk);
  bool success = dofMapper->isActive(dofCard);
  if(success) { dofL = dofMapper->mapDofL(dofCard); }
  
  return(success);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
bool
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
mapDofG(UInt & dofG, const DOFCARD & dofCard) const
{
  assert(startupOk);
  bool success = dofMapper->isActive(dofCard);
  if(success) { dofG = dofMapper->mapDofG(dofCard); }
  
  return(success);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
bool
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
getDof(DOFTYPE & Val, const DOFCARD & dofCard) const
{
  assert(startupOk);
  bool success = dofMapper->isActive(dofCard);
  
  if(success)
  {
    UInt dofL = dofMapper->mapDofL(dofCard);
          Val = getDofL(dofL);
  }
  
  return(success);
}



//_________________________________________________________________________________________________
// EVAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
evalL(const UInt & elL, const point3d & Y, OUTTYPE & val)
{
  //Assert-----------------------------------------------------------
  assert(startupOk);
  assert(elL > 0);
  assert(elL <= grid2d->getNumElements());
  
  //Compute the number of basis--------------------------------------
  DOFCARD dofCard;
  dofCard.setLocalElId(elL);
  
  UInt numBasis = dofMapper->isActive(dofCard) * FETYPE::numBasis;
  
  //Basis evaluation-------------------------------------------------
  sVect<BASETYPE> basisVal(FETYPE::numBasis);
  
  feInterface.setCards(dofMapper->getFeCards().getL(elL), feeder.getCardLocal(elL) );
  feInterface.globalEval(Y,basisVal);
  
  //Evaluation------------------------------------------------------
  UInt dofLid;
  val = traitsBasic<OUTTYPE>::getZero();
  
  for(UInt i=1; i <= numBasis; ++i)
  {
    dofCard = feInterface.getDofCard(i);
    dofCard.setLocalElId(elL);
    
    dofLid = dofMapper->mapDofL(dofCard);
    
    val += traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisVal(i), dofVect(dofLid));
  }
  
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
evalG(const UInt & elG, const point3d & Y, OUTTYPE & val)
{
  assert(startupOk);
  
  UInt elL = grid2d->getElements().getRowMapG(elG).getLid();
  evalL(elL,Y,val);
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
gradL(const UInt & elL, const point3d & Y, OUTTYPE & gradX, OUTTYPE & gradY, OUTTYPE & gradZ)
{
  //Assert-----------------------------------------------------------
  assert(startupOk);
  assert(elL > 0);
  assert(elL <= grid2d->getNumElements());
  
  //Compute the number of basis--------------------------------------
  DOFCARD dofCard;
  dofCard.setLocalElId(elL);
  
  UInt numBasis = dofMapper->isActive(dofCard) * FETYPE::numBasis;
  
  //Basis evaluation-------------------------------------------------
  sVect<BASETYPE> basisGradX(FETYPE::numBasis);
  sVect<BASETYPE> basisGradY(FETYPE::numBasis);
  sVect<BASETYPE> basisGradZ(FETYPE::numBasis);
  
  feInterface.setCards(dofMapper->getFeCards().getL(elL), feeder.getCardLocal(elL) );
  feInterface.globalGrad(Y, basisGradX, basisGradY, basisGradZ);
  
  //Evaluation------------------------------------------------------
  UInt dofLid;
  gradX = traitsBasic<OUTTYPE>::getZero();
  gradY = traitsBasic<OUTTYPE>::getZero();
  gradZ = traitsBasic<OUTTYPE>::getZero();
  
  for(UInt i=1; i <= numBasis; ++i)
  {
    dofCard = feInterface.getDofCard(i);
    dofCard.setLocalElId(elL);
    
    dofLid = dofMapper->mapDofL(dofCard);
    
    gradX += traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisGradX(i), dofVect(dofLid));
    gradY += traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisGradY(i), dofVect(dofLid));
    gradZ += traitsMultiply<BASETYPE,DOFTYPE>::multiply(basisGradZ(i), dofVect(dofLid));
  }
}

template<typename FETYPE, typename DOFTYPE, dms2d_order ORDER, dms2d_mode MODE>
void
feStaticField2d<FETYPE,DOFTYPE,ORDER,MODE>::
gradG(const UInt & elG, const point3d & Y, OUTTYPE & gradX, OUTTYPE & gradY, OUTTYPE & gradZ)
{
  assert(startupOk);
  
  UInt elL = grid2d->getElements().getRowMapG(elG).getLid();
  gradL(elL,Y,gradX,gradY,gradZ);
}
    
#endif
