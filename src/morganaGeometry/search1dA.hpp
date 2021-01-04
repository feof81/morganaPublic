/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef SEARCH1DA_HPP
#define SEARCH1DA_HPP

#include <boost/mpi/collectives.hpp>
#include <iostream>
#include <cstdlib>

#include "pMapItemSendRecv.h"
#include "localSearch1dA.hpp"

using namespace std;
namespace mpi = boost::mpi;


/*! Blobal search function, 1d meshes */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class search1dA : public localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>
{
    /*! @name Typedefs */ //@{
  public:
    typedef geoElement<GEOSHAPE>            GEOELEMENT1D;
    typedef mesh1d<GEOSHAPE,ELMAP,NODEMAP>  MESH1D;
    typedef searchData<ELMAP>               SEARCHDATA;
    //@}

    /*! @name Static variables */ //@{
  public:
    static const UInt nDim = 1;
    //@}
  
    /*! @name Internal flags */ //@{
  public:
    bool startupGlobal;
    bool commDevLoaded;
    //@}
    
    /*! @name Internal structures */ //@{
  public:
    UInt globalQuery;
    UInt numFeasibProc;
    //@}
    
    /*! @name Local internal data */ //@{
  public:
    sVect<UInt>    globalNumElements;
    sVect<point3d> globalBoxMax, globalBoxMin; //! global bounding box structure 
    //@}
    
    /*! @name Links */ //@{
  public:
    Teuchos::RCP<const communicator> commDev;
    //@}

    /*! @name Constructor and setting functions */ //@{
  public:
    search1dA();
    search1dA(const Teuchos::RCP<const communicator> & CommDev);
    search1dA(const communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<const communicator> & CommDev);
    void setCommunicator(const communicator & CommDev);
    //@}
    
    /*! @name Finding functions */ //@{
 public:
    /*! The bounding boxes of all the \c pids are gathered */
    void globalInit();
   
    /*! Computes the minimum distance from the point \c P to the box defined by \c Pmin and \c Pmax */
    Real minDistance(const point3d & P, const point3d & Pmin, const point3d & Pmax) const;
    
    /*! Computes the maximum distance from the point \c P to the box defined by \c Pmin and \c Pmax */
    Real maxDistance(const point3d & P, const point3d & Pmin, const point3d & Pmax) const;
   
    /*! Given a point \c P computes all the pids that may have a minimum distance.
    It computes for each pid the upper and lower bound. The upper bound (dMax) of the \c pid that has the minimum
    lower bound is identified. If a \c pid has a lower bound greated than dMax then is discarded. */
    sVect<UInt> globalMatchingPids(const point3d & P) const;
    
    /*! Every \c pid should specify a nodes vector \c Pvect. The algorithm tries to find them locally.
    If the nearest point in the mesh does not coincides with the requested one the algorithm selects
    the most promising \c pids and send them the point cordinate. Each \c pid returns the research data
    along with the distance. The one with the nearest data is selected. */
    void findGlobal(pVect<point3d,NODEMAP> & Pvect, pVect<SEARCHDATA,NODEMAP> & DataVect);
    //@}
    
    /*! @name Printing functions */ //@{
 public:
    void printPerformances(const UInt & printPid) const;
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
search1dA() : localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::localSearch1dA()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  startupGlobal = false;
  commDevLoaded = false;
  
  globalQuery   = 0;
  numFeasibProc = 0;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
search1dA(const Teuchos::RCP<const communicator> & CommDev) : localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::localSearch1dA()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  startupGlobal = false;
  commDevLoaded = true;
  
  globalQuery   = 0;
  numFeasibProc = 0;
  
  commDev = CommDev;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
search1dA(const communicator & CommDev) : localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::localSearch1dA()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
  
  startupGlobal = false;
  commDevLoaded = true;
  
  globalQuery   = 0;
  numFeasibProc = 0;
  
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const Teuchos::RCP<const communicator> & CommDev)
{
  commDev = CommDev;
  commDevLoaded = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
setCommunicator(const communicator & CommDev)
{
  commDev = Teuchos::rcpFromRef(CommDev);
  commDevLoaded = true;
}



//_________________________________________________________________________________________________
// PRINTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
globalInit()
{
  typedef localSearch1dA<GEOSHAPE,ELMAP,NODEMAP> LOCALSEARCH;
  
  assert(commDevLoaded);
  assert(LOCALSEARCH::startupLocal);  
  
  startupGlobal = true;
  
  UInt numPids = commDev->size();
  
  sVect<Real> maxVectX(numPids), maxVectY(numPids), maxVectZ(numPids);
  sVect<Real> minVectX(numPids), minVectY(numPids), minVectZ(numPids);
  
  all_gather(*commDev, LOCALSEARCH::Xmax, maxVectX);
  all_gather(*commDev, LOCALSEARCH::Ymax, maxVectY);
  all_gather(*commDev, LOCALSEARCH::Zmax, maxVectZ);
  
  all_gather(*commDev, LOCALSEARCH::Xmin, minVectX);
  all_gather(*commDev, LOCALSEARCH::Ymin, minVectY);
  all_gather(*commDev, LOCALSEARCH::Zmin, minVectZ);
  
  globalNumElements.resize(numPids);
  all_gather(*commDev, LOCALSEARCH::grid1d->getNumElements(), globalNumElements);

  globalBoxMax.resize(numPids);
  globalBoxMin.resize(numPids);
  
  for(UInt i=1; i <= numPids; ++i)
  {
    globalBoxMax(i).setX(maxVectX(i));
    globalBoxMax(i).setY(maxVectY(i));
    globalBoxMax(i).setZ(maxVectZ(i));
    
    globalBoxMin(i).setX(minVectX(i));
    globalBoxMin(i).setY(minVectY(i));
    globalBoxMin(i).setZ(minVectZ(i));
  }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Real
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
minDistance(const point3d & P, const point3d & Pmin, const point3d & Pmax) const
{
  Real deltaX = (P.getX() >= Pmax.getX()) * (P.getX() - Pmax.getX())  +  (P.getX() <= Pmin.getX()) * (Pmin.getX() - P.getX());
  Real deltaY = (P.getY() >= Pmax.getY()) * (P.getY() - Pmax.getY())  +  (P.getY() <= Pmin.getY()) * (Pmin.getY() - P.getY());
  Real deltaZ = (P.getZ() >= Pmax.getZ()) * (P.getZ() - Pmax.getZ())  +  (P.getZ() <= Pmin.getZ()) * (Pmin.getZ() - P.getZ());
  
  return( sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ) );
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
Real
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
maxDistance(const point3d & P, const point3d & Pmin, const point3d & Pmax) const
{
  Real deltaX = max( abs(P.getX() - Pmax.getX()) , abs(P.getX() - Pmin.getX()) );
  Real deltaY = max( abs(P.getY() - Pmax.getY()) , abs(P.getY() - Pmin.getY()) );
  Real deltaZ = max( abs(P.getZ() - Pmax.getZ()) , abs(P.getZ() - Pmin.getZ()) );
  
  return( sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ) );
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
sVect<UInt>
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
globalMatchingPids(const point3d & P) const
{
  //Typedefs
  typedef localSearch1dA<GEOSHAPE,ELMAP,NODEMAP> LOCALSEARCH;
  
  //Asserts
  assert(startupGlobal);
  
  //Allocation
  sVect<Real> dMinVect(globalBoxMax.size());
  Real dMin = std::numeric_limits<double>::infinity();
  Real dMax = 0.0;
  Real dTemp;
  point3d Pmax, Pmin;
  
  //Initialization
  assert(globalBoxMin.size() == globalBoxMax.size());
  assert(globalBoxMin.size() >= 1);
  
  /*Pmin = globalBoxMin(1);
  Pmax = globalBoxMax(1);
    
  dMin = minDistance(P,Pmin,Pmax);
  dMax = maxDistance(P,Pmin,Pmax);*/
  
  //Cycle
  for(UInt i=1; i <= globalBoxMax.size(); ++i)
  {
    Pmin = globalBoxMin(i);
    Pmax = globalBoxMax(i);
    
    dTemp       = minDistance(P,Pmin,Pmax);
    dMinVect(i) = dTemp;
    
    if( (dTemp < dMin) && (globalNumElements(i) != 0) )
    {
      dMin = dTemp;
      dMax = maxDistance(P,Pmin,Pmax);
    }
  }
  
  //Matching pids determination
  sVect<UInt> out;
  
  for(UInt i=1; i <= dMinVect.size(); ++i)
  {
    if(dMinVect(i) <= dMax)
    {
      out.push_back(i-1);
    }
  }
  
  return(out);
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
findGlobal(pVect<point3d,NODEMAP> & Pvect, pVect<SEARCHDATA,NODEMAP> & DataVect)
{
  typedef localSearch1dA<GEOSHAPE,ELMAP,NODEMAP> LOCALSEARCH;
  
  assert(startupGlobal);
  
  
  //Local searching________________________________________________________________________________
  DataVect.resize(Pvect.size());
  DataVect.setMap(Pvect.getMapRef());
  
  UInt el;
  UInt pid = commDev->rank();
  
  for(UInt i=1; i <= Pvect.size(); ++i)
  {
    DataVect(i) = localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::findLocal(Pvect(i));
  }
  DataVect.updateFinder();
  
  
  //Global sending map creation____________________________________________________________________
  sVect<UInt> matchingPids;
  pMapItemSendRecv sendItem;
  pMap<pMapItemSendRecv> commMap;
  
  for(UInt i=1; i <= Pvect.size(); ++i)
  {
    if( (DataVect(i).getDistance() >= geoToll) || (!DataVect(i).getFound()))
    {      
      matchingPids = globalMatchingPids(Pvect(i));
      
      for(UInt j=1; j <= matchingPids.size(); ++j)
      {
	if(matchingPids(j) != pid)
	{
	  assert(int(matchingPids(j)) < commDev->size());
	  
	  sendItem.setPid(pid);
	  sendItem.setLid(Pvect.getMapL(i).getLid());
	  sendItem.setGid(Pvect.getMapL(i).getGid());
	  sendItem.setSid(pid);
	  sendItem.setRid(matchingPids(j));
	  
	  commMap.push_back(sendItem);
	}
      }
      
      globalQuery++;
      numFeasibProc += matchingPids.size();
    }
  }
  
  //Communication__________________________________________________________________________________
  pVect<point3d,NODEMAP> Pother;
  
  pVectComm<point3d,NODEMAP> pointsComm(commDev);
  pointsComm.sendRecv(commMap,Pvect,Pother);
  
  
  //Searching the othed pids nodes_________________________________________________________________
  pVect<SEARCHDATA,NODEMAP> dataOther(Pother.size());
  dataOther.setMap(Pother.getMapRef());
  
  for(UInt i=1; i <= Pother.size(); ++i)
  {    
    dataOther(i) = localSearch1dA<GEOSHAPE,ELMAP,NODEMAP>::findLocal(Pother(i));
    
    el = dataOther(i).getElMap().getLid();
    assert(LOCALSEARCH::grid1d->getElements().getRowMapL(el).getPid() == pid);
  }
  
  
  //Transmit back__________________________________________________________________________________
  pVectComm<SEARCHDATA,NODEMAP> otherDataComm(commDev);
  otherDataComm.vectorPid(dataOther);
  
  //Find the minimum distance______________________________________________________________________
  UInt gid;
  Real d;
  
  for(UInt i=1; i <= dataOther.size(); ++i)
  {
    gid = dataOther.getMapL(i).getGid();
    d   = dataOther.getL(i).getDistance();
    
    assert(DataVect.isG(gid));
    
    if((d < DataVect.getG(gid).getDistance()) || (!DataVect.getG(gid).getFound()))
    {
      DataVect.getG(gid) = dataOther(i);
    }
  }
}


//________________________________________________________________________________________________
// PRINTING FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
search1dA<GEOSHAPE,ELMAP,NODEMAP>::
printPerformances(const UInt & printPid) const
{
  typedef localSearch1dA<GEOSHAPE,ELMAP,NODEMAP> LOCALSEARCH;
  
  int totalFeasibProc;
  int totalFeasibleEl;
  int totalLocalQuery;
  int totalGlobalQuery;
  
  reduce(*commDev, int(numFeasibProc),     totalFeasibProc,  std::plus<int>(), printPid);
  reduce(*commDev, int(LOCALSEARCH::numFeasibElements), totalFeasibleEl,  std::plus<int>(), printPid);
  reduce(*commDev, int(LOCALSEARCH::localQuery),        totalLocalQuery,  std::plus<int>(), printPid);
  reduce(*commDev, int(globalQuery),       totalGlobalQuery, std::plus<int>(), printPid);
  
  if(commDev->rank() == int(printPid))
  {
    cout << "Mean feasible elements found for every point : " << Real(totalFeasibleEl) / Real(totalLocalQuery) << endl;
    cout << "Mean pids found for every external point     : " << Real(totalFeasibProc) / Real(totalGlobalQuery) << endl;
  }
}

#endif
