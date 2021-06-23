/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef RESUMEHDF5_H
#define RESUMEHDF5_H

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "EpetraExt_HDF5.h"
#include "EpetraExt_DistArray.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "traitsResumeHDF5.h"

#include "pMap.hpp"
#include "pVect.hpp"
#include "pGraph.hpp"

#include "pMapManip.hpp"
#include "pVectManip.hpp"

#include "pMapComm.hpp"
#include "pVectComm.hpp"

#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"
#include "mesh1dGlobalManip.hpp"
#include "mesh2dGlobalManip.hpp"
#include "mesh3dGlobalManip.hpp"

#include "dofMapStatic1d_options.h"
#include "dofMapStatic2d_options.h"
#include "dofMapStatic3d_options.h"
#include "dofMapStatic1d.hpp"
#include "dofMapStatic2d.hpp"
#include "dofMapStatic3d.hpp"


/*! Class to save to file and load from file run-time data */
class resumeHDF5
{
  /*! @name Links and flags */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructor and functions */ //@{
  public:
    resumeHDF5();
    resumeHDF5(const Teuchos::RCP<communicator> & CommDev);
    resumeHDF5(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
     /*! @name Support functions for arrays */ //@{
  public:
    template<typename DATA>
    void serializeArray(const sVect<sVect<DATA> > & inArray,
                              sVect<UInt>         & outIndex,
                              sVect<DATA>         & outSerial);
    
    template<typename DATA>
    void constructArray(const sVect<UInt>         & inIndex,
                        const sVect<DATA>         & inSerial,
                              sVect<sVect<DATA> > & outArray);
    //@}
    
    /*! @name Convert to distArray */ //@{
  public:
    template<typename DATA>
    EpetraExt::DistArray<Real> convertToDistArray(const sVect<DATA> & morganaVect);
    
    template<typename ITEM>
    EpetraExt::DistArray<Real> convertToDistArray(const pMap<ITEM> & morganaMap);
    
    template<typename DATA, typename MAP>
    EpetraExt::DistArray<Real> convertToDistArray(const pVect<DATA,MAP> & morganaVect);
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    EpetraExt::DistArray<Real> convertToDistArray(const pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph);
    //@}
    
    /*! @name Convert from distArray */ //@{
  public:
    template<typename DATA> 
    void convertFromDistArray(EpetraExt::DistArray<Real> & distArray, sVect<DATA> & morganaVect);
    
    template<typename ITEM> 
    void convertFromDistArray(EpetraExt::DistArray<Real> & distArray, pMap<ITEM> & morganaMap);
    
    template<typename DATA, typename MAP>
    void convertFromDistArray(EpetraExt::DistArray<Real> & distArray, pVect<DATA,MAP> & morganaVect);
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void convertFromDistArray(EpetraExt::DistArray<Real> & distArray, pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph);
    //@}
    
    /*! @name Print to file - Containers */ //@{
  public:
    template<typename DATA>
    void packToFile(const string          & vectName,
                          EpetraExt::HDF5 & hdf5,
                    const sVect<DATA> & morganaVect);
    
    template<typename DATA>
    void printToFile(const string & s,
                     const sVect<DATA> & morganaVect);
    
    template<typename DATA>
    void packToFile(const string & arrayName_index,
                    const string & arrayName_serialized,
                          EpetraExt::HDF5 & hdf5,
                    const sVect<sVect<DATA> > & morganaArray);
    
    template<typename DATA>
    void printToFile(const string & s,
                     const sVect<sVect<DATA> > & morganaArray);
    
    template<typename ITEM> 
    void packToFile(const string          & mapName,
                          EpetraExt::HDF5 & hdf5,
                    const pMap<ITEM> & morganaMap);
    
    template<typename ITEM> 
    void printToFile(const string & s,
                     const pMap<ITEM> & morganaMap);
    
    template<typename DATA, typename MAP>
    void packToFile(const string          & vectName_map,
                    const string          & vectName_data,
                          EpetraExt::HDF5 & hdf5,
                    const pVect<DATA,MAP> & morganaVect);
    
    template<typename DATA, typename MAP>
    void printToFile(const string & s,
                     const pVect<DATA,MAP> & morganaVect);
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void packToFile(const string          & graphName_rowMap,
                    const string          & graphName_colMap,
                    const string          & graphName_data,
                          EpetraExt::HDF5 & hdf5,
                    const pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph);
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void printToFile(const string & s,
                     const pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph);
    //@}
    
    /*! @name Load from file - Containers */ //@{
  public:
    template<typename DATA>
    void extractFromFile(const string          & vectName,
                               EpetraExt::HDF5 & hdf5,
                               sVect<DATA>     & morganaVect);
    
    template<typename DATA>
    void loadFromFile(const string & s,
                            sVect<DATA> & morganaVect);
    
    template<typename DATA>
    void extractFromFile(const string              & arrayName_index,
                         const string              & arrayName_serialized,
                               EpetraExt::HDF5     & hdf5,
                               sVect<sVect<DATA> > & morganaArray);
    
    template<typename DATA>
    void loadFromFile(const string & s,
                            sVect<sVect<DATA> > & morganaArray);
    
    template<typename ITEM> 
    void extractFromFile(const string          & mapName,
                               EpetraExt::HDF5 & hdf5,
                               pMap<ITEM>      & morganaMap);
    
    template<typename ITEM> 
    void loadFromFile(const string & s,
                            pMap<ITEM> & morganaMap);
    
    template<typename DATA, typename MAP>
    void extractFromFile(const string          & vectName_map,
                         const string          & vectName_data,
                               EpetraExt::HDF5 & hdf5,
                               pVect<DATA,MAP> & morganaVect);
    
    template<typename DATA, typename MAP>
    void loadFromFile(const string & s,
                            pVect<DATA,MAP> & morganaVect);
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void extractFromFile(const string          & graphName_rowMap,
                         const string          & graphName_colMap,
                         const string          & graphName_data,
                               EpetraExt::HDF5 & hdf5,
                               pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph);
    
    template<typename ITEM, typename ROWMAP, typename COLMAP>
    void loadFromFile(const string & s,
                            pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph);
    //@}
    
    /*! @name Print to file - Geometry */ //@{
  public:
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void packToFile(const string          & meshName,
                          EpetraExt::HDF5 & hdf5,
                    const mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printToFile(const string & s,
                     const mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void packToFile(const string          & meshName,
                          EpetraExt::HDF5 & hdf5,
                    const mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printToFile(const string & s,
                     const mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void packToFile(const string          & meshName,
                          EpetraExt::HDF5 & hdf5,
                    const mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void printToFile(const string & s,
                     const mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    //@}
    
    /*! @name Load from file - Geometry */ //@{
  public:
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void extractFromFile(const string          & meshName,
                               EpetraExt::HDF5 & hdf5,
                               mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadFromFile(const string & s,
                            mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void extractFromFile(const string          & meshName,
                               EpetraExt::HDF5 & hdf5,
                               mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadFromFile(const string & s,
                            mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void extractFromFile(const string          & meshName,
                               EpetraExt::HDF5 & hdf5,
                               mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    
    template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
    void loadFromFile(const string & s,
                            mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh);
    //@}
    
    /*! @name Load and Write from file - Trilinos */ //@{
  public:
    void printToFile(const string & s,
                     const Epetra_CrsMatrix & A);
    
    void loadFromFile(const string & s,
                            Epetra_CrsMatrix & A);
    //@}
};



//_________________________________________________________________________________________________
// SUPPORT FUNCTIONS FOR ARRAYS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
resumeHDF5::
serializeArray(const sVect<sVect<DATA> > & inArray,
                     sVect<UInt>         & outIndex,
                     sVect<DATA>         & outSerial)
{
  outIndex.resize(inArray.size());
  outSerial.clear();
  
  for(UInt i=1; i <= inArray.size(); ++i)
  { outIndex(i) = inArray(i).size(); }
  
  for(UInt i=1; i <= inArray.size(); ++i)
  {
    for(UInt j=1; j <= inArray(i).size(); ++j)
    {
      outSerial.push_back(inArray(i)(j));
    }
  }
}
    
template<typename DATA>
void
resumeHDF5::
constructArray(const sVect<UInt> & inIndex,
               const sVect<DATA> & inSerial,
                     sVect<sVect<DATA> > & outArray)
{
  UInt k=1;
  outArray.resize(inIndex.size());
  
  for(UInt i=1; i <= inIndex.size(); ++i)
  {
    outArray(i).resize(inIndex(i));
    
    for(UInt j=1; j <= inIndex(i); ++j)
    {
      outArray(i)(j) = inSerial(k);
      ++k;
    }
  }
}


//_________________________________________________________________________________________________
// CONVERT TO DISTRIBUTED ARRAY
//-------------------------------------------------------------------------------------------------
template<typename DATA>
EpetraExt::DistArray<Real>
resumeHDF5::
convertToDistArray(const sVect<DATA> & morganaVect)
{
  assert(commDevLoaded);
  
  //EpetraMap--------------------------------------------------------
  Epetra_MpiComm  epetraComm(*commDev);
  Epetra_Map epetraMap(-1,morganaVect.size(),0, epetraComm);
  
  //Dist array-------------------------------------------------------
  EpetraExt::DistArray<Real> outArray(epetraMap, traitsResumeHDF5<DATA>::size());
  
  for(UInt i=1; i <= morganaVect.size(); ++i)
  {
    for(UInt j=1; j <= traitsResumeHDF5<DATA>::size(); ++j)
    {
      outArray(i-1,j-1) = traitsResumeHDF5<DATA>::getValue(morganaVect(i),j);
    }
  }
  
  return(outArray);
}

template<typename ITEM> 
EpetraExt::DistArray<Real> 
resumeHDF5::
convertToDistArray(const pMap<ITEM> & morganaMap)
{
  assert(commDevLoaded);
  
  //EpetraMap--------------------------------------------------------
  Epetra_MpiComm  epetraComm(*commDev);
  Epetra_Map epetraMap(-1,morganaMap.size(),0, epetraComm);
  
  //Dist array-------------------------------------------------------
  EpetraExt::DistArray<Real> outArray(epetraMap, traitsResumeHDF5<ITEM>::size());
  
  for(UInt i=1; i <= morganaMap.size(); ++i)
  {
    for(UInt j=1; j <= traitsResumeHDF5<ITEM>::size(); ++j)
    {
      outArray(i-1,j-1) = Real( traitsResumeHDF5<ITEM>::getValue(morganaMap(i),j) );
    }
  }
  
  return(outArray);
}

template<typename DATA, typename MAP>
EpetraExt::DistArray<Real>
resumeHDF5::
convertToDistArray(const pVect<DATA,MAP> & morganaVect)
{
  assert(commDevLoaded);
  
  //EpetraMap--------------------------------------------------------
  Epetra_MpiComm  epetraComm(*commDev);
  Epetra_Map epetraMap(-1,morganaVect.size(),0, epetraComm);
  
  //Dist array-------------------------------------------------------
  EpetraExt::DistArray<Real> outArray(epetraMap, traitsResumeHDF5<DATA>::size());
  
  for(UInt i=1; i <= morganaVect.size(); ++i)
  {
    for(UInt j=1; j <= traitsResumeHDF5<DATA>::size(); ++j)
    {
      outArray(i-1,j-1) = traitsResumeHDF5<DATA>::getValue(morganaVect(i),j);
    }
  }
  
  assert(outArray.Map().LinearMap());
  
  return(outArray);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
EpetraExt::DistArray<Real>
resumeHDF5::
convertToDistArray(const pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph)
{
  assert(commDevLoaded);
  
  //EpetraMap--------------------------------------------------------
  Epetra_MpiComm  epetraComm(*commDev);
  Epetra_Map epetraMap(-1,morganaGraph.size(),0, epetraComm);
  
  //Dist array-------------------------------------------------------
  EpetraExt::DistArray<Real> outArray(epetraMap, traitsResumeHDF5<ITEM>::size());
  
  for(UInt i=1; i <= morganaGraph.size(); ++i)
  {
    for(UInt j=1; j <= traitsResumeHDF5<ITEM>::size(); ++j)
    {
      outArray(i-1,j-1) = traitsResumeHDF5<ITEM>::getValue(morganaGraph(i),j);
    }
  }
  
  assert(outArray.Map().LinearMap());
  
  return(outArray);
}


//_________________________________________________________________________________________________
// CONVERT FROM DISTRIBUTED ARRAY
//-------------------------------------------------------------------------------------------------
template<typename DATA> 
void
resumeHDF5::
convertFromDistArray(EpetraExt::DistArray<Real> & distArray, sVect<DATA> & morganaVect)
{
  assert(commDevLoaded);
  
  //Array size-------------------------------------------------------
  UInt row = distArray.MyLength();
  UInt col = distArray.RowSize();
  
  //Create map-------------------------------------------------------
  sVect<DATA> vect(row);
  DATA dof;
  
  for(UInt i=1; i <= row; ++i)
  {
    for(UInt j=1; j <= col; ++j)
    { traitsResumeHDF5<DATA>::setValue(distArray(i-1,j-1),j,dof); }
        
    vect(i) = dof;
  }
  
  morganaVect = vect;
}

template<typename ITEM> 
void
resumeHDF5::
convertFromDistArray(EpetraExt::DistArray<Real> & distArray, pMap<ITEM> & morganaMap)
{  
  assert(commDevLoaded);
  
  //Array size-------------------------------------------------------
  UInt row = distArray.MyLength();
  UInt col = distArray.RowSize();
  
  //Create map-------------------------------------------------------
  pMap<ITEM> map(row);
  ITEM item;
  
  for(UInt i=1; i <= row; ++i)
  {
    for(UInt j=1; j <= col; ++j)
    { traitsResumeHDF5<ITEM>::setValue(distArray(i-1,j-1),j,item); }
        
    map(i) = item;
  }
  
  morganaMap = map;
}

template<typename DATA, typename MAP>
void
resumeHDF5::
convertFromDistArray(EpetraExt::DistArray<Real> & distArray, pVect<DATA,MAP> & morganaVect)
{
  assert(commDevLoaded);
  
  //Array size-------------------------------------------------------
  UInt row = distArray.MyLength();
  UInt col = distArray.RowSize();
  
  //Create map-------------------------------------------------------
  pVect<DATA,MAP> vect(row);
  DATA dof;
  
  for(UInt i=1; i <= row; ++i)
  {
    for(UInt j=1; j <= col; ++j)
    { traitsResumeHDF5<DATA>::setValue(distArray(i-1,j-1),j,dof); }
        
    vect(i) = dof;
  }
  
  morganaVect = vect;
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
resumeHDF5::
convertFromDistArray(EpetraExt::DistArray<Real> & distArray, pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph)
{
  assert(commDevLoaded);
  
  //Array size-------------------------------------------------------
  UInt row = distArray.MyLength();
  UInt col = distArray.RowSize();
  
  //Create map-------------------------------------------------------
  pGraph<ITEM,ROWMAP,COLMAP> graph(row);
  ITEM dof(col-1);
  
  for(UInt i=1; i <= row; ++i)
  {
    for(UInt j=1; j <= col; ++j)
    { traitsResumeHDF5<ITEM>::setValue(distArray(i-1,j-1),j,dof); }
        
    graph(i) = dof;
  }
  
  morganaGraph = graph;
}


//_________________________________________________________________________________________________
// PRINT TO FILE - CONTAINERS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
resumeHDF5::
packToFile(const string          & vectName,
                 EpetraExt::HDF5 & hdf5,
           const sVect<DATA>     & morganaVect)
{
  assert(commDevLoaded);
  
  EpetraExt::DistArray<Real> distArray = convertToDistArray(morganaVect);
  
  commDev->barrier();
  hdf5.Write(vectName,distArray);
}

template<typename DATA> 
void
resumeHDF5::
printToFile(const string & s,
            const sVect<DATA> & morganaVect)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  
  hdf5.Create(s + ".h5");
  packToFile("sVect",hdf5,morganaVect);
  hdf5.Close();
}

template<typename DATA>
void
resumeHDF5::
packToFile(const string          & arrayName_index,
           const string          & arrayName_serialized,
                 EpetraExt::HDF5 & hdf5,
           const sVect<sVect<DATA> > & morganaArray)
{
  assert(commDevLoaded);
  
  sVect<UInt> outIndex;
  sVect<DATA> outSerial;
  serializeArray(morganaArray,outIndex,outSerial);
  
  commDev->barrier();
  packToFile(arrayName_index, hdf5, outIndex);
  
  commDev->barrier();
  packToFile(arrayName_serialized, hdf5, outSerial);
}

template<typename DATA>
void
resumeHDF5::
printToFile(const string & s,
            const sVect<sVect<DATA> > & morganaArray)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  
  hdf5.Create(s + ".h5");  
  packToFile("index",
	     "serialized",
	     hdf5,
	     morganaArray);
  hdf5.Close();
}

template<typename ITEM> 
void
resumeHDF5::
packToFile(const string          & mapName,
                 EpetraExt::HDF5 & hdf5,
           const pMap<ITEM>      & morganaMap)
{
  assert(commDevLoaded);
  
  EpetraExt::DistArray<Real> distArray = convertToDistArray(morganaMap);
  
  commDev->barrier();
  hdf5.Write(mapName,distArray);
}    

template<typename ITEM> 
void
resumeHDF5::
printToFile(const string & s,
            const pMap<ITEM> & morganaMap)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);  
  EpetraExt::HDF5 hdf5(epetraComm);
  
  hdf5.Create(s + ".h5");
  packToFile("map",
             hdf5,
             morganaMap);
  
  hdf5.Close();
}

template<typename DATA, typename MAP>
void
resumeHDF5::
packToFile(const string          & vectName_map,
           const string          & vectName_data,
                 EpetraExt::HDF5 & hdf5,
           const pVect<DATA,MAP> & morganaVect)
{
  assert(commDevLoaded);
  
  EpetraExt::DistArray<Real>  mapArray = convertToDistArray(morganaVect.getMapRef());
  EpetraExt::DistArray<Real> dataArray = convertToDistArray(morganaVect);
  
  hdf5.Write(vectName_map, mapArray);
  hdf5.Write(vectName_data,dataArray);
}

template<typename DATA, typename MAP>
void
resumeHDF5::
printToFile(const string & s,
            const pVect<DATA,MAP> & morganaVect)
{
  assert(commDevLoaded);
  
  EpetraExt::DistArray<Real>  mapArray = convertToDistArray(morganaVect.getMapRef());
  EpetraExt::DistArray<Real> dataArray = convertToDistArray(morganaVect);
  
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  
  commDev->barrier();
  
  hdf5.Create(s + ".h5");
  hdf5.Write("map", mapArray);
  hdf5.Write("data",dataArray);
  hdf5.Close();
  
  commDev->barrier();
}
    
template<typename ITEM, typename ROWMAP, typename COLMAP>
void
resumeHDF5::
packToFile(const string          & graphName_rowMap,
           const string          & graphName_colMap,
           const string          & graphName_data,
                 EpetraExt::HDF5 & hdf5,
           const pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph)
{
  EpetraExt::DistArray<Real>  rowArray = convertToDistArray(morganaGraph.getRowMap());
  EpetraExt::DistArray<Real>  colArray = convertToDistArray(morganaGraph.getColMap());
  EpetraExt::DistArray<Real> dataArray = convertToDistArray(morganaGraph);
  
  commDev->barrier();
  hdf5.Write(graphName_rowMap, rowArray);
  
  commDev->barrier();
  hdf5.Write(graphName_colMap, colArray);
  
  commDev->barrier();
  hdf5.Write(graphName_data,  dataArray);
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
resumeHDF5::
printToFile(const string & s,
            const pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  
  hdf5.Create(s + ".h5");
  packToFile("rowMap",
	     "colMap",
	     "data",
	     hdf5,
	     morganaGraph);
  hdf5.Close();
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - CONTAINERS
//-------------------------------------------------------------------------------------------------
template<typename DATA>
void
resumeHDF5::
extractFromFile(const string          & vectName,
                      EpetraExt::HDF5 & hdf5,
                      sVect<DATA>     & morganaVect)
{
  assert(commDevLoaded);
  
  EpetraExt::DistArray<Real> * dataArray;
  hdf5.Read(vectName, dataArray);
  
  convertFromDistArray(*dataArray,morganaVect);
  delete dataArray;
}

template<typename DATA>
void
resumeHDF5::
loadFromFile(const string & s,
                   sVect<DATA> & morganaVect)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::DistArray<Real> * dataArray;
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Open(s + ".h5");
  hdf5.Read("sVect",dataArray);
  hdf5.Close();
  
  convertFromDistArray(*dataArray,morganaVect);
  delete dataArray;
}

template<typename DATA>
void
resumeHDF5::
extractFromFile(const string              & arrayName_index,
                const string              & arrayName_serialized,
                      EpetraExt::HDF5     & hdf5,
                      sVect<sVect<DATA> > & morganaArray)
{
  assert(commDevLoaded);
  
  sVect<UInt> inIndex;
  sVect<DATA> inSerial;
  
  extractFromFile(arrayName_index,      hdf5, inIndex);
  extractFromFile(arrayName_serialized, hdf5, inSerial);
  
  constructArray(inIndex,inSerial,morganaArray);
}

template<typename DATA>
void
resumeHDF5::
loadFromFile(const string & s,
                   sVect<sVect<DATA> > & morganaArray)
{
  assert(commDevLoaded);
 
  sVect<UInt> inIndex;
  sVect<DATA> inSerial;
  
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Open(s + ".h5");
  
  extractFromFile("index",      hdf5, inIndex);
  extractFromFile("serialized", hdf5, inSerial);
  
  hdf5.Close();
  
  constructArray(inIndex,inSerial,morganaArray);
}

template<typename ITEM> 
void
resumeHDF5::
extractFromFile(const string          & mapName,
                      EpetraExt::HDF5 & hdf5,
                      pMap<ITEM>      & morganaMap)
{
  assert(commDevLoaded);
  
  EpetraExt::DistArray<Real> * distArray;
  hdf5.Read(mapName, distArray);
  
  convertFromDistArray<ITEM>(*distArray,morganaMap);
  delete distArray;
  
  pMapComm<ITEM> mapComm(commDev);
  mapComm.vectorPid(morganaMap);
  
  pMapManip<ITEM> manip;
  manip.setIndexing(morganaMap);
}    

template<typename ITEM> 
void
resumeHDF5::
loadFromFile(const string & s,
                   pMap<ITEM> & morganaMap)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  
  hdf5.Open(s + ".h5");
  extractFromFile("map",hdf5,morganaMap);
  hdf5.Close();
}

template<typename DATA, typename MAP>
void
resumeHDF5::
extractFromFile(const string          & vectName_map,
                const string          & vectName_data,
                      EpetraExt::HDF5 & hdf5,
                      pVect<DATA,MAP> & morganaVect)
{
  assert(commDevLoaded);
  
  EpetraExt::DistArray<Real> *  mapArray;
  EpetraExt::DistArray<Real> * dataArray;
  
  hdf5.Read(vectName_map, mapArray);
  hdf5.Read(vectName_data,dataArray);
  
  pMap<MAP> map;
  convertFromDistArray<MAP>(*mapArray,map);
  convertFromDistArray<DATA,MAP>(*dataArray,morganaVect);
  morganaVect.setMap(map);
  morganaVect.bufferLids();
  morganaVect.updateFinder();
  
  delete mapArray;
  delete dataArray;
  
  pVectComm<DATA,MAP> vectComm(commDev);
  vectComm.vectorPid(morganaVect);
  
  morganaVect.restoreLids();
  pVectManip<DATA,MAP> manip;
  manip.setIndexing(morganaVect);
  
  morganaVect.updateFinder();
}

template<typename DATA, typename MAP>
void
resumeHDF5::
loadFromFile(const string & s,
             pVect<DATA,MAP> & morganaVect)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);  
  EpetraExt::HDF5 hdf5(epetraComm);
  
  hdf5.Open(s + ".h5");
  extractFromFile("map","data",hdf5,morganaVect);
  hdf5.Close();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
resumeHDF5::
extractFromFile(const string          & graphName_rowMap,
                const string          & graphName_colMap,
                const string          & graphName_data,
                      EpetraExt::HDF5 & hdf5,
                      pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph)
{
  assert(commDevLoaded);
  
  //Loading----------------------------------------------------------
  EpetraExt::DistArray<Real> *  rowArray;
  EpetraExt::DistArray<Real> *  colArray;
  EpetraExt::DistArray<Real> * dataArray;
  
  hdf5.Read(graphName_rowMap, rowArray);
  hdf5.Read(graphName_colMap, colArray);
  hdf5.Read(graphName_data,  dataArray);
  
  pMap<ROWMAP> rowMap;
  pMap<COLMAP> colMap;
  convertFromDistArray<ROWMAP>(*rowArray,rowMap);
  convertFromDistArray<COLMAP>(*colArray,colMap);
  convertFromDistArray<ITEM,ROWMAP,COLMAP>(*dataArray,morganaGraph);
    
  delete rowArray;
  delete colArray;
  delete dataArray;
  
  //Post processing rows---------------------------------------------
  morganaGraph.setRowMap(rowMap);
  morganaGraph.bufferLids();
  morganaGraph.updateRowFinder();
  
  pVectComm<ITEM,ROWMAP> rowComm(commDev);
  rowComm.vectorPid(morganaGraph);
  
  morganaGraph.restoreLids();
  pVectManip<ITEM,ROWMAP> manipRow;
  manipRow.setIndexing(morganaGraph);
  
  morganaGraph.updateRowFinder();
  
  //Post processing cols---------------------------------------------
  pMapComm<COLMAP> colComm(commDev);
  colComm.vectorPid(colMap);
  
  pMapManip<COLMAP> colManip;
  colManip.setIndexing(colMap);
  
  morganaGraph.setColMap(colMap);
  morganaGraph.updateColFinder();
}

template<typename ITEM, typename ROWMAP, typename COLMAP>
void
resumeHDF5::
loadFromFile(const string & s,
                   pGraph<ITEM,ROWMAP,COLMAP> & morganaGraph)
{
  assert(commDevLoaded);
  
  Epetra_MpiComm epetraComm(*commDev);  
  EpetraExt::HDF5 hdf5(epetraComm);
  
  hdf5.Open(s + ".h5");  
  extractFromFile("rowMap","colMap","data",hdf5,morganaGraph);  
  hdf5.Close();
}


//_________________________________________________________________________________________________
// PRINT TO FILE - GEOMETRY
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
packToFile(const string & meshName,
                 EpetraExt::HDF5 & hdf5,
           const mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Counts
  mesh1dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> manip(commDev);
  UInt numVertices = manip.getNumGlobalVertices(morganaMesh);
  UInt numElements = manip.getNumGlobalElements(morganaMesh);
  
  //Saving params
  hdf5.Write(meshName + "flags", "standard",    Real(morganaMesh.standard));
  hdf5.Write(meshName + "flags", "numVertices", Real(numVertices));
  hdf5.Write(meshName + "flags", "numElements", Real(numElements));
  
  //Saving data
  if(numVertices != 0) { packToFile("nodes_map",       "nodes_data", hdf5, morganaMesh.nodes);    }
  if(numElements != 0) { packToFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, morganaMesh.elements); }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
printToFile(const string & s,
            const mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Counts
  mesh1dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> manip(commDev);
  UInt numVertices = manip.getNumGlobalVertices(morganaMesh);
  UInt numElements = manip.getNumGlobalElements(morganaMesh);
  
  //Open
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  //Saving params
  hdf5.Write("flags", "standard",    Real(morganaMesh.standard));
  hdf5.Write("flags", "numVertices", Real(numVertices));
  hdf5.Write("flags", "numElements", Real(numElements));
  
  //Saving data
  if(numVertices != 0) { packToFile("nodes_map",       "nodes_data", hdf5, morganaMesh.nodes);    }
  if(numElements != 0) { packToFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, morganaMesh.elements); }
  
  //Close
  hdf5.Close();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
packToFile(const string          & meshName,
                 EpetraExt::HDF5 & hdf5,
           const mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Counts
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> manip(commDev);
  UInt numVertices = manip.getNumGlobalVertices(morganaMesh);
  UInt numEdges    = manip.getNumGlobalEdges(morganaMesh);
  UInt numElements = manip.getNumGlobalElements(morganaMesh);
  
  //Saving params
  hdf5.Write(meshName + "flags", "standard",    Real(morganaMesh.standard));
  hdf5.Write(meshName + "flags", "numVertices", Real(numVertices));
  hdf5.Write(meshName + "flags", "numEdges",    Real(numEdges));
  hdf5.Write(meshName + "flags", "numElements", Real(numElements));
  
  //Saving data
  if(numVertices != 0) { packToFile("nodes_map",       "nodes_data", hdf5, morganaMesh.nodes);    }
  if(numEdges    != 0) { packToFile("edges_rowMap",    "edges_colMap",    "edges_data",    hdf5, morganaMesh.edges);    }
  if(numElements != 0) { packToFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, morganaMesh.elements); }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
printToFile(const string & s,
            const mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Counts
  mesh2dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> manip(commDev);
  UInt numVertices = manip.getNumGlobalVertices(morganaMesh);
  UInt numEdges    = manip.getNumGlobalEdges(morganaMesh);
  UInt numElements = manip.getNumGlobalElements(morganaMesh);
  
  //Open
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  //Saving params
  hdf5.Write("flags", "standard",    Real(morganaMesh.standard));
  hdf5.Write("flags", "numVertices", Real(numVertices));
  hdf5.Write("flags", "numEdges",    Real(numEdges));
  hdf5.Write("flags", "numElements", Real(numElements));
  
  //Saving data
  if(numVertices != 0) { packToFile("nodes_map",       "nodes_data", hdf5, morganaMesh.nodes);    }
  if(numEdges    != 0) { packToFile("edges_rowMap",    "edges_colMap",    "edges_data",    hdf5, morganaMesh.edges);    }
  if(numElements != 0) { packToFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, morganaMesh.elements); }
  
  //Close
  hdf5.Close();
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
packToFile(const string          & meshName,
                 EpetraExt::HDF5 & hdf5,
           const mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Counts
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> manip(commDev);
  UInt numVertices = manip.getNumGlobalVertices(morganaMesh);
  UInt numEdges    = manip.getNumGlobalEdges(morganaMesh);
  UInt numFaces    = manip.getNumGlobalFaces(morganaMesh);
  UInt numElements = manip.getNumGlobalElements(morganaMesh);
  
  //Saving params
  hdf5.Write(meshName + "flags", "standard",    Real(morganaMesh.standard));
  hdf5.Write(meshName + "flags", "numVertices", Real(numVertices));
  hdf5.Write(meshName + "flags", "numEdges",    Real(numEdges));
  hdf5.Write(meshName + "flags", "numFaces",    Real(numFaces));
  hdf5.Write(meshName + "flags", "numElements", Real(numElements));
  
  //Saving data
  if(numVertices != 0) { packToFile("nodes_map",       "nodes_data", hdf5, morganaMesh.nodes);    }
  if(numEdges    != 0) { packToFile("edges_rowMap",    "edges_colMap",    "edges_data",    hdf5, morganaMesh.edges);    }
  if(numFaces    != 0) { packToFile("faces_rowMap",    "faces_colMap",    "faces_data",    hdf5, morganaMesh.faces);    }
  if(numElements != 0) { packToFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, morganaMesh.elements); }
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
printToFile(const string & s,
            const mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Counts
  mesh3dGlobalManip<GEOSHAPE,ELMAP,NODEMAP> manip(commDev);
  UInt numVertices = manip.getNumGlobalVertices(morganaMesh);
  UInt numEdges    = manip.getNumGlobalEdges(morganaMesh);
  UInt numFaces    = manip.getNumGlobalFaces(morganaMesh);
  UInt numElements = manip.getNumGlobalElements(morganaMesh);
  
  //Open
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  //Saving params
  hdf5.Write("flags", "standard",    Real(morganaMesh.standard));
  hdf5.Write("flags", "numVertices", Real(numVertices));
  hdf5.Write("flags", "numEdges",    Real(numEdges));
  hdf5.Write("flags", "numFaces",    Real(numFaces));
  hdf5.Write("flags", "numElements", Real(numElements));
  
  //Saving data
  if(numVertices != 0) { packToFile("nodes_map",       "nodes_data", hdf5, morganaMesh.nodes);    }
  if(numEdges    != 0) { packToFile("edges_rowMap",    "edges_colMap",    "edges_data",    hdf5, morganaMesh.edges);    }
  if(numFaces    != 0) { packToFile("faces_rowMap",    "faces_colMap",    "faces_data",    hdf5, morganaMesh.faces);    }
  if(numElements != 0) { packToFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, morganaMesh.elements); }
  
  //Close
  hdf5.Close();
}


//_________________________________________________________________________________________________
// LOAD FROM FILE - GEOMETRY
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
extractFromFile(const string          & meshName,
                      EpetraExt::HDF5 & hdf5,
                      mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Tyepdefs
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT   BOOLVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D    GRAPH1D;
  
  //Loading params
  Real data, numVertices, numElements;
  hdf5.Read(meshName + "flags", "standard",    data);  morganaMesh.standard = MeshStandards(data);
  hdf5.Read(meshName + "flags", "numVertices", numVertices); 
  hdf5.Read(meshName + "flags", "numElements", numElements);
  
  //Loading nodes
  if(numVertices != 0)
  {
    NODESVECT nodes;
    extractFromFile("nodes_map", "nodes_data", hdf5, nodes);
    morganaMesh.nodes = nodes;
  }

  //Loading elements
  if(numElements != 0)
  {
    GRAPH1D elements;
    extractFromFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, elements);
    morganaMesh.elements = elements;
  }
  
  //Update mesh
  morganaMesh.computeNumVertices();
  morganaMesh.transferMap();
  morganaMesh.verticesComputed = true;
  morganaMesh.mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
loadFromFile(const string & s,
                   mesh1d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Tyepdefs
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT   BOOLVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D    GRAPH1D;
  
  //Open
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Open(s + ".h5");
  
  //Loading params
  Real data, numVertices, numElements;
  hdf5.Read("flags", "standard",    data);  morganaMesh.standard = MeshStandards(data);
  hdf5.Read("flags", "numVertices", numVertices); 
  hdf5.Read("flags", "numElements", numElements);
  
  //Loading nodes
  if(numVertices != 0)
  {
    NODESVECT nodes;
    extractFromFile("nodes_map", "nodes_data", hdf5, nodes);
    morganaMesh.nodes = nodes;
  }

  //Loading elements
  if(numElements != 0)
  {
    GRAPH1D elements;
    extractFromFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, elements);
    morganaMesh.elements = elements;
  }

  //Close
  hdf5.Close();
  
  //Update mesh
  morganaMesh.computeNumVertices();
  morganaMesh.transferMap();
  morganaMesh.verticesComputed = true;
  morganaMesh.mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
extractFromFile(const string          & meshName,
                      EpetraExt::HDF5 & hdf5,
                      mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Tyepdefs
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT   BOOLVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH2D    GRAPH2D;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D    GRAPH1D;
  
  //Loading params
  Real data, numVertices, numEdges, numElements;
  hdf5.Read(meshName + "flags", "standard",    data);  morganaMesh.standard = MeshStandards(data);
  hdf5.Read(meshName + "flags", "numVertices", numVertices); 
  hdf5.Read(meshName + "flags", "numEdges",    numEdges);
  hdf5.Read(meshName + "flags", "numElements", numElements);
  
  //Loading nodes
  if(numVertices != 0)
  {
    NODESVECT nodes;
    extractFromFile("nodes_map", "nodes_data", hdf5, nodes);
    morganaMesh.nodes = nodes;
  }
  
  //Loading edges
  if(numEdges != 0)
  {
    GRAPH1D edges;
    extractFromFile("edges_rowMap", "edges_colMap", "edges_data", hdf5, edges);
    morganaMesh.edges = edges;
  }

  //Loading elements
  if(numElements != 0)
  {
    GRAPH2D elements;
    extractFromFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, elements);
    morganaMesh.elements = elements;
  }
  
  //Update mesh
  morganaMesh.computeNumVertices();
  morganaMesh.transferMap();
  morganaMesh.verticesComputed = true;
  morganaMesh.mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
loadFromFile(const string & s,
                   mesh2d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Tyepdefs
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT   BOOLVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH2D    GRAPH2D;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D    GRAPH1D;
  
  //Open
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Open(s + ".h5");
  
  //Loading params
  Real data, numVertices, numEdges, numElements;
  hdf5.Read("flags", "standard",    data);  morganaMesh.standard = MeshStandards(data);
  hdf5.Read("flags", "numVertices", numVertices); 
  hdf5.Read("flags", "numEdges",    numEdges);
  hdf5.Read("flags", "numElements", numElements);
  
  //Loading nodes
  if(numVertices != 0)
  {
    NODESVECT nodes;
    extractFromFile("nodes_map", "nodes_data", hdf5, nodes);
    morganaMesh.nodes = nodes;
  }
  
  //Loading edges
  if(numEdges != 0)
  {
    GRAPH1D edges;
    extractFromFile("edges_rowMap", "edges_colMap", "edges_data", hdf5, edges);
    morganaMesh.edges = edges;
  }

  //Loading elements
  if(numElements != 0)
  {
    GRAPH2D elements;
    extractFromFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, elements);
    morganaMesh.elements = elements;
  }

  //Close
  hdf5.Close();
  
  //Update mesh
  morganaMesh.computeNumVertices();
  morganaMesh.transferMap();
  morganaMesh.verticesComputed = true;
  morganaMesh.mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
extractFromFile(const string          & meshName,
                      EpetraExt::HDF5 & hdf5,
                      mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Tyepdefs
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT   BOOLVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH3D    GRAPH3D;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH2D    GRAPH2D;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D    GRAPH1D;
  
  //Loading params
  Real data, numVertices, numEdges, numFaces, numElements;
  hdf5.Read(meshName + "flags", "standard",    data);  morganaMesh.standard = MeshStandards(data);
  hdf5.Read(meshName + "flags", "numVertices", numVertices); 
  hdf5.Read(meshName + "flags", "numEdges",    numEdges);
  hdf5.Read(meshName + "flags", "numFaces",    numFaces);
  hdf5.Read(meshName + "flags", "numElements", numElements);
  
  //Loading nodes
  if(numVertices != 0)
  {
    NODESVECT nodes;
    extractFromFile("nodes_map", "nodes_data", hdf5, nodes);
    morganaMesh.nodes = nodes;
  }
  
  //Loading edges
  if(numEdges != 0)
  {
    GRAPH1D edges;
    extractFromFile("edges_rowMap", "edges_colMap", "edges_data", hdf5, edges);
    morganaMesh.edges = edges;
  }
  
  //Loading faces
  if(numFaces != 0)
  {
    GRAPH2D faces;
    extractFromFile("faces_rowMap", "faces_colMap", "faces_data", hdf5, faces);
    morganaMesh.faces = faces;
  }

  //Loading elements
  if(numElements != 0)
  {
    GRAPH3D elements;
    extractFromFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, elements);
    morganaMesh.elements = elements;
  }
  
  //Update mesh
  morganaMesh.computeNumVertices();
  morganaMesh.transferMap();
  morganaMesh.verticesComputed = true;
  morganaMesh.mapTransferred   = true;
}

template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
void
resumeHDF5::
loadFromFile(const string & s,
                   mesh3d<GEOSHAPE,ELMAP,NODEMAP> & morganaMesh)
{
  assert(commDevLoaded);
  
  //Tyepdefs
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::NODESVECT  NODESVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::BOOLVECT   BOOLVECT;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH3D    GRAPH3D;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH2D    GRAPH2D;
  typedef typename mesh3d<GEOSHAPE,ELMAP,NODEMAP>::GRAPH1D    GRAPH1D;
  
  //Open
  Epetra_MpiComm epetraComm(*commDev);
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Open(s + ".h5");
  
  //Loading params
  Real data, numVertices, numEdges, numFaces, numElements;
  hdf5.Read("flags", "standard",    data);  morganaMesh.standard = MeshStandards(data);
  hdf5.Read("flags", "numVertices", numVertices); 
  hdf5.Read("flags", "numEdges",    numEdges);
  hdf5.Read("flags", "numFaces",    numFaces);
  hdf5.Read("flags", "numElements", numElements);
  
  //Loading nodes
  if(numVertices != 0)
  {
    NODESVECT nodes;
    extractFromFile("nodes_map", "nodes_data", hdf5, nodes);
    morganaMesh.nodes = nodes;
  }
  
  //Loading edges
  if(numEdges != 0)
  {
    GRAPH1D edges;
    extractFromFile("edges_rowMap", "edges_colMap", "edges_data", hdf5, edges);
    morganaMesh.edges = edges;
  }
  
  //Loading faces
  if(numFaces != 0)
  {
    GRAPH2D faces;
    extractFromFile("faces_rowMap", "faces_colMap", "faces_data", hdf5, faces);
    morganaMesh.faces = faces;
  }

  //Loading elements
  if(numElements != 0)
  {
    GRAPH3D elements;
    extractFromFile("elements_rowMap", "elements_colMap", "elements_data", hdf5, elements);
    morganaMesh.elements = elements;
  }

  //Close
  hdf5.Close();
  
  //Update mesh
  morganaMesh.computeNumVertices();
  morganaMesh.transferMap();
  morganaMesh.verticesComputed = true;
  morganaMesh.mapTransferred   = true;
}


#endif
