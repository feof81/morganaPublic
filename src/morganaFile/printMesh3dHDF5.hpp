/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef PRINTMESH3DHDF5_HPP
#define PRINTMESH3DHDF5_HPP

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "EpetraExt_HDF5.h"

#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"

#include "pMapGlobalManip.h"
#include "pVectGlobalManip.hpp"
#include "mesh3dGlobalManip.hpp"
#include "printMesh.hpp"

#include "traitsDofType.h"
#include "traitsGeoType.h"

using namespace std;


/*! Print to HDF5 a mesh3d */
template<typename GEOSHAPE, typename PMAPTYPE>
class printMesh3dHDF5
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh3d<GEOSHAPE,PMAPTYPE,PMAPTYPE> MESH3D;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructor and functions */ //@{
  public:
    printMesh3dHDF5();
    printMesh3dHDF5(const Teuchos::RCP<communicator> & CommDev);
    printMesh3dHDF5(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    EpetraExt::DistArray<int>    convertElements(const MESH3D & grid3d) const;
    EpetraExt::DistArray<double> convertNodes(const MESH3D & grid3d) const;
    void print(const string & s, const UInt & pid, const Teuchos::RCP<MESH3D> & grid3d) const;
    void print(const string & s, const UInt & pid, const MESH3D & grid3d) const;
    void printGeoIds(const string & s, const UInt & pid, const Teuchos::RCP<MESH3D> & grid3d) const;
    void printGeoIds(const string & s, const UInt & pid, const MESH3D & grid3d) const;
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
printMesh3dHDF5()
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  commDevLoaded = false;
}

template<typename GEOSHAPE, typename PMAPTYPE>
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
printMesh3dHDF5(const Teuchos::RCP<communicator> & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename PMAPTYPE>
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
printMesh3dHDF5(communicator & CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev       = CommDev;
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev       = Teuchos::rcpFromRef(CommDev);
}

template<typename GEOSHAPE, typename PMAPTYPE>
EpetraExt::DistArray<int>
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
convertElements(const MESH3D & grid3d) const
{
  assert(commDevLoaded);
  
  //Alloc------------------------------------------------------------------------------------------
  Epetra_MpiComm            epetraComm(*commDev);
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  
  
  //Elements conversion----------------------------------------------------------------------------
  typedef pMap<PMAPTYPE>                PMAP;
  typedef typename MESH3D::GRAPH3D      ELEMENTS;
  typedef typename MESH3D::GEOELEMENT3D GEOELEMENT3D;
  
  ELEMENTS elements_morgana = grid3d.getElements();                 //Extract the elements
  elements_morgana.pushToGlobal();
  
  pVectGlobalManip<GEOELEMENT3D,PMAPTYPE> elementsManip(commDev);   //Redistribute
  Epetra_Map elementsMap_epetra = elementsManip.vectorLinear(elements_morgana);
  
  EpetraExt::DistArray<int> elements_epetra(elementsMap_epetra,GEOSHAPE::numPoints);  //Create the elements for epetra
  
  for(UInt i=1; i <= elements_morgana.size(); ++i)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      elements_epetra(i-1,j-1) = elements_morgana.getItemL(i).getCid(j) - 1;
    }
  }
  
  //Linear map-------------------------------------------------------------------------------------
  return(elements_epetra);
}

template<typename GEOSHAPE, typename PMAPTYPE>
EpetraExt::DistArray<double>
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
convertNodes(const MESH3D & grid3d) const
{
  assert(commDevLoaded);
  
  //Alloc------------------------------------------------------------------------------------------
  Epetra_MpiComm            epetraComm(*commDev);
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  
  
  //Nodes conversion-------------------------------------------------------------------------------
  typedef typename MESH3D::NODESVECT NODES;
  
  NODES nodes_morgana = grid3d.getNodes();  //Remap the nodes
  
  pVectGlobalManip<point3d,PMAPTYPE> nodesManip(commDev);
  Epetra_Map nodesMap_epetra = nodesManip.vectorLinear(nodes_morgana);
  
  EpetraExt::DistArray<double> nodes_epetra(nodesMap_epetra,3);
  
  for(UInt i=1; i <= nodes_morgana.size(); ++i)
  {
    nodes_epetra(i-1,0) = nodes_morgana(i).getX();
    nodes_epetra(i-1,1) = nodes_morgana(i).getY();
    nodes_epetra(i-1,2) = nodes_morgana(i).getZ();
  }
  
  //Linear map-------------------------------------------------------------------------------------
  return(nodes_epetra);
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
print(const string & s, const UInt & pid, const Teuchos::RCP<MESH3D> & grid3d) const
{
  assert(commDevLoaded);
  
  //Grid test--------------------------------------------------------------------------------------
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid3d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid3d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid3d);
  
  //Alloc------------------------------------------------------------------------------------------ 
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = convertElements(*grid3d);
  EpetraExt::DistArray<double> nodes_epetra    = convertNodes(*grid3d);
  
  //HDF5 print-------------------------------------------------------------------------------------
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  hdf5.Write("elements",elements_epetra);
  hdf5.Write("nodes",nodes_epetra);
  hdf5.Close();

  
  //Print XDMF-------------------------------------------------------------------------------------
  if(pid == 0)
  {
    string fileName = s + ".xdmf";
    ofstream xdfmFile(fileName.c_str());
    
    xdfmFile << "<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">"  << endl;
    xdfmFile << "  <Domain>"                                                             << endl;
    xdfmFile << "    <Grid Name=\"grid\">"                                               << endl;
    xdfmFile << "      <Topology ";  geoTrait.hdft5String(xdfmFile); xdfmFile << " NumberOfElements=\"" << numElementsG << "\">"         << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << numElementsG << " " << GEOSHAPE::numPoints << "\">" << endl;
    xdfmFile << "          " << s << ".h5:/elements/Values"                              << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Topology>"                                                      << endl;
    xdfmFile << "      <Geometry GeometryType=\"XYZ\">"                                  << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numNodesG << " 3\">" << endl;
    xdfmFile << "          " << s << ".h5:/nodes/Values"                                 << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Geometry>"                                                      << endl;    
    xdfmFile << "    </Grid>"                                                            << endl;
    xdfmFile << "  </Domain>"                                                            << endl;
    xdfmFile << "</Xdmf>"                                                                << endl;  
    
    xdfmFile.close();
  }
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
print(const string & s, const UInt & pid, const MESH3D & grid3d) const
{
  assert(commDevLoaded);
  
  //Grid test--------------------------------------------------------------------------------------
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid3d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid3d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid3d);
  
  //Alloc------------------------------------------------------------------------------------------ 
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = convertElements(grid3d);
  EpetraExt::DistArray<double> nodes_epetra    = convertNodes(grid3d);
  
  //HDF5 print-------------------------------------------------------------------------------------
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  hdf5.Write("elements",elements_epetra);
  hdf5.Write("nodes",nodes_epetra);
  hdf5.Close();

  
  //Print XDMF-------------------------------------------------------------------------------------
  if(pid == 0)
  {
    string fileName = s + ".xdmf";
    ofstream xdfmFile(fileName.c_str());
    
    xdfmFile << "<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">"  << endl;
    xdfmFile << "  <Domain>"                                                             << endl;
    xdfmFile << "    <Grid Name=\"grid\">"                                               << endl;
    xdfmFile << "      <Topology ";  geoTrait.hdft5String(xdfmFile); xdfmFile << " NumberOfElements=\"" << numElementsG << "\">"         << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << numElementsG << " " << GEOSHAPE::numPoints << "\">" << endl;
    xdfmFile << "          " << s << ".h5:/elements/Values"                              << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Topology>"                                                      << endl;
    xdfmFile << "      <Geometry GeometryType=\"XYZ\">"                                  << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numNodesG << " 3\">" << endl;
    xdfmFile << "          " << s << ".h5:/nodes/Values"                                 << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Geometry>"                                                      << endl;    
    xdfmFile << "    </Grid>"                                                            << endl;
    xdfmFile << "  </Domain>"                                                            << endl;
    xdfmFile << "</Xdmf>"                                                                << endl;  
    
    xdfmFile.close();
  }
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
printGeoIds(const string & s, const UInt & pid, const Teuchos::RCP<MESH3D> & grid3d) const
{
  assert(commDevLoaded);
  typedef Real DOFTYPE;
  
  //Grid test--------------------------------------------------------------------------------------
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid3d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid3d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid3d);
  
  //Alloc------------------------------------------------------------------------------------------ 
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = convertElements(*grid3d);
  EpetraExt::DistArray<double> nodes_epetra    = convertNodes(*grid3d);
   
  //GeoIds creation--------------------------------------------------------------------------------
  typedef Real  OUTTYPE;
  typedef pVect<OUTTYPE,PMAPTYPE> OUTVECT; 
   
  sVect<OUTTYPE> tempVect(grid3d->getNumElements());
   
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  { tempVect(i) = Real(grid3d->getElementL(i).getGeoId()); }
   
  OUTVECT data_morgana(grid3d->getElements().getMapRef(), tempVect);
  
  pVectGlobalManip<OUTTYPE,PMAPTYPE> dataManip(commDev); //Remap the vector
  Epetra_Map geoIdMap_epetra = dataManip.vectorLinear(data_morgana);  
   
  //HDF5 print-------------------------------------------------------------------------------------
  traitsDofType<DOFTYPE> dofTrait;
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  hdf5.Write("elements",elements_epetra);                   //Print the mesh
  hdf5.Write("nodes",nodes_epetra);
  
  dofTrait.printHDF5(hdf5, data_morgana, geoIdMap_epetra); //Print the fields
  
  hdf5.Close();

  //Print XDMF-------------------------------------------------------------------------------------
  if(pid == 0)
  {
    string fileName = s + ".xdmf";
    ofstream xdfmFile(fileName.c_str());
    
    xdfmFile << "<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">"  << endl;
    xdfmFile << "  <Domain>"                                                             << endl;
    xdfmFile << "    <Grid Name=\"grid\">"                                               << endl;
    xdfmFile << "      <Topology ";  geoTrait.hdft5String(xdfmFile); xdfmFile << " NumberOfElements=\"" << numElementsG << "\">"         << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << numElementsG << " " << GEOSHAPE::numPoints << "\">" << endl;
    xdfmFile << "          " << s << ".h5:/elements/Values"                              << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Topology>"                                                      << endl;
    xdfmFile << "      <Geometry GeometryType=\"XYZ\">"                                  << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numNodesG << " 3\">" << endl;
    xdfmFile << "          " << s << ".h5:/nodes/Values"                                 << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Geometry>"                                                      << endl;
    
    dofTrait.printXDMF(xdfmFile,"Cell",s,numElementsG);
    
    xdfmFile << "    </Grid>"                                                            << endl;
    xdfmFile << "  </Domain>"                                                            << endl;
    xdfmFile << "</Xdmf>"                                                                << endl;  
    
    xdfmFile.close();
  }
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
printMesh3dHDF5<GEOSHAPE,PMAPTYPE>::
printGeoIds(const string & s, const UInt & pid, const MESH3D & grid3d) const
{
  assert(commDevLoaded);
  typedef Real DOFTYPE; 
  
  //Grid test--------------------------------------------------------------------------------------
  mesh3dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid3d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid3d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid3d);
  
  //Alloc------------------------------------------------------------------------------------------ 
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = convertElements(grid3d);
  EpetraExt::DistArray<double> nodes_epetra    = convertNodes(grid3d);
   
  //GeoIds creation--------------------------------------------------------------------------------
  typedef Real  OUTTYPE;
  typedef pVect<OUTTYPE,PMAPTYPE> OUTVECT; 
   
  sVect<OUTTYPE> tempVect(grid3d.getNumElements());
   
  for(UInt i=1; i <= grid3d.getNumElements(); ++i)
  { tempVect(i) = Real(grid3d.getElementL(i).getGeoId()); }
   
  OUTVECT data_morgana(grid3d.getElements().getMapRef(), tempVect);
  
  pVectGlobalManip<OUTTYPE,PMAPTYPE> dataManip(commDev); //Remap the vector
  Epetra_Map geoIdMap_epetra = dataManip.vectorLinear(data_morgana);  
   
  //HDF5 print-------------------------------------------------------------------------------------
  traitsDofType<DOFTYPE> dofTrait;
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  hdf5.Write("elements",elements_epetra);                     //Print the mesh
  hdf5.Write("nodes",nodes_epetra);
  
  dofTrait.printHDF5(hdf5, data_morgana, geoIdMap_epetra); //Print the fields
  
  hdf5.Close();
  
  //Print XDMF-------------------------------------------------------------------------------------
  if(pid == 0)
  {
    string fileName = s + ".xdmf";
    ofstream xdfmFile(fileName.c_str());
    
    xdfmFile << "<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">"  << endl;
    xdfmFile << "  <Domain>"                                                             << endl;
    xdfmFile << "    <Grid Name=\"grid\">"                                               << endl;
    xdfmFile << "      <Topology ";  geoTrait.hdft5String(xdfmFile); xdfmFile << " NumberOfElements=\"" << numElementsG << "\">"         << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << numElementsG << " " << GEOSHAPE::numPoints << "\">" << endl;
    xdfmFile << "          " << s << ".h5:/elements/Values"                              << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Topology>"                                                      << endl;
    xdfmFile << "      <Geometry GeometryType=\"XYZ\">"                                  << endl;
    xdfmFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numNodesG << " 3\">" << endl;
    xdfmFile << "          " << s << ".h5:/nodes/Values"                                 << endl;
    xdfmFile << "        </DataItem>"                                                    << endl;
    xdfmFile << "      </Geometry>"                                                      << endl;
    
    dofTrait.printXDMF(xdfmFile,"Cell",s,numElementsG);
    
    xdfmFile << "    </Grid>"                                                            << endl;
    xdfmFile << "  </Domain>"                                                            << endl;
    xdfmFile << "</Xdmf>"                                                                << endl;  
    
    xdfmFile.close();
  }
}


#endif
