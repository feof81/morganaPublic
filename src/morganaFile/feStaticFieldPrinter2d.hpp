/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FESTATICFIELDPRINTER2D_HPP
#define FESTATICFIELDPRINTER2D_HPP

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "EpetraExt_HDF5.h"

#include "mesh1d.hpp"
#include "mesh2d.hpp"

#include "printMesh.hpp"
#include "printMesh2dHDF5.hpp"

#include "traitsDofType.h"
#include "traitsGeoType.h"
#include "linearMap_distArray.hpp"

using namespace std;


/*! Print class only for the static finite elements. */
template<typename FIELD>
class feStaticFieldPrinter2d : public printMesh2dHDF5<typename FIELD::GEOSHAPE, typename FIELD::PMAPTYPE>
{ 
    /*! @name Typedefs */ //@{
  public:
    typedef typename FIELD::MESH2D             MESH2D;
    typedef typename FIELD::FIELD_FETYPE       FETYPE;
    typedef typename FIELD::FIELD_DOFTYPE      DOFTYPE;
    typedef typename FIELD::GEOSHAPE           GEOSHAPE;
    typedef typename FIELD::OUTTYPE            OUTTYPE;
    typedef typename FETYPE::FETYPE_PMAPTYPE   PMAPTYPE;
    typedef printMesh2dHDF5<GEOSHAPE,PMAPTYPE> MESHPRINTER;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    static const Real paraFix;
    //@}
    
    /*! @name Constructor and functions */ //@{
  public:
    feStaticFieldPrinter2d();
    feStaticFieldPrinter2d(const Teuchos::RCP<communicator> & CommDev);
    feStaticFieldPrinter2d(communicator & CommDev);
    sVect<OUTTYPE> nodesEval(FIELD & Field) const;
    sVect<OUTTYPE> elementsEval(FIELD & Field) const;
    //@}
    
    /*! @name Print to UCD */ //@{
  public:
    void printUCD_nodes(const string & s, const string & value, const string & unit, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const;
    void printUCD_nodes(const string & s, const string & value, const string & unit, const UInt & pid, FIELD & Field) const;
    void printUCD_elements(const string & s, const string & value, const string & unit, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const;
    void printUCD_elements(const string & s, const string & value, const string & unit, const UInt & pid, FIELD & Field) const;
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    void printHDF5_mesh(const string & s, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const;
    void printHDF5_mesh(const string & s, const UInt & pid, FIELD & Field) const;
    void printHDF5_nodes(const string & s, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const;
    void printHDF5_nodes(const string & s, const UInt & pid, FIELD & Field) const;
    void printHDF5_elements(const string & s, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const;
    void printHDF5_elements(const string & s, const UInt & pid, FIELD & Field) const;
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTOR AND FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FIELD>
const Real feStaticFieldPrinter2d<FIELD>::paraFix = 1.00000000000000001;

template<typename FIELD>
feStaticFieldPrinter2d<FIELD>::
feStaticFieldPrinter2d() : printMesh2dHDF5<typename FIELD::GEOSHAPE, typename FIELD::PMAPTYPE>()
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FIELD>
feStaticFieldPrinter2d<FIELD>::
feStaticFieldPrinter2d(const Teuchos::RCP<communicator> & CommDev) : printMesh2dHDF5<typename FIELD::GEOSHAPE, typename FIELD::PMAPTYPE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FIELD>
feStaticFieldPrinter2d<FIELD>::
feStaticFieldPrinter2d(communicator & CommDev) : printMesh2dHDF5<typename FIELD::GEOSHAPE, typename FIELD::PMAPTYPE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 2>::returnValue);
}

template<typename FIELD>
sVect<typename FIELD::OUTTYPE>
feStaticFieldPrinter2d<FIELD>::
nodesEval(FIELD & Field) const
{ 
  //Alloc
  UInt id;
  OUTTYPE V;
  point3d Y;
  Teuchos::RCP<MESH2D> grid2d = Field.getMesh2d();
  sVect<OUTTYPE> outVect(grid2d->getNumVertices());
  
  for(UInt el=1; el <= grid2d->getNumElements(); el++)
  {
    for(UInt j=1; j <= GEOSHAPE::numPoints; ++j)
    {
      //Evaluate
      Y = GEOSHAPE::getRefNodes(j);
      Field.evalL(el,Y,V);
      
      id = grid2d->getElementL(el).getCid(j);
      outVect(id) = V;
    }
  }
  
  return(outVect);
}

template<typename FIELD>
sVect<typename feStaticFieldPrinter2d<FIELD>::OUTTYPE>
feStaticFieldPrinter2d<FIELD>::
elementsEval(FIELD & Field) const
{
  //Alloc
  OUTTYPE V;
  point3d Y = GEOSHAPE::getBarycenter();
  Teuchos::RCP<MESH2D> grid2d = Field.getMesh2d();
  sVect<OUTTYPE> outVect(grid2d->getNumElements());
  
  for(UInt el=1; el <= grid2d->getNumElements(); el++)
  {    
    Field.evalL(el,Y,V);
    outVect(el) = V;
  }
  
  return(outVect);
}


//_________________________________________________________________________________________________
// PRINT TO UCD
//-------------------------------------------------------------------------------------------------
template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printUCD_nodes(const string & s, const string & value, const string & unit, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const
{
  //Open file
  string fileName = "sf" + num2str(pid) + s + ".inp";
  ofstream out(fileName.c_str());
  
  //Alloc
  Teuchos::RCP<MESH2D> grid2d = Field->getMesh2d();
  traitsDofType<DOFTYPE>  dofTrait;
  traitsGeoType<GEOSHAPE> geoTrait;
  
  //Print general data
  out << grid2d->getNumVertices() << " " << grid2d->getNumElements() << " 1 0 0" << endl;
  
  //Print mesh nodes
  for(UInt i=1; i<= grid2d->getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * grid2d->getNodeL(i).getX() << " ";
    out << paraFix * grid2d->getNodeL(i).getY() << " ";
    out << paraFix * grid2d->getNodeL(i).getZ() << endl;
  }
  
  //Print element nodes
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    out << i << " ";
    out << grid2d->getElementL(i).getGeoId() << "  ";
    out << geoTrait.parariewString() << " ";
    
    for(UInt k=1; k <= grid2d->getElementL(i).size(); ++k)
    {
      out << grid2d->getElementL(i).getCid(k) << " ";
    }
    out << endl;
  }
  
  //Compute the dofs vector
  sVect<OUTTYPE> outVect = nodesEval(*Field);
  
  //Print dofs
  dofTrait.paraviewString(out,value,unit);
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    out << i << " ";
    dofTrait.printDof(out,outVect(i));
  }
}

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printUCD_nodes(const string & s, const string & value, const string & unit, const UInt & pid, FIELD & Field) const
{
  //Open file
  string fileName = "sf" + num2str(pid) + s + ".inp";
  ofstream out(fileName.c_str());
  
  //Alloc
  Teuchos::RCP<MESH2D> grid2d = Field.getMesh2d();
  traitsDofType<DOFTYPE>  dofTrait;
  traitsGeoType<GEOSHAPE> geoTrait;
  
  //Print general data
  out << grid2d->getNumVertices() << " " << grid2d->getNumElements() << " 1 0 0" << endl;
  
  //Print mesh nodes
  for(UInt i=1; i<= grid2d->getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * grid2d->getNodeL(i).getX() << " ";
    out << paraFix * grid2d->getNodeL(i).getY() << " ";
    out << paraFix * grid2d->getNodeL(i).getZ() << endl;
  }
  
  //Print element nodes
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    out << i << " ";
    out << grid2d->getElementL(i).getGeoId() << "  ";
    out << geoTrait.parariewString() << " ";
    
    for(UInt k=1; k <= grid2d->getElementL(i).size(); ++k)
    {
      out << grid2d->getElementL(i).getCid(k) << " ";
    }
    out << endl;
  }
  
  //Compute the dofs vector
  sVect<OUTTYPE> outVect = nodesEval(Field);
  
  //Print dofs
  dofTrait.paraviewString(out,value,unit);
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    out << i << " ";
    dofTrait.printDof(out,outVect(i));
  }
}

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printUCD_elements(const string & s, const string & value, const string & unit, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const
{
  //Open file
  string fileName = "sf" + num2str(pid) + s + ".inp";
  ofstream out(fileName.c_str());
  
  //Alloc
  Teuchos::RCP<MESH2D> grid2d = Field->getMesh2d();
  traitsDofType<DOFTYPE>  dofTrait;
  traitsGeoType<GEOSHAPE> geoTrait;
  
  //Print general data
  out << grid2d->getNumVertices() << " " << grid2d->getNumElements() << " 0 1 0" << endl;
  
  //Print mesh nodes
  for(UInt i=1; i<= grid2d->getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * grid2d->getNodeL(i).getX() << " ";
    out << paraFix * grid2d->getNodeL(i).getY() << " ";
    out << paraFix * grid2d->getNodeL(i).getZ() << endl;
  }
  
  //Print element nodes
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    out << i << " ";
    out << grid2d->getElementL(i).getGeoId() << "  ";
    out << geoTrait.parariewString() << " ";
    
    for(UInt k=1; k <= grid2d->getElementL(i).size(); ++k)
    {
      out << grid2d->getElementL(i).getCid(k) << " ";
    }
    out << endl;
  }
  
  //Compute the dofs vector
  sVect<OUTTYPE> outVect = elementsEval(*Field);
  
  //Print dofs
  dofTrait.paraviewString(out,value,unit);
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    out << i << " ";
    dofTrait.printDof(out,outVect(i));
  }
}

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printUCD_elements(const string & s, const string & value, const string & unit, const UInt & pid, FIELD & Field) const
{
  //Open file
  string fileName = "sf" + num2str(pid) + s + ".inp";
  ofstream out(fileName.c_str());
  
  //Alloc
  Teuchos::RCP<MESH2D> grid2d = Field.getMesh2d();
  traitsDofType<DOFTYPE>  dofTrait;
  traitsGeoType<GEOSHAPE> geoTrait;
  
  //Print general data
  out << grid2d->getNumVertices() << " " << grid2d->getNumElements() << " 0 1 0" << endl;
  
  //Print mesh nodes
  for(UInt i=1; i<= grid2d->getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * grid2d->getNodeL(i).getX() << " ";
    out << paraFix * grid2d->getNodeL(i).getY() << " ";
    out << paraFix * grid2d->getNodeL(i).getZ() << endl;
  }
  
  //Print element nodes
  for(UInt i=1; i <= grid2d->getNumElements(); ++i)
  {
    out << i << " ";
    out << grid2d->getElementL(i).getGeoId() << "  ";
    out << geoTrait.parariewString() << " ";
    
    for(UInt k=1; k <= grid2d->getElementL(i).size(); ++k)
    {
      out << grid2d->getElementL(i).getCid(k) << " ";
    }
    out << endl;
  }
  
  //Compute the dofs vector
  sVect<OUTTYPE> outVect = elementsEval(Field);
  
  //Print dofs
  dofTrait.paraviewString(out,value,unit);
  
  for(UInt i=1; i <= outVect.size(); ++i)
  {
    out << i << " ";
    dofTrait.printDof(out,outVect(i));
  }
}


//_________________________________________________________________________________________________
// PRINT TO HDF5
//-------------------------------------------------------------------------------------------------
template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printHDF5_mesh(const string & s, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const
{
  typedef printMesh2dHDF5<GEOSHAPE,PMAPTYPE> MESHPRINTER;
  Teuchos::RCP<MESH2D>  grid2d = Field->getMesh2d();
  MESHPRINTER::print(s,pid,grid2d);
}

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printHDF5_mesh(const string & s, const UInt & pid, FIELD & Field) const
{
  typedef printMesh2dHDF5<GEOSHAPE,PMAPTYPE> MESHPRINTER;
  Teuchos::RCP<MESH2D>  grid2d = Field->getMesh2d();
  MESHPRINTER::print(s,pid,grid2d);
}

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printHDF5_nodes(const string & s, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const
{
  //Alloc------------------------------------------------------------------------------------------
  typedef printMesh2dHDF5<GEOSHAPE,PMAPTYPE> MESHPRINTER;
  Teuchos::RCP<communicator> commDev = Field->getDofMapper().getCommDev();  
  Teuchos::RCP<MESH2D>        grid2d = Field->getMesh2d();
 
  //Grid test--------------------------------------------------------------------------------------
  mesh2dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid2d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid2d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid2d);
  
  //Mesh-------------------------------------------------------------------------------------------
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = MESHPRINTER::convertElements(*grid2d);
  EpetraExt::DistArray<double> nodes_epetra    = MESHPRINTER::convertNodes(*grid2d);
  
  //Data conversion--------------------------------------------------------------------------------
  typedef pVect<OUTTYPE,PMAPTYPE> OUTVECT;
  
  sVect<OUTTYPE> tempVect = nodesEval(*Field);           //Create the out vector
  OUTVECT data_morgana(grid2d->getNodes().getMapRef(), tempVect);
  
  pVectGlobalManip<OUTTYPE,PMAPTYPE> dataManip(commDev); //Remap the vector
  Epetra_Map dataMap = dataManip.vectorLinear(data_morgana);  
  
  //HDF5 print-------------------------------------------------------------------------------------
  traitsDofType<DOFTYPE>  dofTrait;
  
  commDev->barrier();
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  hdf5.Write("elements",elements_epetra);                   //Print the mesh
  hdf5.Write("nodes",nodes_epetra);
  
  dofTrait.printHDF5(hdf5, data_morgana, dataMap); //Print the fields
  
  hdf5.Close();
  
  commDev->barrier();

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
    
    dofTrait.printXDMF(xdfmFile,"Node",s,numNodesG);
    
    xdfmFile << "    </Grid>"                                                            << endl;
    xdfmFile << "  </Domain>"                                                            << endl;
    xdfmFile << "</Xdmf>"                                                                << endl;  
    
    xdfmFile.close();
  }
}

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printHDF5_nodes(const string & s, const UInt & pid, FIELD & Field) const
{
  //Alloc------------------------------------------------------------------------------------------
  typedef printMesh2dHDF5<GEOSHAPE,PMAPTYPE> MESHPRINTER;
  Teuchos::RCP<communicator> commDev = Field.getDofMapper().getCommDev();  
  Teuchos::RCP<MESH2D>        grid2d = Field.getMesh2d();
 
  //Grid test--------------------------------------------------------------------------------------
  mesh2dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid2d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid2d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid2d);
  
  //Mesh-------------------------------------------------------------------------------------------
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = MESHPRINTER::convertElements(*grid2d);
  EpetraExt::DistArray<double> nodes_epetra    = MESHPRINTER::convertNodes(*grid2d);
  
  //Data conversion--------------------------------------------------------------------------------
  typedef pVect<OUTTYPE,PMAPTYPE> OUTVECT;
  
  sVect<OUTTYPE> tempVect = nodesEval(Field);           //Create the out vector
  OUTVECT data_morgana(grid2d->getNodes().getMapRef(), tempVect);
  
  pVectGlobalManip<OUTTYPE,PMAPTYPE> dataManip(commDev); //Remap the vector
  Epetra_Map dataMap = dataManip.vectorLinear(data_morgana);  
  
  //HDF5 print-------------------------------------------------------------------------------------
  traitsDofType<DOFTYPE>  dofTrait;
  
  commDev->barrier();
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  hdf5.Write("elements",elements_epetra);                   //Print the mesh
  hdf5.Write("nodes",nodes_epetra);
  
  dofTrait.printHDF5(hdf5, data_morgana, dataMap); //Print the fields
  
  hdf5.Close();
  
  commDev->barrier();

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
    
    dofTrait.printXDMF(xdfmFile,"Node",s,numNodesG);
    
    xdfmFile << "    </Grid>"                                                            << endl;
    xdfmFile << "  </Domain>"                                                            << endl;
    xdfmFile << "</Xdmf>"                                                                << endl;  
    
    xdfmFile.close();
  }
}

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printHDF5_elements(const string & s, const UInt & pid, const Teuchos::RCP<FIELD> & Field) const
{
  //Alloc------------------------------------------------------------------------------------------
  typedef printMesh2dHDF5<GEOSHAPE,PMAPTYPE> MESHPRINTER;
  Teuchos::RCP<communicator> commDev = Field->getDofMapper().getCommDev();  
  Teuchos::RCP<MESH2D>        grid2d = Field->getMesh2d();
 
  //Grid test--------------------------------------------------------------------------------------
  mesh2dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid2d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid2d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid2d);
  
  //Mesh-------------------------------------------------------------------------------------------
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = MESHPRINTER::convertElements(*grid2d);
  EpetraExt::DistArray<double> nodes_epetra    = MESHPRINTER::convertNodes(*grid2d);
  
  //Data conversion--------------------------------------------------------------------------------
  typedef pVect<OUTTYPE,PMAPTYPE> OUTVECT;
  
  sVect<OUTTYPE> tempVect = elementsEval(*Field);            //Create the out vector
  OUTVECT data_morgana(grid2d->getElements().getMapRef(), tempVect);
  
  pVectGlobalManip<OUTTYPE,PMAPTYPE> dataManip(commDev); //Remap the vector
  Epetra_Map dataMap = dataManip.vectorLinear(data_morgana);  
  
  //HDF5 print-------------------------------------------------------------------------------------
  traitsDofType<DOFTYPE>  dofTrait;
  
  commDev->barrier();
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  hdf5.Write("elements",elements_epetra);                   //Print the mesh
  hdf5.Write("nodes",nodes_epetra);
  
  dofTrait.printHDF5(hdf5, data_morgana, dataMap); //Print the fields
  
  hdf5.Close();
  
  commDev->barrier();

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

template<typename FIELD>
void
feStaticFieldPrinter2d<FIELD>::
printHDF5_elements(const string & s, const UInt & pid, FIELD & Field) const
{
//Alloc------------------------------------------------------------------------------------------
  typedef printMesh2dHDF5<GEOSHAPE,PMAPTYPE> MESHPRINTER;
  Teuchos::RCP<communicator> commDev = Field.getDofMapper().getCommDev();  
  Teuchos::RCP<MESH2D>        grid2d = Field.getMesh2d();
 
  //Grid test--------------------------------------------------------------------------------------
  mesh2dGlobalManip<GEOSHAPE,PMAPTYPE,PMAPTYPE> gridTester(commDev);
  assert(gridTester.check(grid2d));
  UInt numElementsG = gridTester.getNumGlobalElements(grid2d);
  UInt numNodesG    = gridTester.getNumGlobalVertices(grid2d);
  
  //Mesh-------------------------------------------------------------------------------------------
  Epetra_MpiComm epetraComm(*commDev);
  traitsGeoType<GEOSHAPE> geoTrait;
  
  EpetraExt::DistArray<int>    elements_epetra = MESHPRINTER::convertElements(*grid2d);
  EpetraExt::DistArray<double> nodes_epetra    = MESHPRINTER::convertNodes(*grid2d);
  
  //Data conversion--------------------------------------------------------------------------------
  typedef pVect<OUTTYPE,PMAPTYPE> OUTVECT;
  
  sVect<OUTTYPE> tempVect = elementsEval(Field);            //Create the out vector
  OUTVECT data_morgana(grid2d->getElements().getMapRef(), tempVect);
  
  pVectGlobalManip<OUTTYPE,PMAPTYPE> dataManip(commDev); //Remap the vector
  Epetra_Map dataMap = dataManip.vectorLinear(data_morgana);  
  
  //HDF5 print-------------------------------------------------------------------------------------
  traitsDofType<DOFTYPE>  dofTrait;
  
  commDev->barrier();
  
  EpetraExt::HDF5 hdf5(epetraComm);
  hdf5.Create(s + ".h5");
  
  hdf5.Write("elements",elements_epetra);                   //Print the mesh
  hdf5.Write("nodes",nodes_epetra);
  
  dofTrait.printHDF5(hdf5, data_morgana, dataMap); //Print the fields
  
  hdf5.Close();
  
  commDev->barrier();

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
