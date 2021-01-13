/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef LOADMESH_HPP
#define LOADMESH_HPP

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "load.h"
#include "traitsMapItemFixer.hpp"
#include "../morganaGeometry/geoShapes.h"

#include "pVect.hpp"
#include "pGraph.hpp"

#include "geoShapes.h"
#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"

#include "pVectComm.hpp"
#include "pGraphComm.hpp"


using namespace std;


//_______________________________________________________________________________________________________
// NOT SPECIALIZED
//-------------------------------------------------------------------------------------------------------
/*! Not specialized */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class loadMesh : public load
{
};


//_______________________________________________________________________________________________________
// LINEAR LINE SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
class loadMesh<linearLine,ELMAP,NODEMAP> : public load
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh1d<linearLine,ELMAP,NODEMAP> MESH1D;
    typedef typename MESH1D::NODESVECT       NODESVECT;
    typedef typename MESH1D::GRAPH1D         GRAPH1D;
    //@}
    
  public:
    bool paramSet;
    UInt pid, printPid, numPids;
  
  public:
    loadMesh();
    void setParam(const UInt & Pid,
                  const UInt & PrintPid,
                  const UInt & NumPids);
    
  public:
    void amira(string meshfile,
               NODESVECT & nodes,
               GRAPH1D   & elements);
    
    void gmMesh(string meshfile,
                NODESVECT & nodes,
                GRAPH1D   & elements);
};


template<typename ELMAP, typename NODEMAP>
loadMesh<linearLine,ELMAP,NODEMAP>::
loadMesh() : load()
{
  paramSet = false;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearLine,ELMAP,NODEMAP>::
setParam(const UInt & Pid,
         const UInt & PrintPid,
         const UInt & NumPids)
{
  paramSet = true;
  
  pid      = Pid;
  printPid = PrintPid;
  numPids  = NumPids;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearLine,ELMAP,NODEMAP>::
amira(string meshfile,
      NODESVECT & nodes,
      GRAPH1D   & elements)
{
  //Assert-----------------------------------------------------------------------------------------
  assert(paramSet);
  
  //Typedefs---------------------------------------------------------------------------------------
  typedef linearLine              GEOSHAPE1D;
  typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
  
  //Clear data-------------------------------------------------------------------------------------
  nodes.clear();
  elements.clear();
  
  //Read numbers-----------------------------------------------------------------------------------
  UInt n_nodes, dummy;
  UInt nElements1d = 0;
  
  //File open
  ifstream file(meshfile.c_str());
  
  //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile << endl; }
  
  //Skip of comments
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  
  //Load numbers
  leggiriga(file);
  assegnavalore(2, n_nodes);
  
  leggiriga(file);
  assegnavalore(2, nElements1d);
  
  //Skip of comments
  for(UInt i=1; i<=14; ++i)
  { leggiriga(file); }
  
  //Load nodes-------------------------------------------------------------------------------------
  UInt upperBound = segmentationUpper(pid,numPids,n_nodes);
  UInt lowerBound = segmentationLower(pid,numPids,n_nodes);
  
  UInt k = 1;
  Real X, Y, Z;
  
  for(UInt i=1; i <= n_nodes; ++i)
  {
    leggiriga(file);
    assegnavalore(0,X);     //Assegna coordinata X
    assegnavalore(1,Y);     //Assegna coordinata Y
    assegnavalore(2,Z);     //Assegna coordinata Z
    
    if( (i <= upperBound) && (i >= lowerBound) )
    {
      nodes.push_back(point3d(X,Y,Z),NODEMAP(k,i));
      k++;
    }
  }
  
  //Load elements----------------------------------------------------------------------------------
  leggiriga(file);
  leggiriga(file);
  
  UInt upperBound1d = segmentationUpper(pid,numPids,nElements1d);
  UInt lowerBound1d = segmentationLower(pid,numPids,nElements1d);
  
  UInt id1, id2;
  UInt elIdLoc1d = 1;
  
  ELMAP elMapItem;
  GEOELEMENT1D element1d(true);
  elements.reserve(upperBound1d-lowerBound1d);
  
  for(UInt elId1d=1; elId1d <= nElements1d; ++elId1d)
  {
    leggiriga(file);
    
    if( (elId1d <= upperBound1d) && (elId1d >= lowerBound1d) )
    {
      assegnavalore(0,id1);
      assegnavalore(1,id2);

      element1d(1) = id1 + 1;
      element1d(2) = id2 + 1;
      element1d.setGeoId(0);
      
      elMapItem.setLid(elIdLoc1d);
      elMapItem.setGid(elId1d);
      mapItemFixer(elMapItem);
      
      elements.push_back(element1d,elMapItem);
      elIdLoc1d++;
    }
  }
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearLine,ELMAP,NODEMAP>::
gmMesh(string meshfile,
       NODESVECT & nodes,
       GRAPH1D   & elements)
{
  //Assert----------------------------------------------------------------------------------------
  assert(paramSet);
  
  //Typedefs---------------------------------------------------------------------------------------
  typedef linearLine              GEOSHAPE1D;
  typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;

  //Clear data-------------------------------------------------------------------------------------
  nodes.clear();
  elements.clear();
  
  //Numbers reading--------------------------------------------------------------------------------
  UInt type;
  UInt n_nodes, dummy;
  UInt nTotElements;
  UInt nElements1d = 0;
  
  //File open
  ifstream file(meshfile.c_str());
  
  //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile << endl; }
  
  //Skip comments
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  
  //Skipping nodes
  leggiriga(file);
  assegnavalore(0, n_nodes); //Number of nodes
  
  for(UInt i=1; i <= n_nodes; ++i)
  { leggiriga(file); }

  //Skip comments
  leggiriga(file);
  leggiriga(file);
  
  //Skip corner nodes
  leggiriga(file);
  assegnavalore(0,nTotElements); //Number of elements
  
  for(UInt i=1; i <= nTotElements; ++i)
  {    
    leggiriga(file);
    assegnavalore(0,dummy); // Numerazione gmMesh
    assegnavalore(1,type);  // Tipo di elemento geometrico
    
    switch(type)
    {
      //Is a corner
      case 15 :
      break;
      
      //Is an edge
      case 1 :
      nElements1d++;
      break;

      //Not defined type
      default :
      if(pid == printPid) { cout << "ERROR! data loading failed, gmsh-type  not supported";}
      assert(1 == 2);
    }
  }

  //File closing and primary data output
  file.close();
  
  if(pid == printPid)
  {
    cout << "Estimating mesh size" << endl;
    cout << "Nodes       - global : " << n_nodes << endl;
    cout << "Elements 1d - global : " << nElements1d << endl;
    cout << "Now loading" << endl;
  }
  
  //Now loading------------------------------------------------------------------------------------
  
  //File open
  file.open(meshfile.c_str());
  
   //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Comments skip
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  
  //Read nodes
  leggiriga(file);
  
  UInt upperBound = segmentationUpper(pid,numPids,n_nodes);
  UInt lowerBound = segmentationLower(pid,numPids,n_nodes);
  
  UInt k = 1;
  UInt id1, id2, id3, sub;
  Real X,Y,Z;
  
  for(UInt i=1; i <= n_nodes; ++i)
  {
    leggiriga(file);
    assegnavalore(0,dummy); //Indice nodo ignorato
    assegnavalore(1,X);     //Assegna coordinata X
    assegnavalore(2,Y);     //Assegna coordinata Y
    assegnavalore(3,Z);     //Assegna coordinata Z
    
    if( (i <= upperBound) && (i >= lowerBound) )
    {
      nodes.push_back(point3d(X,Y,Z),NODEMAP(k,i));
      k++;
    }
  }
  
  //Skip comments
  leggiriga(file);
  leggiriga(file);
  
  //Skip corner nodes
  UInt upperBound1d = segmentationUpper(pid,numPids,nElements1d);
  UInt lowerBound1d = segmentationLower(pid,numPids,nElements1d);
  
  UInt numRowEl;
  UInt elId1d = 1, elIdLoc1d = 1;
  
  GEOELEMENT1D element1d(true);
  
  ELMAP elMapItem;
  elements.reserve(upperBound1d-lowerBound1d);
  
  leggiriga(file); //Skip total number of elements
  
  
  //Loading loop
  for(UInt i=1; i <= nTotElements; ++i)
  {    
    numRowEl = leggiriga(file);
    assegnavalore(0,dummy); // Numerazione gmMesh
    assegnavalore(1,type);  // Tipo di elemento geometrico
    
    switch (type)
    {
      //Is a corner
      case 15 :
      break;
      
      //Is 1d element
      case 1 :
        if( (elId1d <= upperBound1d) && (elId1d >= lowerBound1d) )
        {
          assegnavalore(2,dummy); // Numero di flags
          assert((dummy == 3) || (dummy == 2));

          if(dummy == 3)
          {
            assert(numRowEl == 8);
            assegnavalore(3,dummy); // Elimina il primo flag
            assegnavalore(4,sub);   // Sottodominio di appartenenza
            assegnavalore(5,dummy); // Elimina il terzo flag
            assegnavalore(6,id1);   // Primo id
            assegnavalore(7,id2);   // Secondo id
          }

          if(dummy == 2)
          {
            assert(numRowEl == 7);
            assegnavalore(3,dummy); // Elimina il primo flag
            assegnavalore(4,sub);   // Sottodominio di appartenenza
            assegnavalore(5,id1);   // Primo id
            assegnavalore(6,id2);   // Secondo id
          }

          element1d(1) = id1;
          element1d(2) = id2;
          element1d.setGeoId(sub);

          elMapItem.setLid(elIdLoc1d);
          elMapItem.setGid(elId1d);
          mapItemFixer(elMapItem);
    
          elements.push_back(element1d,elMapItem);
          elIdLoc1d++;
        }

        elId1d++;
        break;

      //Not defined type
      default :
      if(pid == printPid) { cout << "ERROR! data loading failed, gmsh-type  not supported";}
      assert(1 == 2);
    }
  }
  
  elements.colIsLocal() = false;

  //Closing
  file.close();
  
  if(pid == printPid)
  { cout << "Mesh loaded successfully" << endl; }
}


//_______________________________________________________________________________________________________
// LINEAR TRIANGLE SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
template<typename ELMAP, typename NODEMAP>
class loadMesh<linearTriangle,ELMAP,NODEMAP> : public load
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh2d<linearTriangle,ELMAP,NODEMAP> MESH2D;
    typedef typename MESH2D::NODESVECT           NODESVECT;
    typedef typename MESH2D::GRAPH2D             GRAPH2D;
    typedef typename MESH2D::GRAPH1D             GRAPH1D;
    //@}
    
  public:
    bool paramSet;
    UInt pid, printPid, numPids;
  
  public:
    loadMesh();
    void setParam(const UInt & Pid,
                  const UInt & PrintPid,
                  const UInt & NumPids);
    
  public:
    void gmMesh(string meshfile,
                NODESVECT & nodes,
                GRAPH2D   & elements,
                GRAPH1D   & Belements);
    
    void gmMeshParallel(const string & meshfile,
                        const Teuchos::RCP<const communicator> & commDev,
                        NODESVECT & nodes,
                        GRAPH2D   & elements,
                        GRAPH1D   & Belements);
    
    void offMesh(string meshfile,
                 NODESVECT & nodes,
                 GRAPH2D   & elements);
};


template<typename ELMAP, typename NODEMAP>
loadMesh<linearTriangle,ELMAP,NODEMAP>::
loadMesh() : load()
{
  paramSet = false;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearTriangle,ELMAP,NODEMAP>::
setParam(const UInt & Pid,
         const UInt & PrintPid,
         const UInt & NumPids)
{
  paramSet = true;
  
  pid      = Pid;
  printPid = PrintPid;
  numPids  = NumPids;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearTriangle,ELMAP,NODEMAP>::
gmMesh(string meshfile,
       NODESVECT & nodes,
       GRAPH2D   & elements,
       GRAPH1D   & Belements)
{
  assert(paramSet);
  
  typedef linearTriangle GEOSHAPE2D;
  typedef typename       GEOSHAPE2D::GEOBSHAPE  GEOSHAPE1D;
    
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
  
  
  //Clear data_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  Belements.clear();
  
  //NUMBERS READING________________________________________________________________________________
  UInt type;
  UInt n_nodes, dummy;
  UInt nTotElements;
  UInt nElements1d = 0;
  UInt nElements2d = 0;
  
  //File open
  ifstream file(meshfile.c_str());
  
  //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile       << endl; }
  
  //Skip of comments
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  
  //Skipping nodes
  leggiriga(file);
  assegnavalore(0, n_nodes); //Number of nodes
  
  for(UInt i=1; i <= n_nodes; ++i)
  {
    leggiriga(file);
  }
  
  //Skip comments
  leggiriga(file);
  leggiriga(file);
  
  //Skip corner nodes
  leggiriga(file);
  assegnavalore(0,nTotElements); //Number of elements
  
  for(UInt i=1; i <= nTotElements; ++i)
  {    
    leggiriga(file);
    assegnavalore(0,dummy); // Numerazione gmMesh
    assegnavalore(1,type);  // Tipo di elemento geometrico
    
    switch (type)
    {
      //Is a corner
      case 15 :
	break;
      
      //Is an edge
      case 1 :
	nElements1d++;
	break;
	
      //Is a 2d element
      case 2 :
	nElements2d++;
	break;
	
        //Not defined type
      default :
	if(pid == printPid)
	{ cout << "ERROR! data loading failed, gmsh-type  not supported";}
	assert(1 == 2);
    }
  }
  
  //File closing and primary data output
  file.close();
  
  if(pid == printPid)
  {
    cout << "Estimating mesh size" << endl;
    cout << "Nodes       - global : " << n_nodes << endl;
    cout << "Elements 1d - global : " << nElements1d << endl;
    cout << "Elements 2d - global : " << nElements2d << endl;
    cout << "Now loading" << endl;
  }
  
  
  
  //NOW LOADING____________________________________________________________________________________
  
  //File open
  file.open(meshfile.c_str());
  
  //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Skip dei commenti______________________________________________________________________________
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  
  //Lettura nodi___________________________________________________________________________________ 
  leggiriga(file);
  
  UInt upperBound = segmentationUpper(pid,numPids,n_nodes);
  UInt lowerBound = segmentationLower(pid,numPids,n_nodes);
  
  UInt k = 1;
  UInt id1, id2, id3, sub;
  Real X,Y,Z;  
  
  for(UInt i=1; i <= n_nodes; ++i)
  {
    leggiriga(file);
    assegnavalore(0,dummy); //Indice nodo ignorato
    assegnavalore(1,X);     //Assegna coordinata X
    assegnavalore(2,Y);     //Assegna coordinata Y
    assegnavalore(3,Z);     //Assegna coordinata Z
    
    if( (i <= upperBound) && (i >= lowerBound) )
    {
      nodes.push_back(point3d(X,Y,Z),NODEMAP(k,i));
      k++;
    }
  }
  
  //Skip comments__________________________________________________________________________________
  leggiriga(file);
  leggiriga(file);
  
  //Skip corner nodes______________________________________________________________________________
  UInt upperBound1d = segmentationUpper(pid,numPids,nElements1d);
  UInt lowerBound1d = segmentationLower(pid,numPids,nElements1d);
  
  UInt upperBound2d = segmentationUpper(pid,numPids,nElements2d);
  UInt lowerBound2d = segmentationLower(pid,numPids,nElements2d);
  
  UInt elId1d = 1, elIdLoc1d = 1;
  UInt elId2d = 1, elIdLoc2d = 1;
  
  GEOELEMENT1D element1d(true);
  GEOELEMENT2D element2d(true);
  
  ELMAP elMapItem;
  Belements.reserve(upperBound1d-lowerBound1d);
  elements.reserve(upperBound2d-lowerBound2d);
  
  UInt numRowEl;
  leggiriga(file); //Skip total number of elements
  
  for(UInt i=1; i <= nTotElements; ++i)
  {    
    numRowEl = leggiriga(file);
    assegnavalore(0,dummy); // Numerazione gmMesh
    assegnavalore(1,type);  // Tipo di elemento geometrico
    
    switch (type)
    {
      //Is a corner
      case 15 :
	break;
      
	
      //Is 1d element
      case 1 :
	if( (elId1d <= upperBound1d) && (elId1d >= lowerBound1d) )
        {
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 8);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 7);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	  }
	  
          element1d(1) = id1;
          element1d(2) = id2;
          element1d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc1d);
	  elMapItem.setGid(elId1d);
	  mapItemFixer(elMapItem);
	  
	  Belements.push_back(element1d,elMapItem);
	  
	  elIdLoc1d++;
	}
	
	elId1d++;
	break;
	
	
      //Is a 2d element
      case 2 :
	if( (elId2d <= upperBound2d) && (elId2d >= lowerBound2d) )
        {
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 9);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	    assegnavalore(8,id3);   // Terzo id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 8);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	    assegnavalore(7,id3);   // Terzo id
	  }
	  
          element2d(1) = id1;
          element2d(2) = id2;
          element2d(3) = id3;
          element2d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc2d);
	  elMapItem.setGid(elId2d);
	  mapItemFixer(elMapItem);
	  
	  elements.push_back(element2d,elMapItem);
	  
	  elIdLoc2d++;
	}
	
	elId2d++;
	break;
	
	
        //Not defined type
      default :
	if(pid == printPid)
	{ cout << "ERROR! data loading failed, gmsh-type  not supported";}
	assert(1 == 2);
    }
  }
  
  elements.colIsLocal() = false;
  Belements.colIsLocal() = false;
  
  
  //Closing
  file.close();
  
  if(pid == printPid)
  {
    cout << "Mesh loaded successfully" << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearTriangle,ELMAP,NODEMAP>::
gmMeshParallel(const string & meshfile,
               const Teuchos::RCP<const communicator> & commDev,
               NODESVECT & nodes,
               GRAPH2D   & elements,
               GRAPH1D   & Belements)
{
  //Typedefs_______________________________________________________________________________________
  typedef linearTriangle GEOSHAPE2D;
  typedef typename       GEOSHAPE2D::GEOBSHAPE  GEOSHAPE1D;
    
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
  
  //Clear data_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  Belements.clear();
  
  //Comm classes___________________________________________________________________________________
  pVectComm<point3d,NODEMAP>             vectorComm(commDev);
  pGraphComm<GEOELEMENT1D,ELMAP,NODEMAP> elComm1d(commDev);
  pGraphComm<GEOELEMENT2D,ELMAP,NODEMAP> elComm2d(commDev);
  
  
  //____________________________________________________________________________________________
  // MASTER PID LOADING
  //-----------------------------------------------------------------------------------------------
  if(commDev->rank() == 0)
  {
    //Elements counting----------------------------------------------
    UInt type;
    UInt n_nodes, dummy;
    UInt nTotElements;
    UInt nElements1d = 0;
    UInt nElements2d = 0;
    
    //File open
    ifstream file(meshfile.c_str());
  
    //File check
    if(!file.good())
    { cout << "ERROR! File not found. Pid: " << pid << endl; }
    assert(file.good());
  
    //Printout file
    cout << "LOADING FROM: " << meshfile << endl;
  
    //Skip of comments
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
    
    //Skipping nodes
    leggiriga(file);
    assegnavalore(0, n_nodes); //Number of nodes
    
    for(UInt i=1; i <= n_nodes; ++i)
    { leggiriga(file); }
  
    //Skip comments
    leggiriga(file);
    leggiriga(file);
  
    //Skip corner nodes
    leggiriga(file);
    assegnavalore(0,nTotElements); //Number of elements
    
    for(UInt i=1; i <= nTotElements; ++i)
    {    
      leggiriga(file);
      assegnavalore(0,dummy); // Numerazione gmMesh
      assegnavalore(1,type);  // Tipo di elemento geometrico
    
      switch (type)
      {
        case 15 : break;                 //Is a corner
        case 1  : nElements1d++; break;  //Is a 1d element
        case 2  : nElements2d++; break;  //Is a 2d element
        default : cout << "ERROR! data loading failed, gmsh-type  not supported"; assert(1 == 2);
      }
    }
    
    //File closing and primary data output
    file.close();
  
    //Printout
    cout << "Estimating mesh size" << endl;
    cout << "Nodes       - global : " << n_nodes << endl;
    cout << "Elements 1d - global : " << nElements1d << endl;
    cout << "Elements 2d - global : " << nElements2d << endl;
    cout << "Now loading" << endl;
    
    
    //Elements loading-----------------------------------------------
    file.open(meshfile.c_str()); //File open
  
    //File check
    if(!file.good())
    { cout << "ERROR! File not found." << endl; }
    assert(file.good());
    
    //Skip dei commenti______________________________________________________________________________
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
  
    //Lettura nodi___________________________________________________________________________________ 
    leggiriga(file);
    
    UInt upperBound = segmentationUpper(0,commDev->size(),n_nodes);
    UInt k = 1, pidLoc = 0;
    UInt id1, id2, id3, sub;
    Real X,Y,Z;
    
    NODESVECT nodesTemp;
    sVect<request> reqs(1);
    
    for(UInt i=1; i <= n_nodes; ++i)
    {
      leggiriga(file);
      assegnavalore(0,dummy); //Indice nodo ignorato
      assegnavalore(1,X);     //Assegna coordinata X
      assegnavalore(2,Y);     //Assegna coordinata Y
      assegnavalore(3,Z);     //Assegna coordinata Z
      
      nodesTemp.push_back(point3d(X,Y,Z),NODEMAP(k,i));
      k++;
      
      if(i == upperBound)
      {
	if(pidLoc == 0)
	{ nodes = nodesTemp; }
	else
	{ vectorComm.send(0,pidLoc,nodesTemp,reqs); }
	
	k = 1;
	pidLoc++;
	upperBound = segmentationUpper(pidLoc,commDev->size(),n_nodes);
	nodesTemp.clear();
      }
    }
    
    //Skip comments__________________________________________________________________________________
    leggiriga(file);
    leggiriga(file);
    
    //Read elements__________________________________________________________________________________
    GRAPH2D elementsTemp2d;
    GRAPH1D elementsTemp1d;
    
    UInt upperBound1d = segmentationUpper(0,commDev->size(),nElements1d);
    UInt upperBound2d = segmentationUpper(0,commDev->size(),nElements2d);

    UInt elIdGlob1d = 1, elIdLoc1d = 1;
    UInt elIdGlob2d = 1, elIdLoc2d = 1;
  
    GEOELEMENT1D element1d(true);
    GEOELEMENT2D element2d(true);
    
    UInt pid1d = 0, pid2d = 0;
    ELMAP elMapItem;
    
    UInt numRowEl;
    leggiriga(file); //Skip total number of elements
    
    //Loop on the elements
    for(UInt i=1; i <= nTotElements; ++i)
    {    
      numRowEl = leggiriga(file);
      assegnavalore(0,dummy); // Numerazione gmMesh
      assegnavalore(1,type);  // Tipo di elemento geometrico
    
      switch (type)
      {
	case 15 : break; //Is a corner

        //Is a 1d element
	case 1 :
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 8);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 7);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	  }
	  
	  element1d(1) = id1;
          element1d(2) = id2;
          element1d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc1d);
	  elMapItem.setGid(elIdGlob1d);
	  mapItemFixer(elMapItem);
	  
	  elementsTemp1d.push_back(element1d,elMapItem);
	  elIdLoc1d++;
	  elIdGlob1d++;
	  
	  if(elIdGlob1d > upperBound1d)
	  {
	    if(pid1d == 0)
	    { Belements = elementsTemp1d; }
	    else
	    { elComm1d.send(0,pid1d,elementsTemp1d,reqs); }
	    
	    elIdLoc1d = 1;
	    pid1d++;
	    upperBound1d = segmentationUpper(pid1d,commDev->size(),nElements1d);
	    elementsTemp1d.clear();
	  }
	  break;
	  
	//Is a 2d element
	case 2 :
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 9);
	    assegnavalore(3,dummy); // Elimina il primo flag 
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	    assegnavalore(8,id3);   // Terzo id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 8);
	    assegnavalore(3,dummy); // Elimina il primo flag 
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	    assegnavalore(7,id3);   // Terzo id
	  }
	  
	
	  element2d(1) = id1;
          element2d(2) = id2;
          element2d(3) = id3;
          element2d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc2d);
	  elMapItem.setGid(elIdGlob2d);
	  mapItemFixer(elMapItem);
	  
	  elementsTemp2d.push_back(element2d,elMapItem);
	  elIdLoc2d++;
	  elIdGlob2d++;
	  
	  if(elIdGlob2d > upperBound2d)
	  {
	    if(pid2d == 0)
	    { elements = elementsTemp2d; }
	    else
	    { elComm2d.send(0,pid2d,elementsTemp2d,reqs); }
	    
	    elIdLoc2d = 1;
	    pid2d++;
	    upperBound2d = segmentationUpper(pid2d,commDev->size(),nElements2d);
	    elementsTemp2d.clear();
	  }
	  break;
      } //End switch
    } //End loop elements
    
    //Closing
    file.close();
  }
  
  //_______________________________________________________________________________________________
  // SLAVES RECEIVING
  //-----------------------------------------------------------------------------------------------
  else
  {
    vectorComm.recv(0,commDev->rank(),nodes);
    elComm1d.recv(0,commDev->rank(),Belements);
    elComm2d.recv(0,commDev->rank(),elements);
  }
  
  
  //Sincronization-----------------------------------------
  commDev->barrier(); 
  
  if(pid == printPid)
  { cout << "Mesh loaded successfully" << endl; }
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearTriangle,ELMAP,NODEMAP>::
offMesh(string meshfile,
        NODESVECT & nodes,
        GRAPH2D   & elements)
{
  //Assert_________________________________________________________________________________________
  assert(paramSet);
  
  //Typedefs_______________________________________________________________________________________
  typedef linearTriangle          GEOSHAPE2D; 
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;

  //Clear data_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  
  //Alloc__________________________________________________________________________________________
  UInt n_nodes, dummy;
  UInt nElements2d = 0;
  ELMAP elMapItem;
  
  //NOW LOADING____________________________________________________________________________________
  
  //File open
  ifstream file(meshfile.c_str());
  
  //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Skip dei comments
  leggiriga(file);
  
  //Main parameters
  leggiriga(file);
  assegnavalore(0,n_nodes);
  assegnavalore(1,nElements2d);
  
  //Load nodes
  UInt upperBound = segmentationUpper(pid,numPids,n_nodes);
  UInt lowerBound = segmentationLower(pid,numPids,n_nodes);
  
  UInt k = 1;
  Real X,Y,Z;  
  
  for(UInt i=1; i <= n_nodes; ++i)
  {
    leggiriga(file);
    assegnavalore(0,X);     //Assegna coordinata X
    assegnavalore(1,Y);     //Assegna coordinata Y
    assegnavalore(2,Z);     //Assegna coordinata Z
    
    if( (i <= upperBound) && (i >= lowerBound) )
    {
      nodes.push_back(point3d(X,Y,Z),NODEMAP(k,i));
      k++;
    }
  }
  
  //Load faces
  UInt upperBound2d = segmentationUpper(pid,numPids,nElements2d);
  UInt lowerBound2d = segmentationLower(pid,numPids,nElements2d);
  
  GEOELEMENT2D element2d(true);
  elements.reserve(upperBound2d-lowerBound2d);
  
  UInt id1, id2, id3;
  k = 1;
  
  for(UInt i=1; i <= nElements2d; ++i)
  {
    leggiriga(file);
    
    if( (i <= upperBound2d) && (i >= lowerBound2d) )
    {
      assegnavalore(1,id1);
      assegnavalore(2,id2);
      assegnavalore(3,id3);
      
      element2d(1) = id1 + 1;
      element2d(2) = id2 + 1;
      element2d(3) = id3 + 1;
      element2d.setGeoId(1);
      
      elMapItem.setLid(k);
      elMapItem.setGid(i);
      mapItemFixer(elMapItem);
  
      elements.push_back(element2d,elMapItem);
  
      k++;
    }
  }
  
  elements.colIsLocal() = false;
  
  //Closing
  file.close();
  
  if(pid == printPid)
  {
    cout << "Mesh loaded successfully" << endl;
  }
}


//_______________________________________________________________________________________________________
// LINEAR QUAD SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Quad 2d */
template<typename ELMAP, typename NODEMAP>
class loadMesh<linearQuad,ELMAP,NODEMAP> : public load
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh2d<linearQuad,ELMAP,NODEMAP>  MESH2D;
    typedef typename MESH2D::NODESVECT        NODESVECT;
    typedef typename MESH2D::GRAPH2D          GRAPH2D;
    typedef pVect<UInt,NODEMAP>               NODESCOLOR;
    //@}
    
  public:
    bool paramSet;
    UInt pid, printPid, numPids;
  
  public:
    loadMesh();
    void setParam(const UInt & Pid,
                  const UInt & PrintPid,
                  const UInt & NumPids);
    
  public:
    void nastran(string meshfile,
                 NODESVECT  & nodes,
                 GRAPH2D    & elements,
                 NODESCOLOR & nodesColor);
    
    void neutral(string meshfile,
                 NODESVECT  & nodes,
                 GRAPH2D    & elements,
                 NODESCOLOR & nodesColor);
};


template<typename ELMAP, typename NODEMAP>
loadMesh<linearQuad,ELMAP,NODEMAP>::
loadMesh() : load()
{
  paramSet = false;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearQuad,ELMAP,NODEMAP>::
setParam(const UInt & Pid,
         const UInt & PrintPid,
         const UInt & NumPids)
{
  paramSet = true;
  
  pid      = Pid;
  printPid = PrintPid;
  numPids  = NumPids;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearQuad,ELMAP,NODEMAP>::
nastran(string meshfile,
        NODESVECT  & nodes,
        GRAPH2D    & elements,
        NODESCOLOR & nodesColor)
{
  assert(paramSet);
  
  typedef linearQuad              GEOSHAPE2D;    
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  
  //File open______________________________________________________________________________________
  ifstream file(meshfile.c_str());
  
  //File check_____________________________________________________________________________________
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file__________________________________________________________________________________
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile       << endl; }
  
  //Clear dati_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  nodesColor.clear();
  
  //Main loading cycle
  UInt numNodes = 0, numElements = 0, blocks = 0;
  UInt dataSet, color;
  UInt id1, id2, id3, id4;
  Real dummy, X, Y, Z;
  GEOELEMENT2D element2d(true);
  
  sVect<point3d>      tempNodes;
  sVect<UInt>         tempColor;
  sVect<GEOELEMENT2D> temp2d;
  
  
  while(!file.eof())
  {
    // -1 di inizio ciclo
    leggiriga(file);
    assegnavalore(0,dummy);
    assert(int(dummy) == -1);
    
    //Data set
    leggiriga(file);
    assegnavalore(0,dataSet);
    
    //Cases
    switch(dataSet)
    {
      //The nodes list
      case 2411 :
	while(!file.eof())
	{
	  leggiriga(file);
	  assegnavalore(0,dummy);
	  
	  if(int(dummy) == -1)
	  {break;}
	  
	  //Assign node color
	  assegnavalore(3,color);
	  tempColor.push_back(color);
	  
	  //Assign coordinates
	  leggiriga(file);
	  
	  assegnavalore(0,X);
	  assegnavalore(1,Y);
	  assegnavalore(2,Z);
	  
	  tempNodes.push_back(point3d(X,Y,Z));
          numNodes++;
	}
	break;
	
	
      //The Elements list
      case 2412 :
	while(!file.eof())
	{
	  leggiriga(file);
	  assegnavalore(0,dummy);
	  
	  if(int(dummy) == -1)
	  {break;}
	  
	  //Color assignment
	  assegnavalore(4,color);
	  assegnavalore(5,dummy);
	  
	  assert(dummy == 4);
	  
	  //Nodes loading
	  leggiriga(file);
	  
	  assegnavalore(0,id1);
	  assegnavalore(1,id2);
	  assegnavalore(2,id3);
	  assegnavalore(3,id4);
	  
	  element2d(1) = id1;
          element2d(2) = id2;
          element2d(3) = id3;
	  element2d(4) = id4;

          element2d.setGeoId(color);
          temp2d.push_back(element2d);
	  
	  numElements++;
	}
	break;
      
      //Skipping
      default :
	do
	{
	  leggiriga(file);
	  assegnavalore(0,dummy);
	}
	while(int(dummy) != -1);
    }
    
    //Increment the number of blocks  
    blocks++;
  }
  
  //Cout number of data blocks
  if(pid == printPid)
  {
    cout << "Read " << blocks << " data blocks" << endl;
    cout << "Segmentation procedure " << endl << endl;
  }
  
  
  //Nodes segmentation_____________________________________________________________________________
  UInt upperBound = segmentationUpper(pid,numPids,numNodes);
  UInt lowerBound = segmentationLower(pid,numPids,numNodes);
  
  UInt  k = 1;
  NODEMAP nodeMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    nodeMapItem.setLid(k);
    nodeMapItem.setGid(i);
    nodeMapItem.setPid(pid);
    
    ++k;
    mapItemFixer(nodeMapItem);
    
    nodes.push_back(tempNodes(i),nodeMapItem);
    nodesColor.push_back(tempColor(i),nodeMapItem);
  }
  
  nodes.updateFinder();
  nodesColor.updateFinder();
  
  //Elements3d segmentation________________________________________________________________________
  upperBound = segmentationUpper(pid,numPids,numElements);
  lowerBound = segmentationLower(pid,numPids,numElements);
  
  k = 1;
  ELMAP elMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    elMapItem.setLid(k);
    elMapItem.setGid(i);
    
    ++k;
    mapItemFixer(elMapItem);
    
    elements.push_back(temp2d(i),elMapItem);
  }
  
  elements.updateRowFinder();
  elements.colIsLocal() = false;
  
  
  //File closing___________________________________________________________________________________
  file.close();
  
  if(pid == printPid)
  {
    cout << "Mesh loaded successfully" << endl;
    cout << "Nodes       - global : " << numNodes << endl;
    cout << "Elements 2d - global : " << numElements<< endl << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearQuad,ELMAP,NODEMAP>::
neutral(string meshfile,
        NODESVECT  & nodes,
        GRAPH2D    & elements,
        NODESCOLOR & nodesColor)
{
  assert(paramSet);
  
  typedef linearQuad              GEOSHAPE2D;    
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  
  //File open______________________________________________________________________________________
  ifstream file(meshfile.c_str());
  
  //File check_____________________________________________________________________________________
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file__________________________________________________________________________________
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile       << endl; }
  
  //Clear dati_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  nodesColor.clear();
  
  //Main loading cycle
  bool globalTerm = false;
  UInt numNodes = 0, numElements = 0, blocks = 0;
  UInt dataSet, color, numRead;
  UInt id1, id2, id3, id4;
  Real dummy, X, Y, Z;
  GEOELEMENT2D element2d(true);
  
  sVect<point3d>      tempNodes;
  sVect<UInt>         tempColor;
  sVect<GEOELEMENT2D> temp2d;
  
  
  while(!file.eof())
  {
    // -1 di inizio ciclo
    numRead = leggiriga(file);
    assegnavalore(0,dummy);
    
    assert((int(dummy) == -1) && (numRead == 1));
    
    //Data set
    leggiriga(file);
    assegnavalore(0,dataSet);
    
    //Cases
    switch(dataSet)
    {
      //The nodes list
      case 403 :
	if(pid == printPid)
	{cout << "Nodes Loading" << endl;}
	
	while(!file.eof())
	{
	  numRead = leggiriga(file);
	  assegnavalore(0,dummy);
	  
	  if( (int(dummy) == -1) && (numRead == 1))
	  {break;}
	  
	  //Assign node color
	  assegnavalore(4,color);
	  tempColor.push_back(color);
	  
	  //Assign coordinates	  
	  assegnavalore(11,X);
	  assegnavalore(12,Y);
	  assegnavalore(13,Z);
	  
	  tempNodes.push_back(point3d(X,Y,Z));
          numNodes++;
	}
	break;
	
	
      //The Elements list
      case 404 :
	if(pid == printPid)
	{cout << "Elements Loading" << endl;}
	
	while(!file.eof())
	{
	  numRead = leggiriga(file);          //Riga 1
	  assegnavalore(0,dummy);
	  
	  if((int(dummy) == -1) && (numRead == 1))
	  {
	    globalTerm = true;
	    break;
	  }
	  
	  //Color assignment
	  assegnavalore(1,color);
	  
	  //Nodes loading
	  leggiriga(file);          //Riga 2
  
	  assegnavalore(0,id4);
	  assegnavalore(1,id1);
	  assegnavalore(2,id2);
	  assegnavalore(3,id3);
	  
	  element2d(1) = id1;
          element2d(2) = id2;
          element2d(3) = id3;
	  element2d(4) = id4;

          element2d.setGeoId(color);
          temp2d.push_back(element2d);
	  
	  numElements++;
	  
	  //Dump lines
	  leggiriga(file);          //Riga 3
	  leggiriga(file);          //Riga 4
	  leggiriga(file);          //Riga 5
	  leggiriga(file);          //Riga 6
	  leggiriga(file);          //Riga 7
	}
	break;
      
      //Skipping
      default :
	do
	{
	  numRead = leggiriga(file);
	  assegnavalore(0,dummy);
	}
	while(! ((int(dummy) == -1) && (numRead == 1)) );
    }
    
    if(globalTerm)
    {
      break;
    }
    
    //Increment the number of blocks  
    blocks++;
  }
  
  //Cout number of data blocks
  if(pid == printPid)
  {
    cout << "Read " << blocks << " data blocks" << endl;
    cout << "Segmentation procedure " << endl << endl;
  }
  
  
  //Nodes segmentation_____________________________________________________________________________
  UInt upperBound = segmentationUpper(pid,numPids,numNodes);
  UInt lowerBound = segmentationLower(pid,numPids,numNodes);
  
  UInt  k = 1;
  NODEMAP nodeMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    nodeMapItem.setLid(k);
    nodeMapItem.setGid(i);
    nodeMapItem.setPid(pid);
    
    ++k;
    mapItemFixer(nodeMapItem);
    
    nodes.push_back(tempNodes(i),nodeMapItem);
    nodesColor.push_back(tempColor(i),nodeMapItem);
  }
  
  nodes.updateFinder();
  nodesColor.updateFinder();
  
  //Elements3d segmentation________________________________________________________________________
  upperBound = segmentationUpper(pid,numPids,numElements);
  lowerBound = segmentationLower(pid,numPids,numElements);
  
  k = 1;
  ELMAP elMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    elMapItem.setLid(k);
    elMapItem.setGid(i);
    
    ++k;
    mapItemFixer(elMapItem);
    
    elements.push_back(temp2d(i),elMapItem);
  }
  
  elements.updateRowFinder();
  elements.colIsLocal() = false;
  
  
  //File closing___________________________________________________________________________________
  file.close();
  
  if(pid == printPid)
  {
    cout << "Mesh loaded successfully" << endl;
    cout << "Nodes       - global : " << numNodes << endl;
    cout << "Elements 2d - global : " << numElements<< endl << endl;
  }
}


//_______________________________________________________________________________________________________
// LINEAR TETRA SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Tetrahedral 3d */
template<typename ELMAP, typename NODEMAP>
class loadMesh<linearTetra,ELMAP,NODEMAP> : public load
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh3d<linearTetra,ELMAP,NODEMAP> MESH3D;
    typedef typename MESH3D::NODESVECT        NODESVECT;
    typedef typename MESH3D::GRAPH3D          GRAPH3D;
    typedef typename MESH3D::GRAPH2D          GRAPH2D;
    //@}
    
  public:
    bool paramSet;
    UInt pid, printPid, numPids;
  
  public:
    loadMesh();
    void setParam(const UInt & Pid,
                  const UInt & PrintPid,
                  const UInt & NumPids);
    
  public:
    void gmMesh(string meshfile,
                NODESVECT & nodes,
                GRAPH3D   & elements,
                GRAPH2D   & Belements);
    
    void gmMeshParallel(const string & meshfile,
                        const Teuchos::RCP<const communicator> & commDev,
                        NODESVECT & nodes,
                        GRAPH3D   & elements,
                        GRAPH2D   & Belements);
};


template<typename ELMAP, typename NODEMAP>
loadMesh<linearTetra,ELMAP,NODEMAP>::
loadMesh() : load()
{
  paramSet = false;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearTetra,ELMAP,NODEMAP>::
setParam(const UInt & Pid,
         const UInt & PrintPid,
         const UInt & NumPids)
{
  paramSet = true;
  
  pid      = Pid;
  printPid = PrintPid;
  numPids  = NumPids;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearTetra,ELMAP,NODEMAP>::
gmMesh(string meshfile,
       NODESVECT & nodes,
       GRAPH3D   & elements,
       GRAPH2D   & Belements)
{
  //Asserts________________________________________________________________________________________
  assert(paramSet);
  
  //Typedefs_______________________________________________________________________________________
  typedef linearTetra GEOSHAPE3D;
  typedef typename    GEOSHAPE3D::GEOBSHAPE  GEOSHAPE2D;
  typedef typename    GEOSHAPE2D::GEOBSHAPE  GEOSHAPE1D;
    
  typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
  
  //Clear data_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  Belements.clear();
  
  //NUMBERS READING________________________________________________________________________________
  UInt type;
  UInt n_nodes, dummy;
  UInt nTotElements;
  UInt nElements2d = 0;
  UInt nElements3d = 0;
  
  //File open
  ifstream file(meshfile.c_str());
  
  //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile       << endl; }
  
  //Skip of comments
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  
  //Skipping nodes
  leggiriga(file);
  assegnavalore(0, n_nodes); //Number of nodes
  
  for(UInt i=1; i <= n_nodes; ++i)
  {
    leggiriga(file);
  }
  
  //Skip comments
  leggiriga(file);
  leggiriga(file);
  
  //Skip corner nodes
  leggiriga(file);
  assegnavalore(0,nTotElements); //Number of elements
  
  for(UInt i=1; i <= nTotElements; ++i)
  {    
    leggiriga(file);
    assegnavalore(0,dummy); // Numerazione gmMesh
    assegnavalore(1,type);  // Tipo di elemento geometrico
    
    switch (type)
    {
      //Is a corner
      case 15 :
	break;
      
      //Is an edge
      case 1 :
	break;
	
      //Is a 2d element
      case 2 :
	nElements2d++;
	break;
	
	//Is a 3d element
      case 4 :
	nElements3d++;
	break;
	
        //Not defined type
      default :
	if(pid == printPid)
	{ cout << "ERROR! data loading failed, gmsh-type  not supported";}
	assert(1 == 2);
    }
  }
  
  //File closing and primary data output
  file.close();
  
  if(pid == printPid)
  {
    cout << "Estimating mesh size" << endl;
    cout << "Nodes       - global : " << n_nodes << endl;
    cout << "Elements 2d - global : " << nElements2d << endl;
    cout << "Elements 3d - global : " << nElements3d << endl;
    cout << "Now loading" << endl;
  }
  
  
  //NOW LOADING____________________________________________________________________________________
  
  //File open
  file.open(meshfile.c_str());
  
  //File check
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Skip dei commenti______________________________________________________________________________
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  leggiriga(file);
  
  //Lettura nodi___________________________________________________________________________________ 
  leggiriga(file);
  
  UInt upperBound = segmentationUpper(pid,numPids,n_nodes);
  UInt lowerBound = segmentationLower(pid,numPids,n_nodes);
  
  UInt k = 1;
  UInt id1, id2, id3, id4, sub;
  Real X,Y,Z;  
  
  for(UInt i=1; i <= n_nodes; ++i)
  {
    leggiriga(file);
    assegnavalore(0,dummy); //Indice nodo ignorato
    assegnavalore(1,X);     //Assegna coordinata X
    assegnavalore(2,Y);     //Assegna coordinata Y
    assegnavalore(3,Z);     //Assegna coordinata Z
    
    if( (i <= upperBound) && (i >= lowerBound) )
    {
      nodes.push_back(point3d(X,Y,Z),NODEMAP(k,i));
      k++;
    }
  }
  
  //Skip comments__________________________________________________________________________________
  leggiriga(file);
  leggiriga(file);
  
  //Skip corner nodes______________________________________________________________________________
  UInt upperBound2d = segmentationUpper(pid,numPids,nElements2d);
  UInt lowerBound2d = segmentationLower(pid,numPids,nElements2d);
  
  UInt upperBound3d = segmentationUpper(pid,numPids,nElements3d);
  UInt lowerBound3d = segmentationLower(pid,numPids,nElements3d);
  
  UInt elId2d = 1, elIdLoc2d = 1;
  UInt elId3d = 1, elIdLoc3d = 1;
  
  GEOELEMENT2D element2d(true);
  GEOELEMENT3D element3d(true);
  
  ELMAP elMapItem;
  Belements.reserve(upperBound2d-lowerBound2d);
  elements.reserve(upperBound3d-lowerBound3d);
  
  UInt numRowEl;
  leggiriga(file); //Skip total number of elements
  
  for(UInt i=1; i <= nTotElements; ++i)
  {    
    numRowEl = leggiriga(file);
    assegnavalore(0,dummy); // Numerazione gmMesh
    assegnavalore(1,type);  // Tipo di elemento geometrico
    
    switch (type)
    {
      //Is a corner
      case 15 :
	break;
      
      //Is an edge
      case 1 :
	break;
	
      //Is a 2d element
      case 2 :
	if( (elId2d <= upperBound2d) && (elId2d >= lowerBound2d) )
        {
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 9);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	    assegnavalore(8,id3);   // Terzo id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 8);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	    assegnavalore(7,id3);   // Terzo id
	  }
	  
          element2d(1) = id1;
          element2d(2) = id2;
          element2d(3) = id3;
          element2d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc2d);
	  elMapItem.setGid(elId2d);
	  mapItemFixer(elMapItem);
	  
	  Belements.push_back(element2d,elMapItem);
	  
	  elIdLoc2d++;
	}
	
	elId2d++;
	break;
	
	//Is a 3d element
      case 4 :
	if( (elId3d <= upperBound3d) && (elId3d >= lowerBound3d) )
        {
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 10);
	    assegnavalore(3,dummy); // Elimina il primo flag 
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	    assegnavalore(8,id3);   // Terzo id
	    assegnavalore(9,id4);   // Quarto id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 9);
	    assegnavalore(3,dummy); // Elimina il primo flag 
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	    assegnavalore(7,id3);   // Terzo id
	    assegnavalore(8,id4);   // Quarto id
	  }
	  
	
	  element3d(1) = id1;
          element3d(2) = id2;
          element3d(3) = id3;
	  element3d(4) = id4;
          element3d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc3d);
	  elMapItem.setGid(elId3d);
	  mapItemFixer(elMapItem);
	  
	  elements.push_back(element3d,elMapItem);
	  
	  elIdLoc3d++;
	}
	
	elId3d++;
	break;
	
        //Not defined type
      default :
	if(pid == printPid)
	{ cout << "ERROR! data loading failed, gmsh-type  not supported";}
	assert(1 == 2);
    }
  }
  
  elements.colIsLocal() = false;
  Belements.colIsLocal() = false;
  
  
  //Closing
  file.close();
  
  if(pid == printPid)
  {
    cout << "Mesh loaded successfully" << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearTetra,ELMAP,NODEMAP>::
gmMeshParallel(const string & meshfile,
               const Teuchos::RCP<const communicator> & commDev,
               NODESVECT & nodes,
               GRAPH3D   & elements,
               GRAPH2D   & Belements)
{  
  //Typedefs_______________________________________________________________________________________
  typedef linearTetra GEOSHAPE3D;
  typedef typename    GEOSHAPE3D::GEOBSHAPE  GEOSHAPE2D;
  typedef typename    GEOSHAPE2D::GEOBSHAPE  GEOSHAPE1D;
    
  typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
  
  //Clear data_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  Belements.clear();
  
  //Comm classes___________________________________________________________________________________
  pVectComm<point3d,NODEMAP>     vectorComm(commDev);
  pVectComm<GEOELEMENT2D,ELMAP>  elComm2d(commDev);
  pVectComm<GEOELEMENT3D,ELMAP>  elComm3d(commDev);
  
  
  //_______________________________________________________________________________________________
  // MASTER PID LOADING
  //-----------------------------------------------------------------------------------------------
  if(commDev->rank() == 0)
  {
    //Elements counting----------------------------------------------
    UInt type;
    UInt n_nodes, dummy;
    UInt nTotElements;
    UInt nElements2d = 0;
    UInt nElements3d = 0;
    
    //File open
    ifstream file(meshfile.c_str());
  
    //File check
    if(!file.good())
    { cout << "ERROR! File not found. Pid: " << pid << endl; }
    assert(file.good());
  
    //Printout file
    cout << "LOADING FROM: " << meshfile << endl;
  
    //Skip of comments
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
    
    //Skipping nodes
    leggiriga(file);
    assegnavalore(0, n_nodes); //Number of nodes
  
    for(UInt i=1; i <= n_nodes; ++i)
    { leggiriga(file); }
  
    //Skip comments
    leggiriga(file);
    leggiriga(file);
  
    //Skip corner nodes
    leggiriga(file);
    assegnavalore(0,nTotElements); //Number of elements
  
    for(UInt i=1; i <= nTotElements; ++i)
    {    
      leggiriga(file);
      assegnavalore(0,dummy); // Numerazione gmMesh
      assegnavalore(1,type);  // Tipo di elemento geometrico
    
      switch (type)
      {
        case 15 : break;                //Is a corner
        case 1  : break;                //Is an edge
        case 2  : nElements2d++; break;  //Is a 2d element
        case 4  : nElements3d++; break;  //Is a 3d element
        default : cout << "ERROR! data loading failed, gmsh-type  not supported"; assert(1 == 2);
      }
    }
  
    //File closing and primary data output
    file.close();
  
    //Printout
    cout << "Estimating mesh size" << endl;
    cout << "Nodes       - global : " << n_nodes << endl;
    cout << "Elements 2d - global : " << nElements2d << endl;
    cout << "Elements 3d - global : " << nElements3d << endl;
    cout << "Now loading" << endl;
    
    
    //Elements loading-----------------------------------------------
       
    //File open
    file.open(meshfile.c_str());
  
    //File check
    if(!file.good())
    { cout << "ERROR! File not found." << endl; }
    assert(file.good());
    
    //Skip dei commenti______________________________________________________________________________
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
    leggiriga(file);
  
    //Lettura nodi___________________________________________________________________________________ 
    leggiriga(file);
    
    UInt upperBound = segmentationUpper(0,commDev->size(),n_nodes);
    UInt k = 1, pidLoc = 0;
    UInt id1, id2, id3, id4, sub;
    Real X,Y,Z;
    NODESVECT nodesTemp;
    
    for(UInt i=1; i <= n_nodes; ++i)
    {
      leggiriga(file);
      assegnavalore(0,dummy); //Indice nodo ignorato
      assegnavalore(1,X);     //Assegna coordinata X
      assegnavalore(2,Y);     //Assegna coordinata Y
      assegnavalore(3,Z);     //Assegna coordinata Z
      
      nodesTemp.push_back(point3d(X,Y,Z),NODEMAP(k,i));
      k++;
      
      if(i == upperBound)
      {
	if(pidLoc == 0)
	{ nodes = nodesTemp; }
	else
	{ vectorComm.send(0,pidLoc,nodesTemp); }
	
	k = 1;
	pidLoc++;
	upperBound = segmentationUpper(pidLoc,commDev->size(),n_nodes);
	nodesTemp.clear();
      }
    }
    
    //Skip comments__________________________________________________________________________________
    leggiriga(file);
    leggiriga(file);
    
    //Read elements__________________________________________________________________________________
    GRAPH3D elementsTemp3d;
    GRAPH2D elementsTemp2d;
    
    UInt upperBound2d = segmentationUpper(0,commDev->size(),nElements2d);
    UInt upperBound3d = segmentationUpper(0,commDev->size(),nElements3d);

    UInt elIdGlob2d = 1, elIdLoc2d = 1;
    UInt elIdGlob3d = 1, elIdLoc3d = 1;
  
    GEOELEMENT2D element2d(true);
    GEOELEMENT3D element3d(true);
    
    UInt pid2d = 0, pid3d = 0;
    ELMAP elMapItem;
    
    UInt numRowEl;
    leggiriga(file); //Skip total number of elements
    
    
    //Loop on the elements
    for(UInt i=1; i <= nTotElements; ++i)
    {    
      numRowEl = leggiriga(file);
      assegnavalore(0,dummy); // Numerazione gmMesh
      assegnavalore(1,type);  // Tipo di elemento geometrico
    
      switch (type)
      {
	case 15 : break; //Is a corner
        case 1  : break; //Is an edge

        //Is a 2d element
	case 2 :
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 9);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	    assegnavalore(8,id3);   // Terzo id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 8);
	    assegnavalore(3,dummy); // Elimina il primo flag
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	    assegnavalore(7,id3);   // Terzo id
	  }
	  
	  element2d(1) = id1;
          element2d(2) = id2;
          element2d(3) = id3;
          element2d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc2d);
	  elMapItem.setGid(elIdGlob2d);
	  mapItemFixer(elMapItem);
	  
	  elementsTemp2d.push_back(element2d,elMapItem);
	  elIdLoc2d++;
	  elIdGlob2d++;
	  
	  if(elIdGlob2d > upperBound2d)
	  {
	    if(pid2d == 0)
	    { Belements = elementsTemp2d; }
	    else
	    { elComm2d.send(0,pid2d,elementsTemp2d); }
	    
	    elIdLoc2d = 1;
	    pid2d++;
	    upperBound2d = segmentationUpper(pid2d,commDev->size(),nElements2d);
	    elementsTemp2d.clear();
	  }
	  break;
	  
	//Is a 3d element
	case 4 :
	  assegnavalore(2,dummy); // Numero di flags
	  assert((dummy == 3) || (dummy == 2));
	  
	  if(dummy == 3)
	  {
	    assert(numRowEl == 10);
	    assegnavalore(3,dummy); // Elimina il primo flag 
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,dummy); // Elimina il terzo flag
	    assegnavalore(6,id1);   // Primo id
	    assegnavalore(7,id2);   // Secondo id
	    assegnavalore(8,id3);   // Terzo id
	    assegnavalore(9,id4);   // Quarto id
	  }
	  
	  if(dummy == 2)
	  {
	    assert(numRowEl == 9);
	    assegnavalore(3,dummy); // Elimina il primo flag 
	    assegnavalore(4,sub);   // Sottodominio di appartenenza
	    assegnavalore(5,id1);   // Primo id
	    assegnavalore(6,id2);   // Secondo id
	    assegnavalore(7,id3);   // Terzo id
	    assegnavalore(8,id4);   // Quarto id
	  }
	  
	
	  element3d(1) = id1;
          element3d(2) = id2;
          element3d(3) = id3;
	  element3d(4) = id4;
          element3d.setGeoId(sub);
	  
	  elMapItem.setLid(elIdLoc3d);
	  elMapItem.setGid(elIdGlob3d);
	  mapItemFixer(elMapItem);
	  
	  elementsTemp3d.push_back(element3d,elMapItem);
	  elIdLoc3d++;
	  elIdGlob3d++;
	  
	  if(elIdGlob3d > upperBound3d)
	  {
	    if(pid3d == 0)
	    { elements = elementsTemp3d; }
	    else
	    { elComm3d.send(0,pid3d,elementsTemp3d); }
	    
	    elIdLoc3d = 1;
	    pid3d++;
	    upperBound3d = segmentationUpper(pid3d,commDev->size(),nElements3d);
	    elementsTemp3d.clear();
	  }
	  break;
      } //End switch
    } //End loop elements
    
    //Closing
    file.close();
  }
  else //RECEIVING-----------------------------------------------------------------------
  {
    vectorComm.recv(0,commDev->rank(),nodes);
    elComm2d.recv(0,commDev->rank(),Belements);
    elComm3d.recv(0,commDev->rank(),elements);
  }
  
  commDev->barrier();
  
  if(commDev->rank() == 0)
  { cout << "Mesh loaded successfully" << endl; }
}


//_______________________________________________________________________________________________________
// LINEAR HEXA SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Hexaedral 3d */
template<typename ELMAP, typename NODEMAP>
class loadMesh<linearHexa,ELMAP,NODEMAP> : public load
{
    /*! @name Typedefs */ //@{
  public:
    typedef mesh3d<linearHexa,ELMAP,NODEMAP>  MESH3D;
    typedef typename MESH3D::NODESVECT        NODESVECT;
    typedef typename MESH3D::GRAPH3D          GRAPH3D;
    typedef pVect<UInt,NODEMAP>               NODESCOLOR;
    //@}
    
  public:
    bool paramSet;
    UInt pid, printPid, numPids;
  
  public:
    loadMesh();
    void setParam(const UInt & Pid,
                  const UInt & PrintPid,
                  const UInt & NumPids);
    
  public:
    void nastran(string meshfile,
                 NODESVECT  & nodes,
                 GRAPH3D    & elements,
                 NODESCOLOR & nodesColor);
    
    void neutral(string meshfile,
                 NODESVECT  & nodes,
                 GRAPH3D    & elements,
                 NODESCOLOR & nodesColor);
};

template<typename ELMAP, typename NODEMAP>
loadMesh<linearHexa,ELMAP,NODEMAP>::
loadMesh() : load()
{
  paramSet = false;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearHexa,ELMAP,NODEMAP>::
setParam(const UInt & Pid,
         const UInt & PrintPid,
         const UInt & NumPids)
{
  paramSet = true;
  
  pid      = Pid;
  printPid = PrintPid;
  numPids  = NumPids;
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearHexa,ELMAP,NODEMAP>::
nastran(string meshfile,
        NODESVECT  & nodes,
        GRAPH3D    & elements,
        NODESCOLOR & nodesColor)
{
  assert(paramSet);
  
  typedef linearHexa              GEOSHAPE3D;    
  typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
  
  //File open______________________________________________________________________________________
  ifstream file(meshfile.c_str());
  
  //File check_____________________________________________________________________________________
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file__________________________________________________________________________________
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile       << endl; }
  
  //Clear dati_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  nodesColor.clear();
  
  //Main loading cycle
  UInt numNodes = 0, numElements = 0, blocks = 0;
  UInt dataSet, color, sizeLoad;
  UInt id1, id2, id3, id4, id5, id6, id7, id8;
  Real dummy, X, Y, Z;
  GEOELEMENT3D element3d(true);
  
  sVect<point3d>      tempNodes;
  sVect<UInt>         tempColor;
  sVect<GEOELEMENT3D> temp3d;
  
  
  while(!file.eof())
  {
    // -1 di inizio ciclo
    sizeLoad = leggiriga(file);
    if(sizeLoad == 0) {break;}
    
    assegnavalore(0,dummy);
    assert(int(dummy) == -1);
    
    //Data set
    leggiriga(file);
    assegnavalore(0,dataSet);
    
    //Cases
    switch(dataSet)
    {
      //The nodes list
      case 2411 :
	while(!file.eof())
	{
	  leggiriga(file);
	  assegnavalore(0,dummy);
	  
	  if(int(dummy) == -1)
	  {break;}
	  
	  //Assign node color
	  assegnavalore(3,color);
	  tempColor.push_back(color);
	  
	  //Assign coordinates
	  leggiriga(file);
	  
	  assegnavalore(0,X);
	  assegnavalore(1,Y);
	  assegnavalore(2,Z);
	  
	  tempNodes.push_back(point3d(X,Y,Z));
          numNodes++;
	}
	break;
	
	
      //The Elements list
      case 2412 :
	while(!file.eof())
	{
	  leggiriga(file);
	  assegnavalore(0,dummy);
	  
	  if(int(dummy) == -1)
	  {break;}
	  
	  //Color assignment
	  assegnavalore(4,color);
	  assegnavalore(5,dummy);
	  
	  assert(dummy == 8);
	  
	  //Nodes loading
	  leggiriga(file);
	  
	  assegnavalore(0,id1);
	  assegnavalore(1,id2);
	  assegnavalore(2,id3);
	  assegnavalore(3,id4);
	  
	  assegnavalore(4,id5);
	  assegnavalore(5,id6);
	  assegnavalore(6,id7);
	  assegnavalore(7,id8);
	  
	  element3d(1) = id1;
          element3d(2) = id2;
          element3d(3) = id3;
	  element3d(4) = id4;
	  
	  element3d(5) = id5;
          element3d(6) = id6;
          element3d(7) = id7;
	  element3d(8) = id8;
	  
          element3d.setGeoId(color);
          temp3d.push_back(element3d);
	  
	  numElements++;
	}
	break;
      
      //Skipping
      default :
	do
	{
	  leggiriga(file);
	  assegnavalore(0,dummy);
	}
	while(int(dummy) != -1);
    }
    
    //Increment the number of blocks  
    blocks++;
  }
  
  //Cout number of data blocks
  if(pid == printPid)
  {
    cout << "Read " << blocks << " data blocks" << endl;
    cout << "Segmentation procedure " << endl << endl;
  }
  
  
  //Nodes segmentation_____________________________________________________________________________
  UInt upperBound = segmentationUpper(pid,numPids,numNodes);
  UInt lowerBound = segmentationLower(pid,numPids,numNodes);
  
  UInt  k = 1;
  NODEMAP nodeMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    nodeMapItem.setLid(k);
    nodeMapItem.setGid(i);
    nodeMapItem.setPid(pid);
    
    ++k;
    mapItemFixer(nodeMapItem);
    
    nodes.push_back(tempNodes(i),nodeMapItem);
    nodesColor.push_back(tempColor(i),nodeMapItem);
  }
  
  nodes.updateFinder();
  nodesColor.updateFinder();
  
  //Elements3d segmentation________________________________________________________________________
  upperBound = segmentationUpper(pid,numPids,numElements);
  lowerBound = segmentationLower(pid,numPids,numElements);
  
  k = 1;
  ELMAP elMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    elMapItem.setLid(k);
    elMapItem.setGid(i);
    
    ++k;
    mapItemFixer(elMapItem);
    
    elements.push_back(temp3d(i),elMapItem);
  }
  
  elements.updateRowFinder();
  elements.colIsLocal() = false;
  
  
  //File closing___________________________________________________________________________________
  file.close();
  
  if(pid == printPid)
  {
    cout << "Mesh loaded successfully" << endl;
    cout << "Nodes       - global : " << numNodes << endl;
    cout << "Elements 3d - global : " << numElements<< endl << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
loadMesh<linearHexa,ELMAP,NODEMAP>::
neutral(string meshfile,
        NODESVECT  & nodes,
        GRAPH3D    & elements,
        NODESCOLOR & nodesColor)
{
  assert(paramSet);
  
  typedef linearHexa              GEOSHAPE3D;    
  typedef geoElement<GEOSHAPE3D>  GEOELEMENT3D;
  
  //File open______________________________________________________________________________________
  ifstream file(meshfile.c_str());
  
  //File check_____________________________________________________________________________________
  if(!file.good())
  { cout << "ERROR! File not found. Pid: " << pid << endl; }
  
  assert(file.good());
  
  //Printout file__________________________________________________________________________________
  if(pid == printPid)
  { cout << "LOADING FROM: " << meshfile       << endl; }
  
  //Clear dati_____________________________________________________________________________________
  nodes.clear();
  elements.clear();
  nodesColor.clear();
  
  //Main loading cycle
  bool globalTerm = false;
  UInt numNodes = 0, numElements = 0, blocks = 0;
  UInt dataSet, color, numRead;
  UInt id1, id2, id3, id4, id5, id6, id7, id8;
  Real dummy, X, Y, Z;
  GEOELEMENT3D element3d(true);
  
  sVect<point3d>      tempNodes;
  sVect<UInt>         tempColor;
  sVect<GEOELEMENT3D> temp3d;
  
  
  while(!file.eof())
  {
    // -1 di inizio ciclo
    numRead = leggiriga(file);
    assegnavalore(0,dummy);
    
    assert((int(dummy) == -1) && (numRead == 1));
    
    //Data set
    leggiriga(file);
    assegnavalore(0,dataSet);
    
    //Cases
    switch(dataSet)
    {
      //The nodes list
      case 403 :
	if(pid == printPid)
	{cout << "Nodes Loading" << endl;}
	
	while(!file.eof())
	{
	  numRead = leggiriga(file);
	  assegnavalore(0,dummy);
	  
	  if( (int(dummy) == -1) && (numRead == 1))
	  {break;}
	  
	  //Assign node color
	  assegnavalore(4,color);
	  tempColor.push_back(color);
	  
	  //Assign coordinates	  
	  assegnavalore(11,X);
	  assegnavalore(12,Y);
	  assegnavalore(13,Z);
	  
	  tempNodes.push_back(point3d(X,Y,Z));
          numNodes++;
	}
	break;
	
	
      //The Elements list
      case 404 :
	if(pid == printPid)
	{cout << "Elements Loading" << endl;}
	
	while(!file.eof())
	{
	  numRead = leggiriga(file);          //Riga 1
	  assegnavalore(0,dummy);
	  
	  if((int(dummy) == -1) && (numRead == 1))
	  {
	    globalTerm = true;
	    break;
	  }
	  
	  //Color assignment
	  assegnavalore(1,color);
	  
	  //Nodes loading
	  leggiriga(file);          //Riga 2
  
	  assegnavalore(0,id1);
	  assegnavalore(1,id2);
	  assegnavalore(2,id3);
	  assegnavalore(3,id4);
	  
	  assegnavalore(4,id5);
	  assegnavalore(5,id6);
	  assegnavalore(6,id7);
	  assegnavalore(7,id8);
	  
	  element3d(1) = id1;
          element3d(2) = id2;
          element3d(3) = id3;
	  element3d(4) = id4;
	  
	  element3d(5) = id5;
          element3d(6) = id6;
          element3d(7) = id7;
	  element3d(8) = id8;

          element3d.setGeoId(color);
          temp3d.push_back(element3d);
	  
	  numElements++;
	  
	  //Dump lines
	  leggiriga(file);          //Riga 3
	  leggiriga(file);          //Riga 4
	  leggiriga(file);          //Riga 5
	  leggiriga(file);          //Riga 6
	  leggiriga(file);          //Riga 7
	}
	break;
      
      //Skipping
      default :
	do
	{
	  numRead = leggiriga(file);
	  assegnavalore(0,dummy);
	}
	while(! ((int(dummy) == -1) && (numRead == 1)) );
    }
    
    if(globalTerm)
    {
      break;
    }
    
    //Increment the number of blocks  
    blocks++;
  }
  
  //Cout number of data blocks
  if(pid == printPid)
  {
    cout << "Read " << blocks << " data blocks" << endl;
    cout << "Segmentation procedure " << endl << endl;
  }
  
  
  //Nodes segmentation_____________________________________________________________________________
  UInt upperBound = segmentationUpper(pid,numPids,numNodes);
  UInt lowerBound = segmentationLower(pid,numPids,numNodes);
  
  UInt  k = 1;
  NODEMAP nodeMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    nodeMapItem.setLid(k);
    nodeMapItem.setGid(i);
    nodeMapItem.setPid(pid);
    
    ++k;
    mapItemFixer(nodeMapItem);
    
    nodes.push_back(tempNodes(i),nodeMapItem);
    nodesColor.push_back(tempColor(i),nodeMapItem);
  }
  
  nodes.updateFinder();
  nodesColor.updateFinder();
  
  //Elements3d segmentation________________________________________________________________________
  upperBound = segmentationUpper(pid,numPids,numElements);
  lowerBound = segmentationLower(pid,numPids,numElements);
  
  k = 1;
  ELMAP elMapItem;
  
  for(UInt i = lowerBound; i <= upperBound; ++i)
  {
    elMapItem.setLid(k);
    elMapItem.setGid(i);
    
    ++k;
    mapItemFixer(elMapItem);
    
    elements.push_back(temp3d(i),elMapItem);
  }
  
  elements.updateRowFinder();
  elements.colIsLocal() = false;
  
  
  //File closing___________________________________________________________________________________
  file.close();
  
  if(pid == printPid)
  {
    cout << "Mesh loaded successfully" << endl;
    cout << "Nodes       - global : " << numNodes << endl;
    cout << "Elements 2d - global : " << numElements<< endl << endl;
  }
}

#endif
