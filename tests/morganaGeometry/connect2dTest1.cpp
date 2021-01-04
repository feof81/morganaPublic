/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <cmath>
#include <iostream>

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "typesInterface.hpp"

#include "pMapItem.h"
#include "pMapItemShare.h"

#include "pGraph.hpp"

#include "geoShapes.h"
#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "mesh1dGlobalManip.hpp"
#include "mesh2dGlobalManip.hpp"
#include "connect2d.hpp"
#include "meshDoctor3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
  typedef linearLine              GEOSHAPE1D;
  typedef linearQuad              GEOSHAPE2D;
  typedef geoElement<GEOSHAPE1D>  GEOELEMENT1D;
  typedef geoElement<GEOSHAPE2D>  GEOELEMENT2D;
  typedef pMapItemShare ELMAP;
  typedef pMapItem      NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP> MESH2D;
  typedef mesh1d<GEOSHAPE1D,ELMAP,NODEMAP> MESH1D;
  
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  Teuchos::RCP<MESH2D> grid2d = Teuchos::rcp(new MESH2D);
  Teuchos::RCP<MESH1D> grid1d = Teuchos::rcp(new MESH1D);
  
  if(world.rank() == 0)
  {
    //Maps
    NODEMAP nodeMapItem;
    ELMAP   elMapItem;
    
    nodeMapItem.setPid(0);
    elMapItem.setPid(0);
    
    //Nodes2d
    pVect<point3d,NODEMAP> nodes2d;
    
    nodes2d.reserve(8);
    nodeMapItem.setLid(1); nodeMapItem.setGid(1);
    nodes2d.push_back(point3d(0.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(2); nodeMapItem.setGid(2);
    nodes2d.push_back(point3d(1.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(3); nodeMapItem.setGid(4);
    nodes2d.push_back(point3d(0.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(4); nodeMapItem.setGid(5);
    nodes2d.push_back(point3d(1.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(5); nodeMapItem.setGid(7);
    nodes2d.push_back(point3d(0.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(6); nodeMapItem.setGid(8);
    nodes2d.push_back(point3d(1.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(7); nodeMapItem.setGid(6);
    nodes2d.push_back(point3d(2.0, 1.0, 0.0),nodeMapItem);

    nodeMapItem.setLid(8); nodeMapItem.setGid(9);
    nodes2d.push_back(point3d(2.0, 2.0, 0.0),nodeMapItem);
    
    
    //Elements2d
    GEOELEMENT2D quad(true);
    pGraph<GEOELEMENT2D,ELMAP,NODEMAP> elList2d;
    
    elList2d.reserve(3);
    quad.setGeoId(1);
    
    elMapItem.setLid(1); elMapItem.setGid(1); elMapItem.setShared(false); elMapItem.setOwned(true);
    quad(1) = 1; quad(2) = 2; quad(3) = 4; quad(4) = 3;
    elList2d.push_back(quad,elMapItem);
    
    elMapItem.setLid(2); elMapItem.setGid(2); elMapItem.setShared(false); elMapItem.setOwned(true);
    quad(1) = 3; quad(2) = 4; quad(3) = 6; quad(4) = 5;
    elList2d.push_back(quad,elMapItem);
    
    elMapItem.setLid(3); elMapItem.setGid(4); elMapItem.setShared(true); elMapItem.setOwned(false);
    quad(1) = 4; quad(2) = 7; quad(3) = 8; quad(4) = 6;
    elList2d.push_back(quad,elMapItem);
    
    elList2d.updateRowFinder();
    elList2d.updateColFinder();
    
    
    //The grid2d
    grid2d->setNodes(nodes2d);
    grid2d->setElements(elList2d);
    
    
    //Nodes1d
    pVect<point3d,NODEMAP> nodes1d;
    
    nodes1d.reserve(7);
    nodeMapItem.setLid(1); nodeMapItem.setGid(8);
    nodes1d.push_back(point3d(1.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(2); nodeMapItem.setGid(1);
    nodes1d.push_back(point3d(0.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(3); nodeMapItem.setGid(2);
    nodes1d.push_back(point3d(0.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(4); nodeMapItem.setGid(3);
    nodes1d.push_back(point3d(0.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(5); nodeMapItem.setGid(4);
    nodes1d.push_back(point3d(1.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(6); nodeMapItem.setGid(5);
    nodes1d.push_back(point3d(2.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(7); nodeMapItem.setGid(6);
    nodes1d.push_back(point3d(2.0, 1.0, 0.0),nodeMapItem);
    
    
    //Elements1d
    GEOELEMENT1D li(true);
    pGraph<GEOELEMENT1D,ELMAP,NODEMAP> elList1d;
    
    elList1d.reserve(6);
    
    elMapItem.setLid(1); elMapItem.setGid(8); elMapItem.setShared(false); elMapItem.setOwned(true); li.setGeoId(1);
    li(1) = 1; li(2) = 2;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(2); elMapItem.setGid(1); elMapItem.setShared(false); elMapItem.setOwned(true); li.setGeoId(1);
    li(1) = 2; li(2) = 3;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(3); elMapItem.setGid(2); elMapItem.setShared(false); elMapItem.setOwned(true); li.setGeoId(1);
    li(1) = 3; li(2) = 4;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(4); elMapItem.setGid(3); elMapItem.setShared(false); elMapItem.setOwned(true); li.setGeoId(1);
    li(1) = 4; li(2) = 5;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(5); elMapItem.setGid(4); elMapItem.setShared(true); elMapItem.setOwned(false); li.setGeoId(2);
    li(1) = 5; li(2) = 6;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(6); elMapItem.setGid(5); elMapItem.setShared(true); elMapItem.setOwned(false); li.setGeoId(2);
    li(1) = 6; li(2) = 7;
    elList1d.push_back(li,elMapItem);
    
    
    //The grid1d
    grid1d->setNodes(nodes1d);
    grid1d->setElements(elList1d);
  }
  else
  {
    //Maps
    NODEMAP nodeMapItem;
    ELMAP   elMapItem;
    
    nodeMapItem.setPid(1);
    elMapItem.setPid(1);
    
    //Nodes2d
    pVect<point3d,NODEMAP> nodes2d;
    
    nodes2d.reserve(6);
    nodeMapItem.setLid(1); nodeMapItem.setGid(2);
    nodes2d.push_back(point3d(1.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(2); nodeMapItem.setGid(3);
    nodes2d.push_back(point3d(2.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(3); nodeMapItem.setGid(5);
    nodes2d.push_back(point3d(1.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(4); nodeMapItem.setGid(6);
    nodes2d.push_back(point3d(2.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(5); nodeMapItem.setGid(8);
    nodes2d.push_back(point3d(1.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(6); nodeMapItem.setGid(9);
    nodes2d.push_back(point3d(2.0, 2.0, 0.0),nodeMapItem);
    
    
    //Elements2d
    GEOELEMENT2D quad(true);
    pGraph<GEOELEMENT2D,ELMAP,NODEMAP> elList2d;
    
    elList2d.reserve(2);
    quad.setGeoId(1);
    
    elMapItem.setLid(1); elMapItem.setGid(3); elMapItem.setShared(false); elMapItem.setOwned(true);
    quad(1) = 1; quad(2) = 2; quad(3) = 4; quad(4) = 3;
    elList2d.push_back(quad,elMapItem);
    
    elMapItem.setLid(2); elMapItem.setGid(4); elMapItem.setShared(true); elMapItem.setOwned(true);
    quad(1) = 3; quad(2) = 4; quad(3) = 6; quad(4) = 5;
    elList2d.push_back(quad,elMapItem);
    
    elList2d.updateRowFinder();
    elList2d.updateColFinder();
    
    
    //The grid2d
    grid2d->setNodes(nodes2d);
    grid2d->setElements(elList2d);
    
    
    //Nodes1d
    pVect<point3d,NODEMAP> nodes1d;
    
    nodes1d.reserve(5);
    nodeMapItem.setLid(1); nodeMapItem.setGid(4);
    nodes1d.push_back(point3d(1.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(2); nodeMapItem.setGid(5);
    nodes1d.push_back(point3d(2.0, 2.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(3); nodeMapItem.setGid(6);
    nodes1d.push_back(point3d(2.0, 1.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(4); nodeMapItem.setGid(7);
    nodes1d.push_back(point3d(2.0, 0.0, 0.0),nodeMapItem);
    
    nodeMapItem.setLid(5); nodeMapItem.setGid(8);
    nodes1d.push_back(point3d(1.0, 0.0, 0.0),nodeMapItem);
    
    
    //Elements1d
    GEOELEMENT1D li(true);
    pGraph<GEOELEMENT1D,ELMAP,NODEMAP> elList1d;
    
    elList1d.reserve(4);
    
    elMapItem.setLid(1); elMapItem.setGid(4); elMapItem.setShared(true); elMapItem.setOwned(true); li.setGeoId(2);
    li(1) = 1; li(2) = 2;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(2); elMapItem.setGid(5); elMapItem.setShared(true); elMapItem.setOwned(true); li.setGeoId(2);
    li(1) = 2; li(2) = 3;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(3); elMapItem.setGid(6); elMapItem.setShared(false); elMapItem.setOwned(true); li.setGeoId(2);
    li(1) = 3; li(2) = 4;
    elList1d.push_back(li,elMapItem);
    
    elMapItem.setLid(4); elMapItem.setGid(7); elMapItem.setShared(false); elMapItem.setOwned(true); li.setGeoId(2);
    li(1) = 4; li(2) = 5;
    elList1d.push_back(li,elMapItem);
    
    
    //The grid1d
    grid1d->setNodes(nodes1d);
    grid1d->setElements(elList1d);
  }
  
  
  //Grid2d checking
  mesh2dGlobalManip<GEOSHAPE2D,ELMAP,NODEMAP> checker2d(world);
  mesh1dGlobalManip<GEOSHAPE1D,ELMAP,NODEMAP> checker1d(world);
  
  bool flag2d = checker2d.check(grid2d);
  bool flag1d = checker1d.check(grid1d);
  
  assert(flag1d);
  assert(flag2d);
  
  
  //Connecting
  connect2d<GEOSHAPE2D,ELMAP,NODEMAP> connectGrid2d(world);
  connectGrid2d.setMesh2d(grid2d);
  connectGrid2d.setMesh1d(grid1d);
  
  connectGrid2d.buildConnectivity();
  connectGrid2d.buildBoundaryConnectivity();
  
  
  //Authomatic cheking  
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(world);
  bool flagD = checkEdges.check(grid2d->getEdges());
 
  world.barrier();
  if(world.rank() == 0)
  { cout << "Check Edges : " << flagD << endl; }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 0)
  { cout << "Check Edges : " << flagD << endl; }
  sleep(1);
  
  
  world.barrier();
  if(world.rank() == 0)
  { 
    cout << "Edges---------------------------------------------------" << endl;
    cout << grid2d->getEdges() << endl << endl;
    
    cout << "Vertex to vertex----------------------------------------" << endl;    
    cout << connectGrid2d.getVertexToVertex() << endl << endl;
    
    cout << "Vertex to element---------------------------------------" << endl;
    cout << connectGrid2d.getVertexToElement() << endl << endl;
    
    cout << "Vertex to edge------------------------------------------" << endl;
    cout << connectGrid2d.getVertexToEdge() << endl << endl;
    
    cout << "Edge to element-----------------------------------------" << endl;
    cout << connectGrid2d.getEdgeToElement() << endl << endl;
    
    cout << "Element to edge-----------------------------------------" << endl;
    cout << connectGrid2d.getElementToEdge() << endl << endl;
    
    cout << "Element to element--------------------------------------" << endl;
    cout << connectGrid2d.getElementToElement() << endl << endl;
    
    cout << "Vertex is boundary--------------------------------------" << endl;
    cout << connectGrid2d.getVertexIsBoundary() << endl << endl;
    
    cout << "Element is boundary-------------------------------------" << endl;
    cout << connectGrid2d.getElementIsBoundary() << endl << endl;
    
    cout << "Edge is boundary----------------------------------------" << endl;
    cout << connectGrid2d.getEdgeIsBoundary() << endl << endl;
    
    cout << "Vertex to Bvertex---------------------------------------" << endl;
    cout << connectGrid2d.getVertexBVertex() << endl << endl;
    
    cout << "Edge to Bedge-------------------------------------------" << endl;
    cout << connectGrid2d.getEdgeBEdge() << endl << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  { 
    cout << "Edges---------------------------------------------------" << endl;
    cout << grid2d->getEdges() << endl << endl;
    
    cout << "Vertex to vertex----------------------------------------" << endl;    
    cout << connectGrid2d.getVertexToVertex() << endl << endl;
    
    cout << "Vertex to element---------------------------------------" << endl;
    cout << connectGrid2d.getVertexToElement() << endl << endl;
    
    cout << "Vertex to edge------------------------------------------" << endl;
    cout << connectGrid2d.getVertexToEdge() << endl << endl;
    
    cout << "Edge to element-----------------------------------------" << endl;
    cout << connectGrid2d.getEdgeToElement() << endl << endl;
    
    cout << "Element to edge-----------------------------------------" << endl;
    cout << connectGrid2d.getElementToEdge() << endl << endl;
    
    cout << "Element to element--------------------------------------" << endl;
    cout << connectGrid2d.getElementToElement() << endl << endl;
    
    cout << "Vertex is boundary--------------------------------------" << endl;
    cout << connectGrid2d.getVertexIsBoundary() << endl << endl;
    
    cout << "Element is boundary-------------------------------------" << endl;
    cout << connectGrid2d.getElementIsBoundary() << endl << endl;
    
    cout << "Edge is boundary----------------------------------------" << endl;
    cout << connectGrid2d.getEdgeIsBoundary() << endl << endl;
    
    cout << "Vertex to Bvertex---------------------------------------" << endl;
    cout << connectGrid2d.getVertexBVertex() << endl << endl;
    
    cout << "Edge to Bedge-------------------------------------------" << endl;
    cout << connectGrid2d.getEdgeBEdge() << endl << endl;
  }
  sleep(1);
}
