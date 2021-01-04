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
#include "mesh3d.hpp"
#include "mesh3dGlobalManip.hpp"
#include "connect3d.hpp"

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;

/*! This test was constructed to test the correctness of the performance upgrades to connect3d */
int main(int argc, char *argv[])
{
  typedef linearTetra            GEOSHAPE3D;
  typedef linearTriangle         GEOSHAPE2D;
  typedef linearLine             GEOSHAPE1D;
  typedef geoElement<GEOSHAPE3D> GEOELEMENT3D;
  typedef geoElement<GEOSHAPE2D> GEOELEMENT2D;
  typedef geoElement<GEOSHAPE1D> GEOELEMENT1D;
  typedef pMapItemShare          ELMAP;
  typedef pMapItemShare          NODEMAP;
  typedef mesh2d<GEOSHAPE2D,ELMAP,NODEMAP> MESH2D;
  typedef mesh3d<GEOSHAPE3D,ELMAP,NODEMAP> MESH3D;
  
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  Teuchos::RCP<MESH2D> grid2d = Teuchos::rcp(new MESH2D);
  Teuchos::RCP<MESH3D> grid3d = Teuchos::rcp(new MESH3D);
  
  if(world.rank() == 0)
  {
    //Nodes 3d
    pVect<point3d,NODEMAP> nodes3d;
    
    nodes3d.reserve(4);
    nodes3d.push_back(point3d(0.0, 0.0, 0.0),pMapItemShare(1,1,0,true,true));
    nodes3d.push_back(point3d(1.0, 0.0, 0.0),pMapItemShare(2,2,0,false,true));
    nodes3d.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(3,3,0,true,true));
    nodes3d.push_back(point3d(0.5, 0.5, 1.0),pMapItemShare(4,5,0,true,true));
    nodes3d.updateFinder();
    
    
    
    //Elements 3d
    GEOELEMENT3D tet(true);
    pGraph<GEOELEMENT3D,ELMAP,NODEMAP> elList3d;
    
    elList3d.reserve(1);
    tet.setGeoId(1); tet(1) = 1; tet(2) = 2; tet(3) = 3; tet(4) = 4;
    elList3d.push_back(tet, pMapItemShare(1,1,false,true)); elList3d.getRowMapL(1).setPid(0);
    
    //The grid3d
    grid3d->setNodes(nodes3d);
    grid3d->setElements(elList3d);
    
    //Nodes 2d
    pVect<point3d,NODEMAP> nodes2d;
    
    nodes2d.reserve(4);
    nodes2d.push_back(point3d(0.0, 0.0, 0.0),pMapItemShare(1,5,0,true,true));
    nodes2d.push_back(point3d(1.0, 0.0, 0.0),pMapItemShare(2,4,0,false,true));
    nodes2d.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(3,3,0,true,true));
    nodes2d.push_back(point3d(0.5, 0.5, 1.0),pMapItemShare(4,1,0,true,true));
    nodes2d.updateFinder();
    
    //Elements 2d
    GEOELEMENT2D tri(true);
    pGraph<GEOELEMENT2D,ELMAP,NODEMAP> elList2d;
    
    elList2d.reserve(3);
    tri.setGeoId(1); tri(1) = 1; tri(2) = 2; tri(3) = 4;
    elList2d.push_back(tri, pMapItemShare(1,1,false,true));
    tri.setGeoId(1); tri(1) = 2; tri(2) = 3; tri(3) = 4;
    elList2d.push_back(tri, pMapItemShare(2,2,false,true));
    tri.setGeoId(1); tri(1) = 1; tri(2) = 3; tri(3) = 2;
    elList2d.push_back(tri, pMapItemShare(3,5,false,true));
    
    //Edges2d
    GEOELEMENT1D ed(true);
    pGraph<GEOELEMENT1D,ELMAP,NODEMAP> elList1d;
    
    elList1d.reserve(6);
    ed.setGeoId(1); ed(1) = 1; ed(2) = 2;
    elList1d.push_back(ed, pMapItemShare(1,1,false,true));
    ed.setGeoId(1); ed(1) = 2; ed(2) = 3;
    elList1d.push_back(ed, pMapItemShare(2,5,false,true));
    ed.setGeoId(1); ed(1) = 1; ed(2) = 3;
    elList1d.push_back(ed, pMapItemShare(3,2,true,false));
    ed.setGeoId(1); ed(1) = 1; ed(2) = 4;
    elList1d.push_back(ed, pMapItemShare(4,4,true,false));
    ed.setGeoId(1); ed(1) = 2; ed(2) = 4;
    elList1d.push_back(ed, pMapItemShare(5,6,false,true));
    ed.setGeoId(1); ed(1) = 3; ed(2) = 4;
    elList1d.push_back(ed, pMapItemShare(6,8,true,true));
    
    //The grid2d
    grid2d->setNodes(nodes2d);
    grid2d->setElements(elList2d);
    grid2d->setEdges(elList1d);
  }
  else
  {
    //Nodes 3d
    pVect<point3d,NODEMAP> nodes3d;
    
    nodes3d.reserve(4);
    nodes3d.push_back(point3d(0.0, 0.0, 0.0),pMapItemShare(1,1,1,true,false));
    nodes3d.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(2,3,1,true,false));
    nodes3d.push_back(point3d(0.0, 1.0, 0.0),pMapItemShare(3,4,1,false,true));
    nodes3d.push_back(point3d(0.5, 0.5, 1.0),pMapItemShare(4,5,1,true,false));
    nodes3d.updateFinder();
    
    //Elements 3d
    GEOELEMENT3D tet(true);
    pGraph<GEOELEMENT3D,ELMAP,NODEMAP> elList3d;
    
    elList3d.reserve(1);
    tet.setGeoId(2); tet(1) = 1; tet(2) = 2; tet(3) = 3; tet(4) = 4;
    elList3d.push_back(tet, pMapItemShare(1,2,false,true)); elList3d.getRowMapL(1).setPid(1);
    
    //The grid3d
    grid3d->setNodes(nodes3d);
    grid3d->setElements(elList3d);
    
    //Nodes 2d
    pVect<point3d,NODEMAP> nodes2d;
    
    nodes2d.reserve(4);
    nodes2d.push_back(point3d(0.0, 0.0, 0.0),pMapItemShare(1,5,1,true,false));
    nodes2d.push_back(point3d(1.0, 1.0, 0.0),pMapItemShare(2,3,1,true,false));
    nodes2d.push_back(point3d(0.0, 1.0, 0.0),pMapItemShare(3,2,1,false,true));
    nodes2d.push_back(point3d(0.5, 0.5, 1.0),pMapItemShare(4,1,1,true,false));
    nodes2d.updateFinder();
    
    //Elements 2d
    GEOELEMENT2D tri(true);
    pGraph<GEOELEMENT2D,ELMAP,NODEMAP> elList2d;
    
    elList2d.reserve(3);
    
    tri.setGeoId(2); tri(1) = 2; tri(2) = 3; tri(3) = 4;
    elList2d.push_back(tri, pMapItemShare(1,3,false,true));
    tri.setGeoId(2); tri(1) = 3; tri(2) = 1; tri(3) = 4;
    elList2d.push_back(tri, pMapItemShare(2,4,false,true));
    tri.setGeoId(2); tri(1) = 2; tri(2) = 1; tri(3) = 3;
    elList2d.push_back(tri, pMapItemShare(3,6,false,true));

    //Edges2d
    GEOELEMENT1D ed(true);
    pGraph<GEOELEMENT1D,ELMAP,NODEMAP> elList1d;
    
    elList1d.reserve(6);
    ed.setGeoId(1); ed(1) = 1; ed(2) = 2;
    elList1d.push_back(ed, pMapItemShare(1,2,true,true));
    ed.setGeoId(1); ed(1) = 2; ed(2) = 3;
    elList1d.push_back(ed, pMapItemShare(2,7,false,true));
    ed.setGeoId(1); ed(1) = 1; ed(2) = 3;
    elList1d.push_back(ed, pMapItemShare(3,3,false,true));
    ed.setGeoId(1); ed(1) = 1; ed(2) = 4;
    elList1d.push_back(ed, pMapItemShare(4,4,true,true));
    ed.setGeoId(1); ed(1) = 2; ed(2) = 4;
    elList1d.push_back(ed, pMapItemShare(5,8,true,false));
    ed.setGeoId(1); ed(1) = 3; ed(2) = 4;
    elList1d.push_back(ed, pMapItemShare(6,9,false,true));
    
    //The grid2d
    grid2d->setNodes(nodes2d);
    grid2d->setElements(elList2d);
    grid2d->setEdges(elList1d);
  }
  
  //Map transfer
  grid2d->transferMap();
  grid3d->transferMap();
  
  //Mesh cheking
  mesh3dGlobalManip<GEOSHAPE3D,ELMAP,NODEMAP> checker(world);
  bool flag = checker.check(grid3d);
  assert(flag);
  
  //Mesh connecting
  connect3d<GEOSHAPE3D,ELMAP,NODEMAP> gridConnect(world);
  gridConnect.setMesh3d(grid3d);
  gridConnect.setMesh2d(grid2d);
  
  gridConnect.buildConnectivity();
  gridConnect.buildBoundaryConnectivity();
  
  //Authomatic cheking
  pGraphGlobalManip<GEOELEMENT2D,ELMAP,NODEMAP> checkFaces(world);
  bool flagF = checkFaces.check(grid3d->getFaces());
  
  pGraphGlobalManip<GEOELEMENT1D,ELMAP,NODEMAP> checkEdges(world);
  bool flagD = checkEdges.check(grid3d->getEdges());
 
  
  //Printout
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Check Faces : " << flagF << endl;
    cout << "Check Edges : " << flagD << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Check Faces : " << flagF << endl;
    cout << "Check Edges : " << flagD << endl;
  }
  sleep(1);
  
  
  //Printout
  world.barrier();
  if(world.rank() == 0)
  {
    cout << "Grid3d--------------------------------------------------" << endl;
    cout << *grid3d << endl << endl;
    
    cout << "Faces---------------------------------------------------" << endl;
    cout << grid3d->getFaces() << endl << endl;
 
    cout << "Edges---------------------------------------------------" << endl;
    cout << grid3d->getEdges() << endl << endl;
    
    cout << "Vertex to vertex----------------------------------------" << endl;    
    cout << gridConnect.getVertexToVertex() << endl << endl;
    
    cout << "Vertex to element---------------------------------------" << endl;
    cout << gridConnect.getVertexToElement() << endl << endl;
    
    cout << "Vertex to edge------------------------------------------" << endl;
    cout << gridConnect.getVertexToEdge() << endl << endl;
    
    cout << "Face to element-----------------------------------------" << endl;
    cout << gridConnect.getFaceToElement() << endl << endl;
    
    cout << "Element to edge-----------------------------------------" << endl;
    cout << gridConnect.getElementToEdge() << endl << endl;
    
    cout << "Element to face-----------------------------------------" << endl;
    cout << gridConnect.getElementToFace() << endl << endl;
    
    cout << "Element to element--------------------------------------" << endl;
    cout << gridConnect.getElementToElement() << endl << endl;
 
    
    cout << "Vertex is boundary--------------------------------------" << endl;
    cout << gridConnect.getVertexIsBoundary() << endl << endl;
    
    cout << "Element is boundary-------------------------------------" << endl;
    cout << gridConnect.getElementIsBoundary() << endl << endl;
    
    cout << "Face is boundary----------------------------------------" << endl;
    cout << gridConnect.getFaceIsBoundary() << endl << endl;
    
    cout << "Edge is boundary----------------------------------------" << endl;
    cout << gridConnect.getEdgeIsBoundary() << endl << endl;
    
    cout << "Vertex to Bvertex---------------------------------------" << endl;
    cout << gridConnect.getVertexBVertex() << endl << endl;
    
    cout << "Face to Bface-------------------------------------------" << endl;
    cout << gridConnect.getFaceBFace() << endl << endl;
    
    cout << "Edge to Bedge-------------------------------------------" << endl;
    cout << gridConnect.getEdgeBEdge() << endl << endl;
  }
  sleep(1);
  
  world.barrier();
  if(world.rank() == 1)
  {
    cout << "Grid3d--------------------------------------------------" << endl;
    cout << *grid3d << endl << endl;
    
    cout << "Faces---------------------------------------------------" << endl;
    cout << grid3d->getFaces() << endl << endl;
 
    cout << "Edges---------------------------------------------------" << endl;
    cout << grid3d->getEdges() << endl << endl;
    
    cout << "Vertex to vertex----------------------------------------" << endl;    
    cout << gridConnect.getVertexToVertex() << endl << endl;
    
    cout << "Vertex to element---------------------------------------" << endl;
    cout << gridConnect.getVertexToElement() << endl << endl;
    
    cout << "Vertex to edge------------------------------------------" << endl;
    cout << gridConnect.getVertexToEdge() << endl << endl;
    
    cout << "Face to element-----------------------------------------" << endl;
    cout << gridConnect.getFaceToElement() << endl << endl;
    
    cout << "Element to edge-----------------------------------------" << endl;
    cout << gridConnect.getElementToEdge() << endl << endl;
    
    cout << "Element to face-----------------------------------------" << endl;
    cout << gridConnect.getElementToFace() << endl << endl;
    
    cout << "Element to element--------------------------------------" << endl;
    cout << gridConnect.getElementToElement() << endl << endl;
 
    
    cout << "Vertex is boundary--------------------------------------" << endl;
    cout << gridConnect.getVertexIsBoundary() << endl << endl;
    
    cout << "Element is boundary-------------------------------------" << endl;
    cout << gridConnect.getElementIsBoundary() << endl << endl;
    
    cout << "Face is boundary----------------------------------------" << endl;
    cout << gridConnect.getFaceIsBoundary() << endl << endl;
    
    cout << "Edge is boundary----------------------------------------" << endl;
    cout << gridConnect.getEdgeIsBoundary() << endl << endl;
    
    cout << "Vertex to Bvertex---------------------------------------" << endl;
    cout << gridConnect.getVertexBVertex() << endl << endl;
    
    cout << "Face to Bface-------------------------------------------" << endl;
    cout << gridConnect.getFaceBFace() << endl << endl;
    
    cout << "Edge to Bedge-------------------------------------------" << endl;
    cout << gridConnect.getEdgeBEdge() << endl << endl;
  }
  sleep(1);
}
