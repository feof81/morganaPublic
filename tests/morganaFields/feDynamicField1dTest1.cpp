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

#include "pMapItem.h"
#include "pMapItemShare.h"

#include "geoShapes.h"
#include "meshInit2d.hpp"
#include "search1dA.hpp"

#include "feSpectralLH1d.hpp"
#include "feDynamicField1d.hpp"

#include "printMesh.hpp"
#include <map>

// mpirun -np 2 ./bin/morgana

using namespace std;
using namespace boost::mpi;
using namespace Teuchos;


int main(int argc, char *argv[])
{
  environment  env(argc,argv);
  communicator world;
  
  assert(world.size() == 2);
  
  UInt pid = world.rank();
  
  typedef linearLine     GEOSHAPE1D;
  typedef linearQuad     GEOSHAPE2D;
  typedef pMapItemShare  PMAPTYPE;  
  
  string meshFile  = "../tests/morganaMeshes/mignonQuad2dB.unv";
  string colorFile = "../tests/morganaMeshes/mignonQuad2dB_color.unv";
  
  meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> init(world);
  init.femap_to_stdB(meshFile, colorFile, false);
  
  
  //Download data
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>::MESH1D     MESH1D;
  typedef meshInit2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>::CONNECT1D  CONNECT1D;
  
  RCP<MESH1D>           grid1d = init.getGrid1d();
  RCP<CONNECT1D> connectGrid1d = init.getConnectGrid1d();
  
  
//Initialize the mapper
  typedef Real                    DOFTYPE;
  typedef feSpectralLH1d<PMAPTYPE> FETYPE;
  
  dofMapDynamic1d_options mapOptions;
  feDynamicField1d<FETYPE,DOFTYPE,dmd1d_vectMajor,dmd1d_standard> field;
  
  field.setCommunicator(world);
  field.setGeometry(grid1d,connectGrid1d);
  field.setOptions(mapOptions); 
  
  
  //Load the feCards
  typedef typename FETYPE::FECARD FECARD;
  FECARD FeCard;
  
  if(world.rank() == 0)
  {
    FeCard.setIsActive(true);
    FeCard.setR(1);
    field.setFeCardG(1,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(1);
    field.setFeCardG(7,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(2);
    field.setFeCardG(2,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(1);
    field.setFeCardG(3,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(2);
    field.setFeCardG(8,FeCard);
  }
  else
  {
    FeCard.setIsActive(true);
    FeCard.setR(2);
    field.setFeCardG(4,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(1);
    field.setFeCardG(9,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(1);
    field.setFeCardG(5,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(2);
    field.setFeCardG(6,FeCard);
    
    FeCard.setIsActive(true);
    FeCard.setR(2);
    field.setFeCardG(10,FeCard);
  }
  
  //Startup
  field.startup();
  
  
  //Dofs loading
  if(world.rank() == 0)
  {
    for(UInt i=1; i <= 12; ++i)
    {
      field.setDofL(i,4.0);
    }
  }
  else
  {
    for(UInt i=1; i <= 13; ++i)
    {
      field.setDofL(i,4.0);
    }
  }
  
  
  //Search 1d Init
  typedef search1dA<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>  SEARCH1D;
  typedef searchData<PMAPTYPE>                   SEARCHDATA;
  
  SEARCH1D search(world);
  search.setMesh(grid1d,connectGrid1d);
  
  search.localInit(3.0);
  search.globalInit();
  
  
  //Global search point
  pVect<point3d,PMAPTYPE>    Pvect;
  pVect<SEARCHDATA,PMAPTYPE> dataVect;
  
  point3d P,Y;
  Real V, Vx, Vy, Vz;
  
  P.setX(1.5); P.setY(0.0); P.setZ(0.0);
  Pvect.push_back(P,PMAPTYPE(1,1,pid,false,false));
  Pvect.updateFinder();
  
  search.findGlobal(Pvect,dataVect);
  
  UInt elL = dataVect(1).getElMap().getLid();
  UInt elG = dataVect(1).getElMap().getGid();
         Y = dataVect(1).getLocCoord();
  
  
  //Evaluate
  if(pid == dataVect(1).getElMap().getPid())
  {
    cout << "elL " << elL << " elG " << elG << endl;
    field.evalL(elL,Y,V);
    field.gradG(elG,Y,Vx,Vy,Vz);
    cout << "pid " << pid << endl << " eval " << V << " gradX " << Vx << " gradY " << Vy << " gradZ " << Vz << endl;
  }
}
