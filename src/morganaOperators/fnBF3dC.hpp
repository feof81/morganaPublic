/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FNBF3DC_HPP
#define FNBF3DC_HPP
#include "../morganaFiniteElements/functionalBF3d.hpp"
#include "../morganaDofs/morganaTypes.hpp"
#include "../morganaDofs/point3d.h"
#include "../morganaDofs/komplex.h"

/*! Functional 3d: (V * N) */
template<typename TEST>
class fnBF3dC : public functionalBF3d<TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalBF3d<TEST>  FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH2D        MESH2D;
    typedef typename FUNCTIONAL::MESH3D        MESH3D;
    typedef typename FUNCTIONAL::CONNECT2D     CONNECT2D;
    typedef typename FUNCTIONAL::CONNECT3D     CONNECT3D;
    
    typedef typename TEST::PMAPTYPE   TEST_PMAPTYPE;
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    set<UInt>            activeIds;
    Teuchos::RCP<MESH2D> grid2d;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fnBF3dC();
    fnBF3dC(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    fnBF3dC(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    void setGeoIds(const set<UInt> & ActiveIds);
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST>
fnBF3dC<TEST>::
fnBF3dC() : functionalBF3d<TEST>()
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
}

template<typename TEST>
fnBF3dC<TEST>::
fnBF3dC(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d) :
functionalBF3d<TEST>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  
  grid2d = Grid2d;
}

template<typename TEST>
fnBF3dC<TEST>::
fnBF3dC(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d) : 
functionalBF3d<TEST>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint3d>::returnValue);
  
  grid2d = Teuchos::rcpFromRef(Grid2d);
}

template<typename TEST>
void
fnBF3dC<TEST>::
setGeoIds(const set<UInt> & ActiveIds)
{
  activeIds = ActiveIds;
}

template<typename TEST>
void
fnBF3dC<TEST>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  typedef functionalBF3d<TEST> FUNCTIONAL;
  
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Yf, eval_field);
  
  //Normal evaluation
  point3d N = FUNCTIONAL::computeNormal(Yf);
  
  //Element active
  UInt el2d  = FUNCTIONAL::getEl2d();
  UInt geoId = grid2d->getElementL(el2d).getGeoId();
  Real val   = Real(activeIds.count(geoId) == 1);

  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  {
    mat(j) = point3d::dot(eval_field(j), N) * val;
  }
}

template<typename TEST>
void
fnBF3dC<TEST>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yf,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
