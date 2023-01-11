/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FNLA3DC_HPP
#define FNLA3DC_HPP
#include "../morganaFiniteElements/functionalLA3d.hpp"
#include "../morganaDofs/morganaTypes.hpp"
#include "../morganaDofs/point3d.h"
#include "../morganaDofs/komplex.h"


/*! Functional: - f * v */
template<typename TEST, typename GEOSHAPE3D>
class fnLA3dC : public functionalLA3d<TEST,GEOSHAPE3D>
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalLA3d<TEST,GEOSHAPE3D>    FUNCTIONAL;
    
    typedef typename FUNCTIONAL::TEST_OUTTYPE  TEST_OUTTYPE;
    typedef typename FUNCTIONAL::MESH2D        MESH2D;
    typedef typename FUNCTIONAL::MESH3D        MESH3D;
    typedef typename FUNCTIONAL::CONNECT2D     CONNECT2D;
    typedef typename FUNCTIONAL::CONNECT3D     CONNECT3D;
    
    typedef typename TEST::PMAPTYPE TEST_PMAPTYPE;
    typedef typename TEST::GEOSHAPE TEST_GEOSHAPE;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    set<UInt>            activeIds;
    Teuchos::RCP<MESH2D> grid2d;
    //@}
    
    /*! @name Functions */ //@{
  public:
    fnLA3dC();
    fnLA3dC(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d);
    fnLA3dC(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d);
    void setGeoIds(const set<UInt> & ActiveIds);
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE3D>
fnLA3dC<TEST,GEOSHAPE3D>::
fnLA3dC() : functionalLA3d<TEST,GEOSHAPE3D>()
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 2>::returnValue);  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename TEST, typename GEOSHAPE3D>
fnLA3dC<TEST,GEOSHAPE3D>::
fnLA3dC(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<MESH3D> & Grid3d, const Teuchos::RCP<CONNECT3D> & ConnectGrid3d, const Teuchos::RCP<MESH2D> & Grid2d, const Teuchos::RCP<CONNECT2D> & ConnectGrid2d) :
functionalLA3d<TEST,GEOSHAPE3D>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 2>::returnValue);  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  
  grid2d  = Grid2d;
}

template<typename TEST, typename GEOSHAPE3D>
fnLA3dC<TEST,GEOSHAPE3D>::
fnLA3dC(communicator & CommDev, MESH3D & Grid3d, CONNECT3D & ConnectGrid3d, MESH2D & Grid2d, CONNECT2D & ConnectGrid2d) :
functionalLA3d<TEST,GEOSHAPE3D>(CommDev,Grid3d,ConnectGrid3d,Grid2d,ConnectGrid2d)
{
  assert(staticAssert<TEST_GEOSHAPE::nDim == 2>::returnValue);  
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
  
  grid2d  = Grid2d;
}

template<typename TEST, typename GEOSHAPE3D>
void
fnLA3dC<TEST,GEOSHAPE3D>::
setGeoIds(const set<UInt> & ActiveIds)
{
  activeIds = ActiveIds;
}

template<typename TEST, typename GEOSHAPE3D>
void
fnLA3dC<TEST,GEOSHAPE3D>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  typedef functionalLA3d<TEST,GEOSHAPE3D> FUNCTIONAL;
  
  assert(mat.size() == FUNCTIONAL::numIndex_test() );
  
  //Field eval
  sVect<TEST_OUTTYPE> eval_field(FUNCTIONAL::numIndex_test());
  FUNCTIONAL::eval_test(Yf, eval_field);
  
  //Element active
  UInt el2d  = FUNCTIONAL::getEl2d();
  UInt geoId = grid2d->getElementL(el2d).getGeoId();
  Real val   = Real(activeIds.count(geoId) == 1);
  
  //Vectoring
  for(UInt j=1; j <= FUNCTIONAL::numIndex_test(); ++j)
  {
    mat(j) = - eval_field(j) * val;
  }
}

template<typename TEST, typename GEOSHAPE3D>
void
fnLA3dC<TEST,GEOSHAPE3D>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yf,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}


#endif
