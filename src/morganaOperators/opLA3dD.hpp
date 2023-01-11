/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPLA3DD_HPP
#define OPLA3DD_HPP

#include "morganaFiniteElements.hpp"
#include "operatorLA3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 2d : U * Q 
 U : 2d vector
 Q : 2d vector */
template<typename FIELD, typename TEST>
class opLA3dD : public operatorLA3d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorLA3d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH2D     MESH2D;
    typedef typename OPERATOR::MESH3D     MESH3D;
    typedef typename OPERATOR::CONNECT2D  CONNECT2D;
    typedef typename OPERATOR::CONNECT3D  CONNECT3D;
    
    typedef typename TEST::GEOSHAPE   GEOSHAPE2D;
    typedef typename TEST::GEOSHAPE   TEST_GEOSHAPE;
    typedef typename TEST::PMAPTYPE   TEST_PMAPTYPE;
    
    typedef typename FIELD::GEOSHAPE  GEOSHAPE3D;
    typedef typename FIELD::GEOSHAPE  FIELD_GEOSHAPE;
    typedef typename FIELD::PMAPTYPE  FIELD_PMAPTYPE;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opLA3dD();
    opLA3dD(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH3D>       & Grid3d,
            const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d,
            const Teuchos::RCP<MESH2D>       & Grid2d,
            const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d);
    
    opLA3dD(communicator & CommDev,
            MESH3D       & Grid3d,
            CONNECT3D    & ConnedGrid3d,
            MESH2D       & Grid2d,
           CONNECT2D     & ConnedGrid2d);
    
    void eval(const point3d & Yf, sVect<komplex> & mat);
    
    void eval(const point3d & Yf,
              sVect<Real> & mat);
    //@}
};


template<typename FIELD, typename TEST>
opLA3dD<FIELD,TEST>::
opLA3dD() : operatorLA3d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint2d>::returnValue);
}

template<typename FIELD, typename TEST>
opLA3dD<FIELD,TEST>::
opLA3dD(const Teuchos::RCP<communicator> & CommDev,
        const Teuchos::RCP<MESH3D>       & Grid3d,
        const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d,
        const Teuchos::RCP<MESH2D>       & Grid2d,
        const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d) 
: operatorLA3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint2d>::returnValue);
}

template<typename FIELD, typename TEST>
opLA3dD<FIELD,TEST>::
opLA3dD(communicator & CommDev,
        MESH3D       & Grid3d,
        CONNECT3D    & ConnedGrid3d,
        MESH2D       & Grid2d, 
        CONNECT2D    & ConnedGrid2d) 
: operatorLA3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint2d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typePoint2d>::returnValue);
}

template<typename FIELD, typename TEST>
void
opLA3dD<FIELD,TEST>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  typedef operatorLA3d<FIELD,TEST>  OPERATOR;
  
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Yf,eval_field);
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Yf,eval_test);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {  
      mat(tot) = ( eval_field(i).getX() * eval_test(j).getX() +
                   eval_field(i).getY() * eval_test(j).getY() );
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST>
void
opLA3dD<FIELD,TEST>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Y,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
