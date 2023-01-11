/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPLA3DE_HPP
#define OPLA3DE_HPP

#include "morganaFiniteElements.hpp"
#include "operatorLA3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator 2d : (U * N) * lam - (cH/cV) * (t*grad(U) * t) * lam*/
template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
class opLA3dI : public operatorLA3d<FIELD,TEST>
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
    
    /*! @name Internal data */ //@{
  public:
    bool coeffOk;
    Teuchos::RCP<COEFFV> coeffV;
    Teuchos::RCP<COEFFH> coeffH;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opLA3dI();
    opLA3dI(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH3D>       & Grid3d,
            const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d,
            const Teuchos::RCP<MESH2D>       & Grid2d,
            const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d);
    
    opLA3dI(communicator & CommDev,
            MESH3D       & Grid3d,
            CONNECT3D    & ConnedGrid3d,
            MESH2D       & Grid2d,
            CONNECT2D    & ConnedGrid2d);
    
    void setCoeff(const Teuchos::RCP<COEFFV> & CoeffV,
                  const Teuchos::RCP<COEFFH> & CoeffH);
    
    void setCoeff(const COEFFV & CoeffV,
                  const COEFFH & CoeffH);
    
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};


template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
opLA3dI<FIELD,TEST,COEFFV,COEFFH>::
opLA3dI() : operatorLA3d<FIELD,TEST>()
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
opLA3dI<FIELD,TEST,COEFFV,COEFFH>::
opLA3dI(const Teuchos::RCP<communicator> & CommDev,
        const Teuchos::RCP<MESH3D>       & Grid3d,
        const Teuchos::RCP<CONNECT3D>    & ConnedGrid3d,
        const Teuchos::RCP<MESH2D>       & Grid2d,
        const Teuchos::RCP<CONNECT2D>    & ConnedGrid2d) : operatorLA3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
opLA3dI<FIELD,TEST,COEFFV,COEFFH>::
opLA3dI(communicator & CommDev,
        MESH3D       & Grid3d,
        CONNECT3D    & ConnedGrid3d,
        MESH2D       & Grid2d,
        CONNECT2D    & ConnedGrid2d) : operatorLA3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d,Grid2d,ConnedGrid2d)
{
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType  == typeReal>::returnValue);
}

template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
void
opLA3dI<FIELD,TEST,COEFFV,COEFFH>::
setCoeff(const Teuchos::RCP<COEFFV> & CoeffV,
         const Teuchos::RCP<COEFFH> & CoeffH)

{
  coeffOk = true;
  coeffV  = CoeffV;
  coeffH  = CoeffH;
}    

template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
void
opLA3dI<FIELD,TEST,COEFFV,COEFFH>::
setCoeff(const COEFFV & CoeffV,
         const COEFFH & CoeffH)
{
  coeffOk = true;
  coeffV  = Teuchos::rcpFromRef(CoeffV);
  coeffH  = Teuchos::rcpFromRef(CoeffH);
} 

template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
void
opLA3dI<FIELD,TEST,COEFFV,COEFFH>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  //Asserts
  assert(coeffOk);
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
    
  //Geo eval
  point3d Yv = OPERATOR::mapVolumeY(Yf);
  point3d N  = OPERATOR::computeNormal(Yf);
  
  std::pair<point3d,point3d> t = N.orthoBasis();
  point3d t1 = t.first;
  point3d t2 = t.second;
  
  //Field eval
  sVect<FIELD_OUTTYPE> eval_field(OPERATOR::numIndex_field());
  OPERATOR::eval_field(Yf,eval_field);
  
  //Grad eval
  sVect<FIELD_OUTTYPE> gradX_field(OPERATOR::numIndex_field());
  sVect<FIELD_OUTTYPE> gradY_field(OPERATOR::numIndex_field());
  sVect<FIELD_OUTTYPE> gradZ_field(OPERATOR::numIndex_field());
  OPERATOR::evalGrad_field(Yf,gradX_field,gradY_field,gradZ_field);  
  
  //Test eval
  sVect<TEST_OUTTYPE> eval_test(OPERATOR::numIndex_test());
  OPERATOR::eval_test(Yf,eval_test);
  
  //Coeff eval
  Real cH, cV;
  
  coeffH->evalL(OPERATOR::getEl2d(),Yf,cH);
  coeffV->evalL(OPERATOR::getEl3d(),Yv,cV);
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {
      mat(tot)  = ( point3d::dot(eval_field(i), N) * eval_test(j) );
      
      mat(tot) -= point3d::dot(gradX_field(i) * t1.getX() + 
                               gradY_field(i) * t1.getY() +
                               gradZ_field(i) * t1.getZ(), t1) * (cH/cV) * eval_test(j);
       
      mat(tot) -= point3d::dot(gradX_field(i) * t2.getX() + 
                               gradY_field(i) * t2.getY() +
                               gradZ_field(i) * t2.getZ(), t2) * (cH/cV) * eval_test(j);
    
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST, typename COEFFV, typename COEFFH>
void
opLA3dI<FIELD,TEST,COEFFV,COEFFH>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yf,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
