/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef OPHY3DB_HPP
#define OPHY3DB_HPP

#include "morganaFiniteElements.hpp"
#include "operatorHY3d.hpp"
#include "../morganaDofs/komplex.h"


/*! Operator hybrid 3d : (U * N) * lam - (1/2) * (cH/cV) * (t*grad(U) * t) * lam */
template<typename FIELD, typename TEST, typename COEFFH, typename COEFFV>
class opHY3dB : public operatorHY3d<FIELD,TEST>
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorHY3d<FIELD,TEST>  OPERATOR;
    
    typedef typename OPERATOR::FIELD_OUTTYPE FIELD_OUTTYPE;
    typedef typename OPERATOR::TEST_OUTTYPE  TEST_OUTTYPE;
    
    typedef typename OPERATOR::MESH3D     MESH3D;
    typedef typename OPERATOR::CONNECT3D  CONNECT3D;
    
    typedef typename COEFFH::FIELD_DOFTYPE COEFFH_DOFTYPE;
    typedef typename COEFFV::FIELD_DOFTYPE COEFFV_DOFTYPE;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool coeffOk;
    Teuchos::RCP<COEFFH> coeffH;
    Teuchos::RCP<COEFFV> coeffV;
    //@}
    
    /*! @name Functions */ //@{
  public:
    opHY3dB();
    opHY3dB(const Teuchos::RCP<communicator> & CommDev,
            const Teuchos::RCP<MESH3D> & Grid3d,
            const Teuchos::RCP<CONNECT3D> & ConnedGrid3);
    
    opHY3dB(communicator & CommDev,
            MESH3D & Grid3d,
            CONNECT3D & ConnedGrid3d);
    
    void setCoeff(const Teuchos::RCP<COEFFH> & CoeffH,
                  const Teuchos::RCP<COEFFV> & CoeffV);
    
    void eval(const point3d & Yf, sVect<komplex> & mat);
    void eval(const point3d & Yf, sVect<Real> & mat);
    //@}
};


template<typename FIELD, typename TEST, typename COEFFH, typename COEFFV>
opHY3dB<FIELD,TEST,COEFFH,COEFFV>::
opHY3dB() : operatorHY3d<FIELD,TEST>()
{
  coeffOk = false;
  
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType   == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFFH_DOFTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFFV_DOFTYPE>::myType == typeReal>::returnValue);
}

template<typename FIELD, typename TEST, typename COEFFH, typename COEFFV>
opHY3dB<FIELD,TEST,COEFFH,COEFFV>::
opHY3dB(const Teuchos::RCP<communicator> & CommDev,
        const Teuchos::RCP<MESH3D> & Grid3d,
        const Teuchos::RCP<CONNECT3D> & ConnedGrid3d) : operatorHY3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  coeffOk = false;
  
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType   == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFFH_DOFTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFFV_DOFTYPE>::myType == typeReal>::returnValue);
}

template<typename FIELD, typename TEST, typename COEFFH, typename COEFFV>
opHY3dB<FIELD,TEST,COEFFH,COEFFV>::
opHY3dB(communicator & CommDev,
        MESH3D & Grid3d,
        CONNECT3D & ConnedGrid3d) : operatorHY3d<FIELD,TEST>(CommDev,Grid3d,ConnedGrid3d)
{
  coeffOk = false;
  
  assert(staticAssert<traitsBasic<FIELD_OUTTYPE>::myType  == typePoint3d>::returnValue);
  assert(staticAssert<traitsBasic<TEST_OUTTYPE>::myType   == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFFH_DOFTYPE>::myType == typeReal>::returnValue);
  assert(staticAssert<traitsBasic<COEFFV_DOFTYPE>::myType == typeReal>::returnValue);
}

template<typename FIELD, typename TEST, typename COEFFH, typename COEFFV>
void
opHY3dB<FIELD,TEST,COEFFH,COEFFV>::
setCoeff(const Teuchos::RCP<COEFFH> & CoeffH,
         const Teuchos::RCP<COEFFV> & CoeffV)
{
  coeffOk = true;
  
  coeffH = CoeffH;
  coeffV = CoeffV;
}

template<typename FIELD, typename TEST, typename COEFFH, typename COEFFV>
void
opHY3dB<FIELD,TEST,COEFFH,COEFFV>::
eval(const point3d & Yf, sVect<Real> & mat)
{
  //Typedefs
  typedef operatorHY3d<FIELD,TEST>  OPERATOR;
  
  //Asserts
  assert(coeffOk);
  assert(mat.size() == (OPERATOR::numIndex_field() * OPERATOR::numIndex_test()) );
  
  //Geo eval
  point3d Yv = OPERATOR::mapVolumeY(Yf);
  point3d N = OPERATOR::computeNormal(Yf);
  
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
  
  coeffH->evalL(OPERATOR::el3d,Yv,cH);
  coeffV->evalL(OPERATOR::el3d,Yv,cV);  
  
  //Matrixing
  UInt tot = 1;
  
  for(UInt j=1; j <= OPERATOR::numIndex_test(); ++j)
  {
    for(UInt i=1; i <= OPERATOR::numIndex_field(); ++i)
    {
      mat(tot)  = ( point3d::dot(eval_field(i), N) * eval_test(j) );
      
      mat(tot) -= point3d::dot(gradX_field(i) * t1.getX() + 
                               gradY_field(i) * t1.getY() +
                               gradZ_field(i) * t1.getZ(), t1) * 0.5 * (cH/cV) * eval_test(j);
       
      mat(tot) -= point3d::dot(gradX_field(i) * t2.getX() + 
                               gradY_field(i) * t2.getY() +
                               gradZ_field(i) * t2.getZ(), t2) * 0.5 * (cH/cV) * eval_test(j);
    
      ++tot;
    }
  }
}

template<typename FIELD, typename TEST, typename COEFFH, typename COEFFV>
void
opHY3dB<FIELD,TEST,COEFFH,COEFFV>::
eval(const point3d & Yf, sVect<komplex> & mat)
{
  sVect<Real> matReal(mat.size());
  eval(Yf,matReal);
  
  for(UInt i=1; i <= mat.size(); ++i)
  { mat(i) = komplex(matReal(i),0.0); }
}

#endif
