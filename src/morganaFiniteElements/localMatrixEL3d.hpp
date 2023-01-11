/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef LOCALMATRIXEL3D_HPP
#define LOCALMATRIXEL3D_HPP

#include "morganaIntegrator.hpp"
#include "intTableEL3d.hpp"
#include "intSwitch.hpp"
#include "integratorsList.hpp"

#include "operatorEL3d.hpp"
#include "morganaFiniteElements.hpp"


/*! Local matrix -> Operarator for standard integration 3d */
template<typename FIELD, typename TEST, intClass INT_CLASS = intDefault, intTypes INT_TYPE = STANDARD, UInt INT_PRECISION = 1>
class localMatrixEL3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef operatorEL3d<FIELD,TEST>      OPERATOR;
    typedef typename FIELD::PMAPTYPE      PMAPTYPE;
    typedef typename OPERATOR::MESH3D     MESH3D;
    
    typedef typename FIELD::FIELD_FETYPE  FIELD_FETYPE;
    typedef typename FIELD::FIELD_DOFTYPE FIELD_DOFTYPE;
    
    typedef typename TEST::FIELD_FETYPE   TEST_FETYPE;
    typedef typename TEST::FIELD_DOFTYPE  TEST_DOFTYPE;
    //@}
    
    /*! @name Integrators Typedefs */ //@{
  public:
    static const intClass DEFAULTCLASS = intTableEL3d<FIELD_FETYPE::feBaseLabel,TEST_FETYPE::feBaseLabel>::intFlag;
    static const intClass     FINCLASS = intSwitch<INT_CLASS,DEFAULTCLASS>::intFlag;
    
    typedef typename integratorsList<FINCLASS, OPERATOR, INT_TYPE, INT_PRECISION>::INTEGRATOR INTEGRATOR;
    typedef typename INTEGRATOR::INTCARD   INTCARD;
    typedef pVect<INTCARD,PMAPTYPE>        INTCARDS;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    INTEGRATOR             integrator;
    INTCARDS               intCards;
    Teuchos::RCP<OPERATOR> op;
    bool                   opLoaded;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    localMatrixEL3d();
    localMatrixEL3d(const Teuchos::RCP<OPERATOR> & Op);
    void setOperator(const Teuchos::RCP<OPERATOR> & Op);
    void setIntCardL(const UInt & lid, const INTCARD & IntCard);
    void setIntCardG(const UInt & gid, const INTCARD & IntCard);
    void setIntCards(const INTCARDS & INTCards);
    //@}
    
    /*! @name Matrix building */ //@{
  public:    
    UInt numIndex_col();
    void indexL_col(sVect<UInt> & indices, const UInt & base = 0);
    void indexG_col(sVect<UInt> & indices, const UInt & base = 0);
    
    UInt numIndex_row();
    void indexL_row(sVect<UInt> & indices, const UInt & base = 0);
    void indexG_row(sVect<UInt> & indices, const UInt & base = 0);
    
    sVect<Real> matrix();
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
localMatrixEL3d()
{
  assert(staticAssert<INTEGRATOR::opFlag == opEL3d>::returnValue);
  opLoaded = false;
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
localMatrixEL3d(const Teuchos::RCP<OPERATOR> & Op)
{
  assert(staticAssert<INTEGRATOR::opFlag == opEL3d>::returnValue);
  opLoaded = true;
  op       = Op;
  
  //Sizing the integration cards
  Teuchos::RCP<MESH3D> grid3d = op->getIntegrationGrid();
  
  intCards.resize(grid3d->getNumElements());
  intCards.setMap(grid3d->getElements().getRowMap());
  intCards.updateFinder();
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setOperator(const Teuchos::RCP<OPERATOR> & Op)
{
  opLoaded = true;
  op       = Op;
  
  //Sizing the integration cards
  Teuchos::RCP<MESH3D> grid3d = op->getIntegrationGrid();
  
  intCards.resize(grid3d->getNumElements());
  intCards.setMap(grid3d->getElements().getRowMap());
  intCards.updateFinder();
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCardL(const UInt & lid, const INTCARD & IntCard)
{
  assert(opLoaded);
  assert(lid <= intCards.sizeL());
  intCards.getL(lid) = IntCard;
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCardG(const UInt & gid, const INTCARD & IntCard)
{
  assert(opLoaded);
  assert(intCards.isG(gid));
  intCards.getG(gid) = IntCard;
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCards(const INTCARDS & INTCards)
{
  assert(opLoaded);
  assert(INTCards.size() == intCards.size());
  intCards = INTCards;
}



//_________________________________________________________________________________________________
// MATRIX BUILD
//-------------------------------------------------------------------------------------------------
template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
UInt
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
numIndex_col()
{
  assert(opLoaded);
  return(op->numIndex_field());
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexL_col(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  
  op->indexL_field(indices);
  
  for(UInt i=1; i <= indices.size(); ++i)
  { indices(i) = indices(i) - 1 + base; }
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexG_col(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  
  op->indexG_field(indices);
  
  for(UInt i=1; i <= indices.size(); ++i)
  { indices(i) = indices(i) - 1 + base; }
}
    
template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
UInt
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
numIndex_row()
{
  assert(opLoaded);
  return(op->numIndex_test());
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexL_row(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  
  op->indexL_test(indices);
  
  for(UInt i=1; i <= indices.size(); ++i)
  { indices(i) = indices(i) - 1 + base; }
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexG_row(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  
  op->indexG_test(indices);
  
  for(UInt i=1; i <= indices.size(); ++i)
  { indices(i) = indices(i) - 1 + base; }
}

template<typename FIELD, typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
sVect<Real>
localMatrixEL3d<FIELD,TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
matrix()
{
  assert(opLoaded);
  
  UInt el = op->getElement();
  integrator.setIntCard(intCards(el));
  
  return(integrator.integration(op));
}



#endif
