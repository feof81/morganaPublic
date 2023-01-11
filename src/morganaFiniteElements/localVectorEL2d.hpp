/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef LOCALVECTOREL2D_HPP
#define LOCALVECTOREL2D_HPP

#include "morganaIntegrator.hpp"
#include "intTableFEL2d.hpp"
#include "intSwitch.hpp"
#include "integratorsList.hpp"

#include "functionalEL2d.hpp"
#include "morganaFiniteElements.hpp"


/*! Local vector -> Functional for standard integration 2d */
template<typename TEST, intClass INT_CLASS = intDefault, intTypes INT_TYPE = STANDARD, UInt INT_PRECISION = 1>
class localVectorEL2d
{
    /*! @name Typedefs */ //@{
  public:
    typedef functionalEL2d<TEST>          FUNCTIONAL;
    typedef typename TEST::PMAPTYPE       PMAPTYPE;
    typedef typename FUNCTIONAL::MESH2D   MESH2D;
    
    typedef typename TEST::FIELD_FETYPE   TEST_FETYPE;
    typedef typename TEST::FIELD_DOFTYPE  TEST_DOFTYPE;
    //@}
    
    /*! @name Integrators Typedefs */ //@{
  public:
    static const intClass DEFAULTCLASS = intTableFEL2d< TEST_FETYPE::feBaseLabel >::intFlag;
    static const intClass     FINCLASS = intSwitch< INT_CLASS,DEFAULTCLASS >::intFlag;
    
    typedef typename integratorsList<FINCLASS, FUNCTIONAL, INT_TYPE, INT_PRECISION>::INTEGRATOR INTEGRATOR;
    typedef typename INTEGRATOR::INTCARD   INTCARD;
    typedef pVect<INTCARD,PMAPTYPE>        INTCARDS;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    INTEGRATOR      integrator;
    INTCARDS        intCards;
    Teuchos::RCP<FUNCTIONAL> op;
    bool            opLoaded;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    localVectorEL2d();
    localVectorEL2d(const Teuchos::RCP<FUNCTIONAL> & Op);
    void setFunctional(const Teuchos::RCP<FUNCTIONAL> & Op);
    void setIntCardL(const UInt & lid, const INTCARD & IntCard);
    void setIntCardG(const UInt & gid, const INTCARD & IntCard);
    void setIntCards(const INTCARDS & INTCards);
    //@}
    
    /*! @name Matrix building */ //@{
  public:
    UInt numIndex_row();
    void indexL_row(sVect<UInt> & indices, const UInt & base = 0);
    void indexG_row(sVect<UInt> & indices, const UInt & base = 0);
    sVect<Real> vector();
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET
//-------------------------------------------------------------------------------------------------
template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
localVectorEL2d()
{
  assert(staticAssert<INTEGRATOR::opFlag == fnEL2d>::returnValue);
  opLoaded = false;
}

template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
localVectorEL2d(const Teuchos::RCP<FUNCTIONAL> & Op)
{
  assert(staticAssert<INTEGRATOR::opFlag == fnEL2d>::returnValue);
  opLoaded = true;
  op       = Op;
  
  //Sizing the integration cards
  Teuchos::RCP<MESH2D> grid2d = op->getIntegrationGrid();
  
  intCards.resize(grid2d->getNumElements());
  intCards.setMap(grid2d->getElements().getRowMap());
  intCards.updateFinder();
}

template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setFunctional(const Teuchos::RCP<FUNCTIONAL> & Op)
{
  opLoaded = true;
  op       = Op;
  
  //Sizing the integration cards
  Teuchos::RCP<MESH2D> grid2d = op->getIntegrationGrid();
  
  intCards.resize(grid2d->getNumElements());
  intCards.setMap(grid2d->getElements().getRowMap());
  intCards.updateFinder();
}

template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCardL(const UInt & lid, const INTCARD & IntCard)
{
  assert(opLoaded);
  assert(lid <= intCards.sizeL());
  intCards.getL(lid) = IntCard;
}

template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCardG(const UInt & gid, const INTCARD & IntCard)
{
  assert(opLoaded);
  assert(intCards.isG(gid));
  intCards.getG(gid) = IntCard;
}

template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCards(const INTCARDS & INTCards)
{
  assert(opLoaded);
  assert(INTCards.size() == intCards.size());
  intCards = INTCards;
}



//_________________________________________________________________________________________________
// MATRIX BUILD
//-------------------------------------------------------------------------------------------------
template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
UInt
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
numIndex_row()
{
  assert(opLoaded);
  return(op->numIndex_test());
}

template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexL_row(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  
  op->indexL_test(indices);
  
  for(UInt i=1; i <= indices.size(); ++i)
  { indices(i) = indices(i) - 1 + base; }
}

template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexG_row(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  
  op->indexG_test(indices);
  
  for(UInt i=1; i <= indices.size(); ++i)
  { indices(i) = indices(i) - 1 + base; }
}
    
template<typename TEST, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
sVect<Real>
localVectorEL2d<TEST,INT_CLASS,INT_TYPE,INT_PRECISION>::
vector()
{
  assert(opLoaded);
  
  UInt el = op->getElement();
  integrator.setIntCard(intCards(el));
  
  return(integrator.integration(op));
}


#endif
