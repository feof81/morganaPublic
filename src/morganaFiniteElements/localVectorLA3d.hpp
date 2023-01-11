/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef LOCALVECTORLA3D_HPP
#define LOCALVECTORLA3D_HPP

#include "morganaIntegrator.hpp"
#include "intTableFLA3d.hpp"
#include "intSwitch.hpp"
#include "integratorsList.hpp"
#include "functionalLA3d.hpp"
#include "localVectorBF3d.hpp"
#include "morganaFiniteElements.hpp"


/*! Local vector -> Operarator for the imposition of the Lagrange multipliers 3d */
template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS = intDefault, intTypes INT_TYPE = STANDARD, UInt INT_PRECISION = 1>
class localVectorLA3d
{
  /*! @name Typedefs */ //@{
  public:
    typedef functionalLA3d<TEST,GEOSHAPE3D> FUNCTIONAL;
    typedef typename TEST::PMAPTYPE         PMAPTYPE;
    typedef typename FUNCTIONAL::MESH3D     MESH3D;
    
    typedef typename TEST::FIELD_FETYPE   TEST_FETYPE;
    typedef typename TEST::FIELD_DOFTYPE  TEST_DOFTYPE;
    //@}
    
    /*! @name Integrators Typedefs */ //@{
  public:
    static const intClass DEFAULTCLASS = intTableFLA3d< TEST_FETYPE::feBaseLabel >::intFlag;
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
    localVectorLA3d();
    localVectorLA3d(const Teuchos::RCP<FUNCTIONAL> & Op);
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
template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
localVectorLA3d()
{
  assert(staticAssert<INTEGRATOR::opFlag == fnLA3d>::returnValue);
  opLoaded = false;
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
localVectorLA3d(const Teuchos::RCP<FUNCTIONAL> & Op)
{
  assert(staticAssert<INTEGRATOR::opFlag == fnLA3d>::returnValue);
  opLoaded = true;
  op       = Op;
  
  //Sizing the integration cards
  Teuchos::RCP<MESH3D> grid3d = op->getIntegrationGrid();
  
  intCards.resize(grid3d->getNumElements());
  intCards.setMap(grid3d->getElements().getRowMap());
  intCards.updateFinder();
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
setFunctional(const Teuchos::RCP<FUNCTIONAL> & Op)
{
  opLoaded = true;
  op       = Op;
  
  //Sizing the integration cards
  Teuchos::RCP<MESH3D> grid3d = op->getIntegrationGrid();
  
  intCards.resize(grid3d->getNumElements());
  intCards.setMap(grid3d->getElements().getRowMap());
  intCards.updateFinder();
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCardL(const UInt & lid, const INTCARD & IntCard)
{
  assert(opLoaded);
  assert(lid <= intCards.sizeL());
  intCards.getL(lid) = IntCard;
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCardG(const UInt & gid, const INTCARD & IntCard)
{
  assert(opLoaded);
  assert(intCards.isG(gid));
  intCards.getG(gid) = IntCard;
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
setIntCards(const INTCARDS & INTCards)
{
  assert(opLoaded);
  assert(INTCards.size() == intCards.size());
  intCards = INTCards;
}



//_________________________________________________________________________________________________
// MATRIX BUILD
//-------------------------------------------------------------------------------------------------
template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
UInt
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
numIndex_row()
{
  assert(opLoaded);
  UInt num = 0;
  
  for(UInt f=1; f <= GEOSHAPE3D::numFaces; ++f) 
  {
    op->setLocFace(f);
    
    if( op->isBoundary() )
    { num += op->numIndex_test(); }
  }
  
  return(num);
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexL_row(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  assert(indices.size() == numIndex_row());
  
  UInt k = 1;
  sVect<UInt> tempVect;
  
  for(UInt f=1; f <= GEOSHAPE3D::numFaces; ++f) 
  {
    op->setLocFace(f);
    
    if( op->isBoundary() )
    {
      tempVect.resize(op->numIndex_test());
      op->indexL_test(tempVect);
      
      for(UInt i=1; i <= op->numIndex_test(); ++i)
      {
	indices(k) = tempVect(i) -1 + base;
	++k;
      }
    }
  }
  
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
void
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
indexG_row(sVect<UInt> & indices, const UInt & base)
{
  assert(opLoaded);
  assert(indices.size() == numIndex_row());
  
  UInt k = 1;
  sVect<UInt> tempVect;
  
  for(UInt f=1; f <= GEOSHAPE3D::numFaces; ++f) 
  {
    op->setLocFace(f);
    
    if( op->isBoundary() )
    {
      tempVect.resize(op->numIndex_test());
      op->indexG_test(tempVect);
      
      for(UInt i=1; i <= op->numIndex_test(); ++i)
      {
	indices(k) = tempVect(i) -1 + base;
	++k;
      }
    }
  }
}

template<typename TEST, typename GEOSHAPE3D, intClass INT_CLASS, intTypes INT_TYPE, UInt INT_PRECISION>
sVect<Real>
localVectorLA3d<TEST,GEOSHAPE3D,INT_CLASS,INT_TYPE,INT_PRECISION>::
vector()
{
  assert(opLoaded);
  
  UInt k = 1;
  UInt el = op->getElement();
  sVect<Real> tempVect;
  sVect<Real> vect(numIndex_row());
  
  for(UInt f=1; f <= GEOSHAPE3D::numFaces; ++f) 
  {
    op->setLocFace(f);
    
    if( op->isBoundary() )
    {
      integrator.setIntCard(intCards(el));
      tempVect = integrator.integration(op);
      
      for(UInt i=1; i <= tempVect.size(); ++i)
      {
	vect(k) = tempVect(i);
	++k;
      }
    }
  }
  
  return(vect);
}

#endif
