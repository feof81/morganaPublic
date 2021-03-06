/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FIELDSALGEBRA_HPP
#define FIELDSALGEBRA_HPP

#include "typesInterface.hpp"


/*! Basic algebraic operations for the fields */
template<typename FIELD>
class fieldsAlgebra
{
    /*! @name Typedefs */ //@{
  public:
    typedef Teuchos::RCP<FIELD>           FIELD_RCP;
    typedef Teuchos::RCP<const FIELD>     CONSTFIELD_RCP;
    typedef typename FIELD::FIELD_DOFTYPE DOFTYPE;
    //@}
  
    /*! @name Combination functions */ //@{
  public:
    static FIELD_RCP linearCombination(const CONSTFIELD_RCP & fieldA,
				       const Real           & coeffA,
				       const CONSTFIELD_RCP & fieldB,
				       const Real           & coeffB );
    
    static FIELD     linearCombination(const FIELD & fieldA,
		                       const Real  & coeffA,
				       const FIELD & fieldB,
				       const Real  & coeffB );
    
    static FIELD_RCP linearCombination(const CONSTFIELD_RCP & fieldA,
				       const Real           & coeffA,
				       const CONSTFIELD_RCP & fieldB,
				       const Real           & coeffB,
				       const CONSTFIELD_RCP & fieldC,
				       const Real           & coeffC );
    
    static FIELD     linearCombination(const FIELD & fieldA,
				       const Real  & coeffA,
				       const FIELD & fieldB,
				       const Real  & coeffB,
				       const FIELD & fieldC,
				       const Real  & coeffC );
    
    static FIELD_RCP productCombination(const CONSTFIELD_RCP & fieldA1,
					 const CONSTFIELD_RCP & fieldA2,
					 const Real           & coeffA );
    
    static FIELD     productCombination(const FIELD & fieldA1,
					 const FIELD & fieldA2,
					 const Real  & coeffA );
    
    static FIELD_RCP productCombination(const CONSTFIELD_RCP & fieldA1,
					 const CONSTFIELD_RCP & fieldA2,
					 const Real           & coeffA,
                                        const CONSTFIELD_RCP & fieldB1,
					 const CONSTFIELD_RCP & fieldB2,
					 const Real           & coeffB );
    
    static FIELD     productCombination(const FIELD & fieldA1,
					 const FIELD & fieldA2,
					 const Real  & coeffA,
                                        const FIELD & fieldB1,
					 const FIELD & fieldB2,
					 const Real  & coeffB );
    
    static FIELD_RCP productCombination(const CONSTFIELD_RCP & fieldA1,
					 const CONSTFIELD_RCP & fieldA2,
					 const Real           & coeffA,
                                        const CONSTFIELD_RCP & fieldB1,
					 const CONSTFIELD_RCP & fieldB2,
					 const Real           & coeffB,
                                        const CONSTFIELD_RCP & fieldC1,
					 const CONSTFIELD_RCP & fieldC2,
					 const Real & coeffC );
    
    static FIELD     productCombination(const FIELD & fieldA1,
					 const FIELD & fieldA2,
					 const Real  & coeffA,
                                        const FIELD & fieldB1,
					 const FIELD & fieldB2,
					 const Real  & coeffB,
                                        const FIELD & fieldC1,
					 const FIELD & fieldC2,
					 const Real  & coeffC );
    //@}
    
    /*! @name Modification functions */ //@{
  public:
    static void  scale(FIELD_RCP & field, const Real & coeff);
    static void  scale(FIELD & field, const Real & coeff);
    static void  sumDof(FIELD_RCP & field, const DOFTYPE & dof);
    static void  sumDof(FIELD & field, const DOFTYPE & dof);    
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename FIELD>
typename fieldsAlgebra<FIELD>::FIELD_RCP
fieldsAlgebra<FIELD>::
linearCombination(const CONSTFIELD_RCP & fieldA, const Real & coeffA, const CONSTFIELD_RCP & fieldB, const Real & coeffB)
{
  FIELD outField = linearCombination(*fieldA,coeffA,*fieldB,coeffB);
  return(Teuchos::rcp(new FIELD(outField)));
}

template<typename FIELD>
FIELD
fieldsAlgebra<FIELD>::
linearCombination(const FIELD & fieldA, const Real & coeffA, const FIELD & fieldB, const Real & coeffB)
{
  typedef typename FIELD::DOFVECT DOFVECT;
  
  DOFVECT dofVectA = fieldA.getDofVect();
  DOFVECT dofVectB = fieldB.getDofVect();
  
  assert(dofVectA.size() == dofVectB.size());
  
  for(UInt i=1; i <= dofVectA.size(); ++i)
  { dofVectA(i) = dofVectA(i) * coeffA + dofVectB(i) * coeffB; }
  
  FIELD outField(fieldA);
  outField.setDofVect(dofVectA);
  
  return(outField);
}

template<typename FIELD>
typename fieldsAlgebra<FIELD>::FIELD_RCP
fieldsAlgebra<FIELD>::
linearCombination(const CONSTFIELD_RCP & fieldA, const Real & coeffA, const CONSTFIELD_RCP & fieldB, const Real & coeffB, const CONSTFIELD_RCP & fieldC, const Real & coeffC)
{
  FIELD outField = linearCombination(*fieldA,coeffA,*fieldB,coeffB,*fieldC,coeffC);
  return(Teuchos::rcp(new FIELD(outField)));
}

template<typename FIELD>
FIELD
fieldsAlgebra<FIELD>::
linearCombination(const FIELD & fieldA, const Real & coeffA, const FIELD & fieldB, const Real & coeffB, const FIELD & fieldC, const Real & coeffC)
{
  typedef typename FIELD::DOFVECT DOFVECT;
  
  DOFVECT dofVectA = fieldA.getDofVect();
  DOFVECT dofVectB = fieldB.getDofVect();
  DOFVECT dofVectC = fieldC.getDofVect();
  
  assert(dofVectA.size() == dofVectB.size());
  assert(dofVectA.size() == dofVectC.size());
  
  for(UInt i=1; i <= dofVectA.size(); ++i)
  { 
    dofVectA(i) = dofVectA(i) * coeffA +
                  dofVectB(i) * coeffB +
                  dofVectC(i) * coeffC;
  }
  
  FIELD outField(fieldA);
  outField.setDofVect(dofVectA);
  
  return(outField);
}

template<typename FIELD>
typename fieldsAlgebra<FIELD>::FIELD_RCP
fieldsAlgebra<FIELD>::
productCombination(const CONSTFIELD_RCP & fieldA1,
		   const CONSTFIELD_RCP & fieldA2,
		   const Real           & coeffA )
{
  FIELD outField = productCombination(*fieldA1,*fieldA2,coeffA);
  return(Teuchos::rcp(new FIELD(outField)));
}

template<typename FIELD>
FIELD
fieldsAlgebra<FIELD>::
productCombination(const FIELD & fieldA1,
		   const FIELD & fieldA2,
		   const Real  & coeffA )
{
  typedef typename FIELD::DOFVECT DOFVECT;
  
  DOFVECT dofVectA1 = fieldA1.getDofVect();
  DOFVECT dofVectA2 = fieldA2.getDofVect();
  
  assert(dofVectA1.size() == dofVectA2.size());
  
  for(UInt i=1; i <= dofVectA1.size(); ++i)
  { 
    dofVectA1(i) = dofVectA1(i) * dofVectA2(i) * coeffA;
  }
  
  FIELD outField(fieldA1);
  outField.setDofVect(dofVectA1);
  
  return(outField);
}

template<typename FIELD>
typename fieldsAlgebra<FIELD>::FIELD_RCP
fieldsAlgebra<FIELD>::
productCombination(const CONSTFIELD_RCP & fieldA1,
		   const CONSTFIELD_RCP & fieldA2,
		   const Real           & coeffA,
                   const CONSTFIELD_RCP & fieldB1,
		   const CONSTFIELD_RCP & fieldB2,
		   const Real           & coeffB )
{
  FIELD outField = productCombination(*fieldA1,*fieldA2,coeffA, *fieldB1,*fieldB2,coeffB);
  return(Teuchos::rcp(new FIELD(outField)));
}

template<typename FIELD>
FIELD
fieldsAlgebra<FIELD>::
productCombination(const FIELD & fieldA1,
		   const FIELD & fieldA2,
		   const Real  & coeffA,
                   const FIELD & fieldB1,
		   const FIELD & fieldB2,
		   const Real  & coeffB )
{
  typedef typename FIELD::DOFVECT DOFVECT;
  
  DOFVECT dofVectA1 = fieldA1.getDofVect();
  DOFVECT dofVectA2 = fieldA2.getDofVect();
  
  DOFVECT dofVectB1 = fieldB1.getDofVect();
  DOFVECT dofVectB2 = fieldB2.getDofVect();
  
  assert(dofVectA1.size() == dofVectA2.size());
  assert(dofVectA1.size() == dofVectB1.size());
  assert(dofVectA1.size() == dofVectB2.size());
  
  for(UInt i=1; i <= dofVectA1.size(); ++i)
  { 
    dofVectA1(i) = dofVectA1(i) * dofVectA2(i) * coeffA +
                   dofVectB1(i) * dofVectB2(i) * coeffB;
  }
  
  FIELD outField(fieldA1);
  outField.setDofVect(dofVectA1);
  
  return(outField);
}

template<typename FIELD>
typename fieldsAlgebra<FIELD>::FIELD_RCP
fieldsAlgebra<FIELD>::
productCombination(const CONSTFIELD_RCP & fieldA1,
		   const CONSTFIELD_RCP & fieldA2,
		   const Real           & coeffA,
                   const CONSTFIELD_RCP & fieldB1,
		   const CONSTFIELD_RCP & fieldB2,
		   const Real           & coeffB,
                   const CONSTFIELD_RCP & fieldC1,
		   const CONSTFIELD_RCP & fieldC2,
		   const Real & coeffC )
{
  FIELD outField = productCombination(*fieldA1,*fieldA2,coeffA, *fieldB1,*fieldB2,coeffB, *fieldC1,*fieldC2,coeffC);
  return(Teuchos::rcp(new FIELD(outField)));
}

template<typename FIELD>
FIELD
fieldsAlgebra<FIELD>::
productCombination(const FIELD & fieldA1,
		   const FIELD & fieldA2,
		   const Real  & coeffA,
                   const FIELD & fieldB1,
		   const FIELD & fieldB2,
		   const Real  & coeffB,
                   const FIELD & fieldC1,
		   const FIELD & fieldC2,
		   const Real  & coeffC )
{
  typedef typename FIELD::DOFVECT DOFVECT;
  
  DOFVECT dofVectA1 = fieldA1.getDofVect();
  DOFVECT dofVectA2 = fieldA2.getDofVect();
  
  DOFVECT dofVectB1 = fieldB1.getDofVect();
  DOFVECT dofVectB2 = fieldB2.getDofVect();
  
  DOFVECT dofVectC1 = fieldC1.getDofVect();
  DOFVECT dofVectC2 = fieldC2.getDofVect();
  
  assert(dofVectA1.size() == dofVectA2.size());
  assert(dofVectA1.size() == dofVectB1.size());
  assert(dofVectA1.size() == dofVectB2.size());
  assert(dofVectA1.size() == dofVectC1.size());
  assert(dofVectA1.size() == dofVectC2.size());
  
  for(UInt i=1; i <= dofVectA1.size(); ++i)
  { 
    dofVectA1(i) = dofVectA1(i) * dofVectA2(i) * coeffA +
                   dofVectB1(i) * dofVectB2(i) * coeffB +
                   dofVectC1(i) * dofVectC2(i) * coeffC;
  }
  
  FIELD outField(fieldA1);
  outField.setDofVect(dofVectA1);
  
  return(outField);
}

template<typename FIELD>
void
fieldsAlgebra<FIELD>::
scale(FIELD_RCP & field, const Real & coeff)
{
  scale(*field,coeff);
}

template<typename FIELD>
void
fieldsAlgebra<FIELD>::
scale(FIELD & field, const Real & coeff)
{
  typedef typename FIELD::DOFVECT DOFVECT;
  
  DOFVECT dofVect = field.getDofVect();
  
  for(UInt i=1; i <= dofVect.size(); ++i)
  { dofVect(i) = dofVect(i) * coeff; }
  
  field.setDofVect(dofVect);
}

template<typename FIELD>
void
fieldsAlgebra<FIELD>::
sumDof(FIELD_RCP & field, const DOFTYPE & dof)
{
  sumDof(*field,dof);
}

template<typename FIELD>
void
fieldsAlgebra<FIELD>::
sumDof(FIELD & field, const DOFTYPE & dof)
{
  typedef typename FIELD::DOFVECT DOFVECT;
  
  DOFVECT dofVect = field.getDofVect();
  
  for(UInt i=1; i <= dofVect.size(); ++i)
  { dofVect(i) = dofVect(i) + dof; }
  
  field.setDofVect(dofVect);
}


#endif
