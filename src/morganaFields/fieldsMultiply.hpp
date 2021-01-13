/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FIELDSMULTIPLY_HPP
#define FIELDSMULTIPLY_HPP

#include "traitsBasic.h"
#include "traitsMultiply.h"


//_________________________________________________________________________________________________
// TRAITS CLASS
//-------------------------------------------------------------------------------------------------
enum fieldsMultiplyOp {fmMultiply};

/*! Trait class for binary operations */
template<typename DOFTYPEA, typename DOFTYPEB, fieldsMultiplyOp OP>
class traitsFieldsMultiply
{ };


/*! Direct multiplication */
template<typename DOFTYPEA, typename DOFTYPEB>
class traitsFieldsMultiply<DOFTYPEA,DOFTYPEB,fmMultiply>
{
  public:
    typedef typename traitsMultiply<DOFTYPEA,DOFTYPEB>::DIRECTTYPE OUTTYPE;
  
  public:
    static OUTTYPE compute(const DOFTYPEA & dofA, const DOFTYPEB & dofB);
};

template<typename DOFTYPEA, typename DOFTYPEB>
typename traitsFieldsMultiply<DOFTYPEA,DOFTYPEB,fmMultiply>::OUTTYPE
traitsFieldsMultiply<DOFTYPEA,DOFTYPEB,fmMultiply>::
compute(const DOFTYPEA & dofA, const DOFTYPEB & dofB)
{
  return(traitsMultiply<DOFTYPEA,DOFTYPEB>::multiply(dofA,dofB));
}



//_________________________________________________________________________________________________
// FIELDS MULTIPLY
//-------------------------------------------------------------------------------------------------
/*! Class for binary operation on similar fields (same length of dofVect) */
template<typename FIELDA, typename FIELDB, typename FIELDO, fieldsMultiplyOp OP>
class fieldsMultiply
{
    /*! @name Typedefs */ //@{
  public:    
    typedef Teuchos::RCP<const FIELDA>  CONSTFIELDA_RCP;
    typedef Teuchos::RCP<const FIELDB>  CONSTFIELDB_RCP;
    typedef Teuchos::RCP<FIELDO>        FIELDO_RCP;
    
    typedef typename FIELDA::FIELD_DOFTYPE  DOFTYPEA;
    typedef typename FIELDB::FIELD_DOFTYPE  DOFTYPEB;
    typedef typename FIELDO::FIELD_DOFTYPE  DOFTYPEO;
    
    typedef traitsFieldsMultiply<DOFTYPEA,DOFTYPEB,OP> TRAITSMULTIPLY;
    //@}
    
    /*! @name Multiply */ //@{
  public:
    static void multiply(const FIELDA & fieldA, const FIELDB & fieldB, FIELDO & fieldO);
    static void multiply(const CONSTFIELDA_RCP & fieldA, const CONSTFIELDB_RCP & fieldB, const FIELDO_RCP & fieldO);
    //@}
};


template<typename FIELDA, typename FIELDB, typename FIELDO, fieldsMultiplyOp OP>
void
fieldsMultiply<FIELDA,FIELDB,FIELDO,OP>::
multiply(const FIELDA & fieldA, const FIELDB & fieldB, FIELDO & fieldO)
{
  typedef typename TRAITSMULTIPLY::OUTTYPE  OUTDOFTYPE;
  assert(staticAssert<traitsBasic<OUTDOFTYPE>::myType == traitsBasic<DOFTYPEO>::myType>::returnValue);
  assert(fieldA.getDofVect().size() == fieldB.getDofVect().size());
  assert(fieldA.getDofVect().size() == fieldO.getDofVect().size());
  
  OUTDOFTYPE dofType;
  
  for(UInt i=1; i <= fieldA.getDofVect().size(); ++i)
  {
    dofType = TRAITSMULTIPLY::compute(fieldA.getDofL(i), fieldB.getDofL(i));
    fieldO.setDofL(i,dofType);
  }
}

template<typename FIELDA, typename FIELDB, typename FIELDO, fieldsMultiplyOp OP>
void
fieldsMultiply<FIELDA,FIELDB,FIELDO,OP>::
multiply(const CONSTFIELDA_RCP & fieldA, const CONSTFIELDB_RCP & fieldB, const FIELDO_RCP & fieldO)
{
  assert(fieldA.total_count() != 0);
  assert(fieldB.total_count() != 0);
  assert(fieldO.total_count() != 0);
  
  multiply(*fieldA,*fieldB,*fieldO);
}

#endif
