/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSVECTORBUILDER_HPP
#define TRAITSVECTORBUILDER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_Comm.hpp"

#include "Kokkos_DefaultNode.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_DataAccess.h"

#include "morganaFiniteElements.hpp"
#include "traitsBasic.h"
#include "morganaGeometry.hpp"

using namespace std;
using namespace boost::mpi;



//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------

/*! The vector building class. Is split in two specializations. */
template<typename LOCVECTOR, typename INTPMAPTYPE> class traitsVectorBuilder;




//_________________________________________________________________________________________________
// SPECIALIZATION PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! The matrix building class. The \c pMapItem specialization */
template<typename LOCVECTOR>
class traitsVectorBuilder<LOCVECTOR,pMapItem>
{
    /*! @name Epetra/Tpetra Typedefs */ //@{
  public:
    typedef Teuchos::MpiComm<int>    COMM;
    typedef Tpetra::Map<>            TPETRA_MAP;
    typedef Tpetra::MultiVector<>    TPETRA_VECTOR;
    
    typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
    typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef typename TPETRA_VECTOR::scalar_type      TPETRA_SCALAR;
      
    typedef Epetra_FEVector  EPETRA_VECTOR;
    //@}
    
    /*! @name Morgana Typedefs */ //@{
  public:
    typedef typename LOCVECTOR::FUNCTIONAL   FUNCTIONAL;
    typedef typename LOCVECTOR::INTCARD      INTCARD;
    typedef typename LOCVECTOR::INTCARDS     INTCARDS;
    typedef typename LOCVECTOR::PMAPTYPE     PMAPTYPE;
    
    typedef typename FUNCTIONAL::INTGRID       INTGRID;
    typedef typename FUNCTIONAL::TEST_DOFTYPE  TEST_DOFTYPE;
    typedef typename FUNCTIONAL::TEST_LISTMAP  ROW_MAP;
    //@}
    
    /*! @name Max values */ //@{
  public:
    static const UInt maxI_test  = traitsBasic<TEST_DOFTYPE>::numI;
    static const UInt maxJ_test  = traitsBasic<TEST_DOFTYPE>::numJ;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool operatorLoaded;
    bool commDevLoaded;
    Teuchos::RCP<FUNCTIONAL>   op;
    Teuchos::RCP<communicator> commDev;
    LOCVECTOR                  localVector;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    traitsVectorBuilder();
    traitsVectorBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    traitsVectorBuilder(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    void setFunctional(const Teuchos::RCP<FUNCTIONAL> & Op);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    void setIntCardL(const UInt & lid, const INTCARD & IntCard);
    void setIntCardG(const UInt & gid, const INTCARD & IntCard);
    void setIntCards(const INTCARDS & INTCards);
    //@}
    
    /*! @name Build */ //@{
  public:
    void buildEpetraVector(Teuchos::RCP<EPETRA_VECTOR> & vector);
    void buildTpetraVector(Teuchos::RCP<TPETRA_VECTOR> & vector);
    //@}
};


template<typename LOCVECTOR>
traitsVectorBuilder<LOCVECTOR,pMapItem>::
traitsVectorBuilder()
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapPlain>::returnValue);

  operatorLoaded = false;
  commDevLoaded  = false;
}

template<typename LOCVECTOR>
traitsVectorBuilder<LOCVECTOR,pMapItem>::
traitsVectorBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapPlain>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  op      = Op;
  commDev = CommDev;
  
  localVector.setFunctional(Op);
}

template<typename LOCVECTOR>
traitsVectorBuilder<LOCVECTOR,pMapItem>::
traitsVectorBuilder(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapPlain>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  op      = Op;
  commDev = Teuchos::rcpFromRef(CommDev);
  
  localVector.setFunctional(Op);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
setFunctional(const Teuchos::RCP<FUNCTIONAL> & Op)
{
  operatorLoaded = true;
  op = Op;
  localVector.setFunctional(Op);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
setIntCardL(const UInt & lid, const INTCARD & IntCard)
{
  localVector.setIntCardL(lid,IntCard);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
setIntCardG(const UInt & gid, const INTCARD & IntCard)
{
  localVector.setIntCardG(gid,IntCard);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
setIntCards(const INTCARDS & INTCards)
{
  localVector.setIntCards(INTCards);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
buildEpetraVector(Teuchos::RCP<EPETRA_VECTOR> & vector)
{
  //Alloc and assert
  assert(operatorLoaded);
  assert(commDevLoaded);
  
  Teuchos::RCP<INTGRID> intGrid = op->getIntegrationGrid();
  assert(intGrid->getMeshStandard() == STDB);  
  
  Epetra_MpiComm  epetraComm(*commDev);
  
  UInt sizeRow;
  int    * pointerRow;
  double * pointerVect;
  sVect<UInt> indexRow;
  sVect<Real> vect;
  
  
  //Maps and matrix creation
  ROW_MAP morganaRowMap = op->getListMap_test();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  
  Teuchos::RCP<Epetra_Map> epetraRowMap;
  mapManip.exportEpetraMap(morganaRowMap,epetraRowMap,epetraComm,0);
  
  vector = Teuchos::rcp(new Epetra_FEVector(*epetraRowMap));
  
  
  //Main loop
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      op->setIJ(Itest, Jtest);
      
      for(UInt i=1; i <= intGrid->getNumElements(); ++i)
      {
        op->setElement(i);

        sizeRow = localVector.numIndex_row();
        indexRow.resize(sizeRow);

        localVector.indexG_row(indexRow);
        vect = localVector.vector();

        pointerRow  = (int*)    & indexRow[0];
        pointerVect = (double*) & vect[0];

        vector->SumIntoGlobalValues(sizeRow, pointerRow, pointerVect);
      }
    }
  }
  
  //Vector assemble
  vector->GlobalAssemble();
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItem>::
buildTpetraVector(Teuchos::RCP<TPETRA_VECTOR> & vector)
{
  //Alloc and assert-------------------------------------------------------------------------------
  assert(operatorLoaded);
  assert(commDevLoaded);
  
  Teuchos::RCP<INTGRID> intGrid = op->getIntegrationGrid();
  assert(intGrid->getMeshStandard() == STDB);  
  
  Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::rcp(new Teuchos::MpiComm<int> (*commDev));
  
  //Maps and matrix creation-----------------------------------------------------------------------
  ROW_MAP morganaRowMap = op->getListMap_test();
  pVect<Real,pMapItem> sumVect(morganaRowMap.size());
  sumVect.setMap(morganaRowMap);
  sumVect.updateFinder();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  
  Teuchos::RCP<const TPETRA_MAP> tpetraRowMap;
  mapManip.exportTpetraMap(morganaRowMap,tpetraRowMap);
  
  vector = Teuchos::rcp(new TPETRA_VECTOR(tpetraRowMap,1) );
  
  //Alloc------------------------------------------------------------------------------------------
  UInt sizeRow;
  sVect<UInt> indexRow;
  sVect<Real> vect;
  
  //Main loop--------------------------------------------------------------------------------------
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      op->setIJ(Itest, Jtest);
      
      for(UInt i=1; i <= intGrid->getNumElements(); ++i)
      {
        if(intGrid->getElements().getRowMapL(i).getOwned())
        {
          op->setElement(i);

          sizeRow = localVector.numIndex_row();
          indexRow.resize(sizeRow);
          
          localVector.indexL_row(indexRow,1);
          vect = localVector.vector();          
          
          for(UInt k=1; k <= indexRow.size(); ++k)
          { sumVect(indexRow(k)) += vect(k); }
        }
      } //End cycle elements
    }
  }
  
  //Vector assemble--------------------------------------------------------------------------------  
  pVectGlobalManip<Real,pMapItem> vectorManip(commDev);
  std::plus<Real> op;
  vectorManip.allReduce(sumVect,op);
  
  //Translation------------------------------------------------------------------------------------
  UInt gid;
  Real value;
  
  for(UInt i=1; i <= morganaRowMap.size(); ++i)
  {
    gid   = morganaRowMap(i).getGid();
    value = sumVect.getDataG(gid);
    vector->sumIntoGlobalValue(TPETRA_GLOBAL_ORDINAL(gid-1),
                               0,
                               TPETRA_SCALAR(value));
  }
}



//_________________________________________________________________________________________________
// SPECIALIZATION PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------


/*! The matrix building class. The \c pMapItemShare specialization */
template<typename LOCVECTOR>
class traitsVectorBuilder<LOCVECTOR,pMapItemShare>
{
  /*! @name Epetra/Tpetra Typedefs */ //@{
  public:
    typedef Teuchos::MpiComm<int>    COMM;
    typedef Tpetra::Map<>            TPETRA_MAP;
    typedef Tpetra::MultiVector<>    TPETRA_VECTOR;
    
    typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
    typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef typename TPETRA_VECTOR::scalar_type      TPETRA_SCALAR;
      
    typedef Epetra_FEVector  EPETRA_VECTOR;
    //@}
    
    /*! @name Morgana Typedefs */ //@{
  public:
    typedef typename LOCVECTOR::FUNCTIONAL   FUNCTIONAL;
    typedef typename LOCVECTOR::INTCARD      INTCARD;
    typedef typename LOCVECTOR::INTCARDS     INTCARDS;
    typedef typename LOCVECTOR::PMAPTYPE     PMAPTYPE;
    
    typedef typename FUNCTIONAL::INTGRID       INTGRID;
    typedef typename FUNCTIONAL::TEST_DOFTYPE  TEST_DOFTYPE;
    typedef typename FUNCTIONAL::TEST_LISTMAP  ROW_MAP;
    //@}
    
    /*! @name Max values */ //@{
  public:
    static const UInt maxI_test  = traitsBasic<TEST_DOFTYPE>::numI;
    static const UInt maxJ_test  = traitsBasic<TEST_DOFTYPE>::numJ;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool operatorLoaded;
    bool commDevLoaded;
    Teuchos::RCP<FUNCTIONAL>   op;
    Teuchos::RCP<communicator> commDev;
    LOCVECTOR         localVector;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    traitsVectorBuilder();
    traitsVectorBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    traitsVectorBuilder(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op);
    void setFunctional(const Teuchos::RCP<FUNCTIONAL> & Op);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    void setIntCardL(const UInt & lid, const INTCARD & IntCard);
    void setIntCardG(const UInt & gid, const INTCARD & IntCard);
    void setIntCards(const INTCARDS & INTCards);
    //@}
    
    /*! @name Build */ //@{
  public:
    void buildEpetraVector(Teuchos::RCP<EPETRA_VECTOR> & vector);
    void buildTpetraVector(Teuchos::RCP<TPETRA_VECTOR> & vector);
    //@}
};


template<typename LOCVECTOR>
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
traitsVectorBuilder()
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapShare>::returnValue);

  operatorLoaded = false;
  commDevLoaded  = false;
}

template<typename LOCVECTOR>
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
traitsVectorBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapShare>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  op      = Op;
  commDev = CommDev;
  
  localVector.setFunctional(Op);
}

template<typename LOCVECTOR>
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
traitsVectorBuilder(communicator & CommDev, const Teuchos::RCP<FUNCTIONAL> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapShare>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  op      = Op;
  commDev = Teuchos::rcpFromRef(CommDev);
  
  localVector.setFunctional(Op);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
setFunctional(const Teuchos::RCP<FUNCTIONAL> & Op)
{
  operatorLoaded = true;
  op = Op;
  localVector.setFunctional(Op);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
setIntCardL(const UInt & lid, const INTCARD & IntCard)
{
  localVector.setIntCardL(lid,IntCard);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
setIntCardG(const UInt & gid, const INTCARD & IntCard)
{
  localVector.setIntCardG(gid,IntCard);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
setIntCards(const INTCARDS & INTCards)
{
  localVector.setIntCards(INTCards);
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
buildEpetraVector(Teuchos::RCP<EPETRA_VECTOR> & vector)
{
  //Alloc and assert
  assert(operatorLoaded);
  assert(commDevLoaded);
  
  Teuchos::RCP<INTGRID> intGrid = op->getIntegrationGrid();  
  Epetra_MpiComm  epetraComm(*commDev);
  
  UInt sizeRow;
  int    * pointerRow;
  double * pointerVect;
  sVect<UInt> indexRow;
  sVect<Real> vect;
  
  
  //Maps and matrix creation
  ROW_MAP morganaRowMap = op->getListMap_test();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  
  Teuchos::RCP<Epetra_Map> epetraRowMap;
  mapManip.exportEpetraMap(morganaRowMap,epetraRowMap,epetraComm,0);
  
  vector = Teuchos::rcp(new Epetra_FEVector(*epetraRowMap));
  

  //Main loop
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      op->setIJ(Itest, Jtest);
      
      for(UInt i=1; i <= intGrid->getNumElements(); ++i)
      {
        if(intGrid->getElements().getRowMapL(i).getOwned())
        {
          op->setElement(i);

          sizeRow = localVector.numIndex_row();
          indexRow.resize(sizeRow);

          localVector.indexG_row(indexRow);
          vect = localVector.vector();

          pointerRow  = (int*)    & indexRow[0];
          pointerVect = (double*) & vect[0];

          vector->SumIntoGlobalValues(sizeRow, pointerRow, pointerVect);
        }
      } //End cycle elements
    }
  }
  
  //Vector assemble
  vector->GlobalAssemble();
}

template<typename LOCVECTOR>
void
traitsVectorBuilder<LOCVECTOR,pMapItemShare>::
buildTpetraVector(Teuchos::RCP<TPETRA_VECTOR> & vector)
{
  //Alloc and assert-------------------------------------------------------------------------------
  assert(operatorLoaded);
  assert(commDevLoaded);

  Teuchos::RCP<INTGRID> intGrid = op->getIntegrationGrid();
  Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::rcp(new Teuchos::MpiComm<int> (*commDev));
  
  //Maps and matrix creation-----------------------------------------------------------------------
  ROW_MAP morganaRowMap = op->getListMap_test();
  pVect<Real,pMapItemShare> sumVect(morganaRowMap.size());
  sumVect.setMap(morganaRowMap);
  sumVect.updateFinder();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  
  Teuchos::RCP<const TPETRA_MAP> tpetraRowMap;
  mapManip.exportTpetraMap(morganaRowMap,tpetraRowMap);
  
  vector = Teuchos::rcp(new TPETRA_VECTOR(tpetraRowMap,1) );
  
  //Alloc------------------------------------------------------------------------------------------
  UInt sizeRow;
  sVect<UInt> indexRow;
  sVect<Real> vect;

  //Main loop--------------------------------------------------------------------------------------
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      op->setIJ(Itest, Jtest);
      
      for(UInt i=1; i <= intGrid->getNumElements(); ++i)
      {
        if(intGrid->getElements().getRowMapL(i).getOwned())
        {
          op->setElement(i);

          sizeRow = localVector.numIndex_row();
          indexRow.resize(sizeRow);
          
          localVector.indexL_row(indexRow,1);
          vect = localVector.vector();          
          
          for(UInt k=1; k <= indexRow.size(); ++k)
          { sumVect(indexRow(k)) += vect(k); }
        }
      } //End cycle elements
    }
  }
  
  //Vector assemble--------------------------------------------------------------------------------  
  pVectGlobalManip<Real,pMapItemShare> vectorManip(commDev);
  std::plus<Real> op;
  vectorManip.allReduce(sumVect,op);
  
  //Translation------------------------------------------------------------------------------------
  UInt gid;
  Real value;
  
  for(UInt i=1; i <= morganaRowMap.size(); ++i)
  {
    gid   = morganaRowMap(i).getGid();
    value = sumVect.getDataG(gid);
    vector->sumIntoGlobalValue(TPETRA_GLOBAL_ORDINAL(gid-1),
                               0,
                               TPETRA_SCALAR(value));
  }
}


#endif
