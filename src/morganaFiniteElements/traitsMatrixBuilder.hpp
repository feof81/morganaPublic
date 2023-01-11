/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef TRAITSMATRIXBUILDER_HPP
#define TRAITSMATRIXBUILDER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ArrayViewDecl.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Kokkos_DefaultNode.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_DataAccess.h"

#include "morganaFiniteElements.hpp"
#include "traitsBasic.h"
#include "morganaGeometry.hpp"

using namespace std;
using namespace boost::mpi;



//_________________________________________________________________________________________________
// EPETRA SWITCH
//-------------------------------------------------------------------------------------------------

/*! Traits for epetra matrix */
template<OPClass CLASS>
struct traitsMajor;

/*! Traits for epetra matrix: the rowMajor specialization */
template<>
struct traitsMajor<opRowMajor>
{
  static const int epetraMajor = Epetra_FECrsMatrix::ROW_MAJOR;
};

/*! Traits for epetra matrix: the colMajor specialization */
template<>
struct traitsMajor<opColMajor>
{
  static const int epetraMajor = Epetra_FECrsMatrix::COLUMN_MAJOR;
};



//_________________________________________________________________________________________________
// TPETRA SWITCH
//-------------------------------------------------------------------------------------------------

/*! Traits for tpetra matrix */
template<OPClass CLASS>
struct traitsMajorTpetra;

/*! Traits for tpetra matrix: the rowMajor specialization */
template<>
struct traitsMajorTpetra<opRowMajor>
{
  typedef Teuchos::MpiComm<int>    COMM;
  typedef Tpetra::Map<>            TPETRA_MAP;
  typedef Tpetra::CrsMatrix<>      TPETRA_CRS;
  typedef Teuchos::RCP<TPETRA_CRS> TPETRA_CRS_RCP;
   
  typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
  typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
  typedef typename TPETRA_CRS::scalar_type         TPETRA_SCALAR;
  
  void tpetraCrsInsert(const TPETRA_CRS_RCP & globMat,
                       const sVect<UInt> & indexRow,
                       const sVect<UInt> & indexCol,
                       const sVect<Real> & locMat)
  {
    UInt k = 1;
    sVect<TPETRA_GLOBAL_ORDINAL> dummyCols(indexCol.size());
    sVect<TPETRA_SCALAR>         dummyVals(indexCol.size());
    
    Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosCols(dummyCols);
    Teuchos::ArrayView<TPETRA_SCALAR>         teuchosVals(dummyVals);
    
    //Translate columns
    for(UInt j=1; j <= indexCol.size(); ++j)
    { teuchosCols[j-1] = indexCol(j); };
    
    //Insert values
    for(UInt i=1; i <= indexRow.size(); ++i)
    {
      //Translate the 
      for(UInt j=1; j <= indexCol.size(); ++j)
      { 
        teuchosVals[j-1] = locMat(k);
        ++k;
      }
      
      //Insertion in the matrix
      globMat->insertGlobalValues(indexRow(i), teuchosCols, teuchosVals);
    }
  }
};

/*! Traits for tpetra matrix: the colMajor specialization */
template<>
struct traitsMajorTpetra<opColMajor>
{
  typedef Teuchos::MpiComm<int>    COMM;
  typedef Tpetra::Map<>            TPETRA_MAP;
  typedef Tpetra::CrsMatrix<>      TPETRA_CRS;
  typedef Teuchos::RCP<TPETRA_CRS> TPETRA_CRS_RCP;
    
  typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
  typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
  typedef typename TPETRA_CRS::scalar_type         TPETRA_SCALAR;
  
  void tpetraCrsInsert(const TPETRA_CRS_RCP               & globMat,
                       const sVect<UInt> & indexRow,
                       const sVect<UInt> & indexCol,
                       const sVect<Real> & locMat)
  {
    UInt k = 1;
    sVect<TPETRA_GLOBAL_ORDINAL> dummyCols(indexCol.size());
    sVect<TPETRA_SCALAR>         dummyVals(indexCol.size());
    
    Teuchos::ArrayView<TPETRA_GLOBAL_ORDINAL> teuchosCols(dummyCols);
    Teuchos::ArrayView<TPETRA_SCALAR>         teuchosVals(dummyVals);
    
    //Translate columns
    for(UInt j=1; j <= indexCol.size(); ++j)
    { teuchosCols[j-1] = indexCol(j); };
    
    //Insert values
    for(UInt i=1; i <= indexRow.size(); ++i)
    {
      //Translate the 
      for(UInt j=1; j <= indexCol.size(); ++j)
      { 
        k = (j-1) * indexRow.size() + i;
        teuchosVals[j-1] = locMat(k);
      }
      
      //Insertion in the matrix
      globMat->insertGlobalValues(indexRow(i), teuchosCols, teuchosVals);
    }
  }
};




//_________________________________________________________________________________________________
// UNSPECIALIZED
//-------------------------------------------------------------------------------------------------

/*! The matrix building class. Is split in two specializations. */
template<typename LOCMATRIX, typename INTPMAPTYPE> class traitsMatrixBuilder;



//_________________________________________________________________________________________________
// SPECIALIZATION PMAPITEM
//-------------------------------------------------------------------------------------------------

/*! The matrix building class. The \c pMapItem specialization */
template<typename LOCMATRIX>
class traitsMatrixBuilder<LOCMATRIX,pMapItem>
{
    /*! @name Epetra/Tpetra Typedefs */ //@{
  public:
    typedef Teuchos::MpiComm<int> COMM;
    typedef Tpetra::Map<>         TPETRA_MAP;
    typedef Tpetra::CrsMatrix<>   TPETRA_CRS;
    
    typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
    typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef typename TPETRA_CRS::scalar_type         TPETRA_SCALAR;

    typedef Epetra_FECrsMatrix EPETRA_CRS;
    //@}
    
    /*! @name Morgana Typedefs */ //@{
  public:
    typedef typename LOCMATRIX::OPERATOR     OPERATOR;
    typedef typename LOCMATRIX::INTCARD      INTCARD;
    typedef typename LOCMATRIX::INTCARDS     INTCARDS;
    typedef typename LOCMATRIX::PMAPTYPE     PMAPTYPE;
    
    typedef typename OPERATOR::LOOPGRID      LOOPGRID;
    typedef typename OPERATOR::FIELD_DOFTYPE FIELD_DOFTYPE;
    typedef typename OPERATOR::TEST_DOFTYPE  TEST_DOFTYPE;
    typedef typename OPERATOR::TEST_LISTMAP  ROW_MAP;
    typedef typename OPERATOR::FIELD_LISTMAP COL_MAP;
    
    typedef pVect<UInt,pMapItem>            PVECT_NUMROWS;
    typedef pVectComm<UInt,pMapItem>        PCOMM_NUMROWS;
    typedef pVectGlobalManip<UInt,pMapItem> PGMANIP_NUMROWS;
    //@}
    
    /*! @name Max values */ //@{
  public:
    static const UInt maxI_field = traitsBasic<FIELD_DOFTYPE>::numI;
    static const UInt maxJ_field = traitsBasic<FIELD_DOFTYPE>::numJ;
    static const UInt maxI_test  = traitsBasic<TEST_DOFTYPE>::numI;
    static const UInt maxJ_test  = traitsBasic<TEST_DOFTYPE>::numJ;
    //@}
  
    /*! @name Internal data */ //@{
  public:
    bool operatorLoaded;
    bool commDevLoaded;
    Teuchos::RCP<OPERATOR> op;
    Teuchos::RCP<communicator> commDev;
    LOCMATRIX  localMatrix;
    UInt allocSize;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    traitsMatrixBuilder();
    traitsMatrixBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<OPERATOR> & Op);
    traitsMatrixBuilder(communicator & CommDev, const Teuchos::RCP<OPERATOR> & Op);
    void setOperator(const Teuchos::RCP<OPERATOR> & Op);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    void setIntCardL(const UInt & lid, const INTCARD & IntCard);
    void setIntCardG(const UInt & gid, const INTCARD & IntCard);
    void setIntCards(const INTCARDS & INTCards);
    void setAllocSize(const UInt & AllocSize);
    //@}
    
    /*! @name Build */ //@{
  public:
    sVect<UInt> numEntPerRow();
    void buildEpetraCrs(Teuchos::RCP<EPETRA_CRS> & matrix);
    void buildTpetraCrs(Teuchos::RCP<TPETRA_CRS> & matrix);
    //@}
};


//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename LOCMATRIX>
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
traitsMatrixBuilder()
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapPlain>::returnValue);

  operatorLoaded = false;
  commDevLoaded  = false;
  
  allocSize = 1;
}

template<typename LOCMATRIX>
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
traitsMatrixBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<OPERATOR> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapPlain>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  allocSize = 1;
  op      = Op;
  commDev = CommDev;
  
  localMatrix.setOperator(Op);
}

template<typename LOCMATRIX>
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
traitsMatrixBuilder(communicator & CommDev, const Teuchos::RCP<OPERATOR> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapPlain>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  allocSize = 1;
  op      = Op;
  commDev = Teuchos::rcpFromRef(CommDev);
  
  localMatrix.setOperator(Op);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
setOperator(const Teuchos::RCP<OPERATOR> & Op)
{
  operatorLoaded = true;
  op = Op;
  localMatrix.setOperator(Op);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
setIntCardL(const UInt & lid, const INTCARD & IntCard)
{
  localMatrix.setIntCardL(lid,IntCard);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
setIntCardG(const UInt & gid, const INTCARD & IntCard)
{
  localMatrix.setIntCardG(gid,IntCard);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
setIntCards(const INTCARDS & INTCards)
{
  localMatrix.setIntCards(INTCards);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
setAllocSize(const UInt & AllocSize)
{
  allocSize = AllocSize;
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
buildEpetraCrs(Teuchos::RCP<EPETRA_CRS> & matrix)
{
  //Alloc and assert
  assert(operatorLoaded);
  assert(commDevLoaded);
  
  Teuchos::RCP<LOOPGRID> loopGrid = op->getLoopGrid();
  assert(loopGrid->getMeshStandard() == STDB);  
  
  Epetra_MpiComm  epetraComm(*commDev);
  
  UInt major = traitsMajor<OPERATOR::opClass>::epetraMajor;
  
  UInt sizeRow, sizeCol;
  int    * pointerRow, * pointerCol;
  double * pointerMat;
  sVect<UInt> indexRow, indexCol;
  sVect<Real> mat;
  
  
  //Maps and matrix creation
  ROW_MAP morganaRowMap = op->getListMap_test();
  COL_MAP morganaColMap = op->getListMap_field();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  mapManip.destroyOverlap(morganaColMap);
  
  Teuchos::RCP<Epetra_Map> epetraRowMap;
  Teuchos::RCP<Epetra_Map> epetraColMap;
  
  mapManip.exportEpetraMap(morganaRowMap,epetraRowMap,epetraComm,0);
  mapManip.exportEpetraMap(morganaColMap,epetraColMap,epetraComm,0);
  
  matrix = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*epetraRowMap, allocSize));
  
  
  //Main loop
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      for(UInt Ifield = 1; Ifield <= maxI_field; ++Ifield)
      {
        for(UInt Jfield = 1; Jfield <= maxJ_field; ++Jfield)
        {
          op->setIJ(Itest, Jtest, Ifield, Jfield);
    
          for(UInt i=1; i <= loopGrid->getNumElements(); ++i)
          {    
            op->setElement(i);
    
            sizeRow = localMatrix.numIndex_row();
            sizeCol = localMatrix.numIndex_col();
    
            indexRow.resize(sizeRow);
            indexCol.resize(sizeCol);
    
            localMatrix.indexG_row(indexRow);
            localMatrix.indexG_col(indexCol);
            mat = localMatrix.matrix();
    
            pointerRow = (int*)    & indexRow[0];
            pointerCol = (int*)    & indexCol[0];
            pointerMat = (double*) & mat[0];
    
            matrix->InsertGlobalValues(sizeRow, pointerRow, sizeCol, pointerCol, pointerMat, major);
          }
        }
      }
    }
  }
  
  //Matrix assemble
  matrix->GlobalAssemble(*epetraColMap, *epetraRowMap);
}

template<typename LOCMATRIX>
sVect<UInt>
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
numEntPerRow()
{
  //Alloc------------------------------------------------------------
  Teuchos::RCP<LOOPGRID> loopGrid = op->getLoopGrid();
  UInt sizeRow, sizeCol;
  sVect<UInt> indexRow, indexCol;
    
  //Maps and vectors-------------------------------------------------
  ROW_MAP morganaRowMap = op->getListMap_test();
  sVect<UInt>       voidData(morganaRowMap.size());
  sVect<set<UInt> > rowIndexes(morganaRowMap.size());
  PVECT_NUMROWS     numRowsVect(morganaRowMap,voidData);
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  
  //Compute local num entries----------------------------------------
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      for(UInt Ifield = 1; Ifield <= maxI_field; ++Ifield)
      {
        for(UInt Jfield = 1; Jfield <= maxJ_field; ++Jfield)
        {
          op->setIJ(Itest, Jtest, Ifield, Jfield);
    
          for(UInt i=1; i <= loopGrid->getNumElements(); ++i)
          {
            op->setElement(i);
    
            sizeRow = localMatrix.numIndex_row();
            sizeCol = localMatrix.numIndex_col();
    
            indexRow.resize(sizeRow);
            indexCol.resize(sizeCol);
    
            localMatrix.indexL_row(indexRow);
            localMatrix.indexL_col(indexCol);
            
            for(UInt i=1; i <= sizeRow; ++i)
            {
              for(UInt j=1; j <= sizeCol; ++j)
              { rowIndexes(indexRow(i)+1).insert(indexCol(j)+1); }
            }
          }
        }
      }
    }
  }
  
  //Parallel sum-----------------------------------------------------------------------------------
  UInt gid;
  UInt minGid = std::numeric_limits<UInt>::max();
  UInt maxGid = 0;
  
  //Fill vect
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  { numRowsVect(i) = rowIndexes(i).size(); }
  
  //Normal distribution
  PCOMM_NUMROWS numRowsComm(commDev);
  numRowsComm.vectorNormal(numRowsVect);
  
  //Compute the sum
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  {
    minGid = std::min(minGid, numRowsVect.getMapL(i).getGid());
    maxGid = std::max(maxGid, numRowsVect.getMapL(i).getGid());
  }
  
  sVect<UInt> sumVect(maxGid -minGid +1);
  
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  {
    gid = numRowsVect.getMapL(i).getGid();
    sumVect(gid -minGid +1) += numRowsVect(i);
  }
  
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  {
    gid = numRowsVect.getMapL(i).getGid();
    numRowsVect(i) = sumVect(gid -minGid +1);
  }
  
  //Comm back
  numRowsComm.vectorPid(numRowsVect);
  
  //Destroy overlap--------------------------------------------------------------------------------
  sVect<UInt> outSize(morganaRowMap.size());
  
  for(UInt i=1; i <= morganaRowMap.size(); ++i)
  {
    gid        = morganaRowMap(i).getGid();
    outSize(i) = numRowsVect.getDataG(gid);
  }

  return(outSize);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItem>::
buildTpetraCrs(Teuchos::RCP<TPETRA_CRS> & matrix)
{
  //Alloc and assert-------------------------------------------------------------------------------
  assert(operatorLoaded);
  assert(commDevLoaded);
  
  Teuchos::RCP<LOOPGRID> loopGrid = op->getLoopGrid();
  assert(loopGrid->getMeshStandard() == STDB);
  
  static const OPClass major = OPERATOR::opClass;
  traitsMajorTpetra<major> tpetraManip;
  
  UInt sizeRow, sizeCol;
  sVect<UInt> indexRow, indexCol;
  sVect<Real> mat;
  
  //Maps and matrix creation-----------------------------------------------------------------------
  ROW_MAP morganaRowMap = op->getListMap_test();
  COL_MAP morganaColMap = op->getListMap_field();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  mapManip.destroyOverlap(morganaColMap);
  
  Teuchos::RCP<const TPETRA_MAP> tpetraRowMap;
  Teuchos::RCP<const TPETRA_MAP> tpetraColMap;
  
  mapManip.exportTpetraMap(morganaRowMap, tpetraRowMap);
  mapManip.exportTpetraMap(morganaColMap, tpetraColMap);
  
  //Estimate matrix size and alloc-----------------------------------------------------------------
  sVect<UInt> rowIndexes = numEntPerRow();
  
  sVect<TPETRA_GLOBAL_TYPE> totalSize(rowIndexes.size());
  for(UInt i=1; i <= rowIndexes.size(); ++i)
  { totalSize(i) = TPETRA_GLOBAL_TYPE(rowIndexes(i)); }
  
  Teuchos::ArrayView<const TPETRA_GLOBAL_TYPE> numEntPerRow(totalSize);
  matrix = rcp( new TPETRA_CRS(tpetraRowMap, numEntPerRow) );
  
  //Main loop--------------------------------------------------------------------------------------
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      for(UInt Ifield = 1; Ifield <= maxI_field; ++Ifield)
      {
        for(UInt Jfield = 1; Jfield <= maxJ_field; ++Jfield)
        {
          op->setIJ(Itest, Jtest, Ifield, Jfield);
    
          for(UInt i=1; i <= loopGrid->getNumElements(); ++i)
          {
            op->setElement(i);
    
            sizeRow = localMatrix.numIndex_row();
            sizeCol = localMatrix.numIndex_col();
    
            indexRow.resize(sizeRow);
            indexCol.resize(sizeCol);
    
            localMatrix.indexG_row(indexRow);
            localMatrix.indexG_col(indexCol);
            mat = localMatrix.matrix();
    
            tpetraManip.tpetraCrsInsert(matrix, indexRow, indexCol, mat);
          }
        }
      }
    }
  }
  
  //Matrix assemble--------------------------------------------------------------------------------
  matrix->fillComplete(tpetraColMap, tpetraRowMap);
}
    


//_________________________________________________________________________________________________
// SPECIALIZATION PMAPITEMSHARE
//-------------------------------------------------------------------------------------------------

/*! The matrix building class. The \c pMapItemShare specialization */
template<typename LOCMATRIX>
class traitsMatrixBuilder<LOCMATRIX,pMapItemShare>
{
  /*! @name Epetra/Tpetra Typedefs */ //@{
  public:
    typedef Teuchos::MpiComm<int> COMM;
    typedef Tpetra::Map<>         TPETRA_MAP;
    typedef Tpetra::CrsMatrix<>   TPETRA_CRS;
    
    typedef Tpetra::global_size_t                    TPETRA_GLOBAL_TYPE;
    typedef typename TPETRA_MAP::global_ordinal_type TPETRA_GLOBAL_ORDINAL;
    typedef typename TPETRA_CRS::scalar_type         TPETRA_SCALAR;
    
    typedef Epetra_FECrsMatrix EPETRA_CRS;
    //@}
    
    /*! @name Morgana Typedefs */ //@{
  public:
    typedef typename LOCMATRIX::OPERATOR     OPERATOR;
    typedef typename LOCMATRIX::INTCARD      INTCARD;
    typedef typename LOCMATRIX::INTCARDS     INTCARDS;
    typedef typename LOCMATRIX::PMAPTYPE     PMAPTYPE;
    
    typedef typename OPERATOR::LOOPGRID      LOOPGRID;
    typedef typename OPERATOR::FIELD_DOFTYPE FIELD_DOFTYPE;
    typedef typename OPERATOR::TEST_DOFTYPE  TEST_DOFTYPE;
    typedef typename OPERATOR::TEST_LISTMAP  ROW_MAP;
    typedef typename OPERATOR::FIELD_LISTMAP COL_MAP;
    
    typedef pVect<UInt,pMapItemShare>            PVECT_NUMROWS;
    typedef pVectComm<UInt,pMapItemShare>        PCOMM_NUMROWS;
    typedef pVectGlobalManip<UInt,pMapItemShare> PGMANIP_NUMROWS;
    //@}
    
    /*! @name Max values */ //@{
  public:
    static const UInt maxI_field = traitsBasic<FIELD_DOFTYPE>::numI;
    static const UInt maxJ_field = traitsBasic<FIELD_DOFTYPE>::numJ;
    static const UInt maxI_test  = traitsBasic<TEST_DOFTYPE>::numI;
    static const UInt maxJ_test  = traitsBasic<TEST_DOFTYPE>::numJ;
    //@}
  
    /*! @name Internal data */ //@{
  public:
    bool operatorLoaded;
    bool commDevLoaded;
    Teuchos::RCP<OPERATOR> op;
    Teuchos::RCP<communicator> commDev;
    LOCMATRIX  localMatrix;
    UInt allocSize;
    //@}
    
    /*! @name Constructors and set */ //@{
  public:
    traitsMatrixBuilder();
    traitsMatrixBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<OPERATOR> & Op);
    traitsMatrixBuilder(communicator & CommDev, const Teuchos::RCP<OPERATOR> & Op);
    void setOperator(const Teuchos::RCP<OPERATOR> & Op);
    void setCommDev(const Teuchos::RCP<communicator> & CommDev);
    void setCommDev(communicator & CommDev);
    void setIntCardL(const UInt & lid, const INTCARD & IntCard);
    void setIntCardG(const UInt & gid, const INTCARD & IntCard);
    void setIntCards(const INTCARDS & INTCards);
    void setAllocSize(const UInt & AllocSize);
    //@}
    
    /*! @name Build */ //@{
  public:
    sVect<UInt> numEntPerRow();
    void buildEpetraCrs(Teuchos::RCP<EPETRA_CRS> & matrix);
    void buildTpetraCrs(Teuchos::RCP<TPETRA_CRS> & matrix);
    //@}
};


template<typename LOCMATRIX>
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
traitsMatrixBuilder()
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapShare>::returnValue);

  operatorLoaded = false;
  commDevLoaded  = false;
  
  allocSize = 1;
}

template<typename LOCMATRIX>
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
traitsMatrixBuilder(const Teuchos::RCP<communicator> & CommDev, const Teuchos::RCP<OPERATOR> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapShare>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  allocSize = 1;
  op      = Op;
  commDev = CommDev;
  
  localMatrix.setOperator(Op);
}

template<typename LOCMATRIX>
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
traitsMatrixBuilder(communicator & CommDev, const Teuchos::RCP<OPERATOR> & Op)
{
  assert(staticAssert<PMAPTYPE::parallelType == pMapShare>::returnValue);
  
  operatorLoaded = true;
  commDevLoaded  = true;
  
  allocSize = 1;
  op      = Op;
  commDev = Teuchos::rcpFromRef(CommDev);
  
  localMatrix.setOperator(Op);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
setOperator(const Teuchos::RCP<OPERATOR> & Op)
{
  operatorLoaded = true;
  op = Op;
  localMatrix.setOperator(Op);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
setCommDev(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
setCommDev(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
setIntCardL(const UInt & lid, const INTCARD & IntCard)
{
  localMatrix.setIntCardL(lid,IntCard);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
setIntCardG(const UInt & gid, const INTCARD & IntCard)
{
  localMatrix.setIntCardG(gid,IntCard);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
setIntCards(const INTCARDS & INTCards)
{
  localMatrix.setIntCards(INTCards);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
setAllocSize(const UInt & AllocSize)
{
  allocSize = AllocSize;
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
buildEpetraCrs(Teuchos::RCP<EPETRA_CRS> & matrix)
{
  //Alloc and assert
  assert(operatorLoaded);
  assert(commDevLoaded);
  
  Teuchos::RCP<LOOPGRID> loopGrid = op->getLoopGrid(); 
  
  Epetra_MpiComm  epetraComm(*commDev);
  
  int major = traitsMajor<OPERATOR::opClass>::epetraMajor;
  
  UInt sizeRow, sizeCol;
  int    * pointerRow, * pointerCol;
  double * pointerMat;
  sVect<UInt> indexRow, indexCol;
  sVect<Real> mat;
  
  
  //Maps and matrix creation
  ROW_MAP morganaRowMap = op->getListMap_test();
  COL_MAP morganaColMap = op->getListMap_field();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  mapManip.destroyOverlap(morganaColMap);
  
  Teuchos::RCP<Epetra_Map> epetraRowMap;
  Teuchos::RCP<Epetra_Map> epetraColMap;
  
  mapManip.exportEpetraMap(morganaRowMap,epetraRowMap,epetraComm,0);
  mapManip.exportEpetraMap(morganaColMap,epetraColMap,epetraComm,0);
  
  matrix = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*epetraRowMap, allocSize));
  
  //Main loop
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      for(UInt Ifield = 1; Ifield <= maxI_field; ++Ifield)
      {
        for(UInt Jfield = 1; Jfield <= maxJ_field; ++Jfield)
        {
	  op->setIJ(Itest, Jtest, Ifield, Jfield);
	  
	  for(UInt i=1; i <= loopGrid->getNumElements(); ++i)
          {
	    if(loopGrid->getElements().getRowMapL(i).getOwned())
	    {  
	      op->setElement(i);
	    
	      sizeRow = localMatrix.numIndex_row();
	      sizeCol = localMatrix.numIndex_col();
	    
	      indexRow.resize(sizeRow);
	      indexCol.resize(sizeCol);
	    
	      localMatrix.indexG_row(indexRow);
	      localMatrix.indexG_col(indexCol);
	      mat = localMatrix.matrix();
	     
	      pointerRow = (int*)    & indexRow[0];
	      pointerCol = (int*)    & indexCol[0];
	      pointerMat = (double*) & mat[0];
	    
	      matrix->InsertGlobalValues(sizeRow, pointerRow, sizeCol, pointerCol, pointerMat, major);
	    }
	  }  //End cycle elements
	  
	}  //End cycle types
      }
    }
  }
  
  //Matrix assemble
  matrix->GlobalAssemble(*epetraColMap, *epetraRowMap);
}

template<typename LOCMATRIX>
sVect<UInt>
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
numEntPerRow()
{
  //Alloc------------------------------------------------------------
  Teuchos::RCP<LOOPGRID> loopGrid = op->getLoopGrid();
  UInt sizeRow, sizeCol;
  sVect<UInt> indexRow, indexCol;
    
  //Maps and vectors-------------------------------------------------
  ROW_MAP morganaRowMap = op->getListMap_test();
  sVect<UInt>       voidData(morganaRowMap.size());
  sVect<set<UInt> > rowIndexes(morganaRowMap.size());
  PVECT_NUMROWS     numRowsVect(morganaRowMap,voidData);
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  
  //Compute local num entries----------------------------------------
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      for(UInt Ifield = 1; Ifield <= maxI_field; ++Ifield)
      {
        for(UInt Jfield = 1; Jfield <= maxJ_field; ++Jfield)
        {
          op->setIJ(Itest, Jtest, Ifield, Jfield);
    
          for(UInt i=1; i <= loopGrid->getNumElements(); ++i)
          {
            op->setElement(i);
    
            sizeRow = localMatrix.numIndex_row();
            sizeCol = localMatrix.numIndex_col();
    
            indexRow.resize(sizeRow);
            indexCol.resize(sizeCol);
    
            localMatrix.indexL_row(indexRow);
            localMatrix.indexL_col(indexCol);
            
            for(UInt i=1; i <= sizeRow; ++i)
            {
              for(UInt j=1; j <= sizeCol; ++j)
              { rowIndexes(indexRow(i)+1).insert(indexCol(j)+1); }
            }
          }
        }
      }
    }
  }
  
  //Parallel sum-----------------------------------------------------------------------------------
  UInt gid;
  UInt minGid = std::numeric_limits<UInt>::max();
  UInt maxGid = 0;
  
  //Fill vect
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  { numRowsVect(i) = rowIndexes(i).size(); }
  
  //Normal distribution
  PCOMM_NUMROWS numRowsComm(commDev);
  numRowsComm.vectorNormal(numRowsVect);
  
  //Compute the sum
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  {
    minGid = std::min(minGid, numRowsVect.getMapL(i).getGid());
    maxGid = std::max(maxGid, numRowsVect.getMapL(i).getGid());
  }
  
  sVect<UInt> sumVect(maxGid -minGid +1);
  
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  {
    gid = numRowsVect.getMapL(i).getGid();
    sumVect(gid -minGid +1) += numRowsVect(i);
  }
  
  for(UInt i=1; i <= numRowsVect.size(); ++i)
  {
    gid = numRowsVect.getMapL(i).getGid();
    numRowsVect(i) = sumVect(gid -minGid +1);
  }
  
  //Comm back
  numRowsComm.vectorPid(numRowsVect);
  
  //Destroy overlap--------------------------------------------------------------------------------
  sVect<UInt> outSize(morganaRowMap.size());
  
  for(UInt i=1; i <= morganaRowMap.size(); ++i)
  {
    gid        = morganaRowMap(i).getGid();
    outSize(i) = numRowsVect.getDataG(gid);
  }

  return(outSize);
}

template<typename LOCMATRIX>
void
traitsMatrixBuilder<LOCMATRIX,pMapItemShare>::
buildTpetraCrs(Teuchos::RCP<TPETRA_CRS> & matrix)
{
  //Alloc and assert-------------------------------------------------------------------------------
  assert(operatorLoaded);
  assert(commDevLoaded);
  
  Teuchos::RCP<LOOPGRID> loopGrid = op->getLoopGrid();
  
  static const OPClass major = OPERATOR::opClass;
  traitsMajorTpetra<major> tpetraManip;
  
  UInt sizeRow, sizeCol;
  sVect<UInt> indexRow, indexCol;
  sVect<Real> mat;
  
  //Maps and matrix creation-----------------------------------------------------------------------
  ROW_MAP morganaRowMap = op->getListMap_test();
  COL_MAP morganaColMap = op->getListMap_field();
  
  pMapGlobalManip<PMAPTYPE> mapManip(commDev);
  mapManip.destroyOverlap(morganaRowMap);
  mapManip.destroyOverlap(morganaColMap);
  
  Teuchos::RCP<const TPETRA_MAP> tpetraRowMap;
  Teuchos::RCP<const TPETRA_MAP> tpetraColMap;
  
  mapManip.exportTpetraMap(morganaRowMap, tpetraRowMap);
  mapManip.exportTpetraMap(morganaColMap, tpetraColMap);
  
  //Estimate matrix size and alloc-----------------------------------------------------------------
  sVect<UInt> rowIndexes = numEntPerRow();
  
  sVect<TPETRA_GLOBAL_TYPE> totalSize(rowIndexes.size());
  for(UInt i=1; i <= rowIndexes.size(); ++i)
  { totalSize(i) = TPETRA_GLOBAL_TYPE(rowIndexes(i)); }
  
  Teuchos::ArrayView<const TPETRA_GLOBAL_TYPE> numEntPerRow(totalSize);
  matrix = rcp( new TPETRA_CRS(tpetraRowMap, numEntPerRow) );
  
  //Main loop--------------------------------------------------------------------------------------
  for(UInt Itest = 1; Itest <= maxI_test; ++Itest)
  {
    for(UInt Jtest = 1; Jtest <= maxJ_test; ++Jtest)
    {
      for(UInt Ifield = 1; Ifield <= maxI_field; ++Ifield)
      {
        for(UInt Jfield = 1; Jfield <= maxJ_field; ++Jfield)
        {
          op->setIJ(Itest, Jtest, Ifield, Jfield);
    
          for(UInt i=1; i <= loopGrid->getNumElements(); ++i)
          {
            if(loopGrid->getElements().getRowMapL(i).getOwned())
            {
              op->setElement(i);
    
              sizeRow = localMatrix.numIndex_row();
              sizeCol = localMatrix.numIndex_col();
    
              indexRow.resize(sizeRow);
              indexCol.resize(sizeCol);

              localMatrix.indexG_row(indexRow);
              localMatrix.indexG_col(indexCol);
              mat = localMatrix.matrix();
    
              tpetraManip.tpetraCrsInsert(matrix, indexRow, indexCol, mat);
            }
          }  //End cycle elements
     
        }  //End cycle dofs
      }
    }
  }
  
  //Matrix assemble
  matrix->fillComplete(tpetraColMap, tpetraRowMap);
}


#endif
