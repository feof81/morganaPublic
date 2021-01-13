/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FESTATICFIELD3DGLOBALMANIP_HPP
#define FESTATICFIELD3DGLOBALMANIP_HPP

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "exprtk.hpp"

#include "pMapItem.h"
#include "pMapItemShare.h"
#include "pVectGlobalManip.hpp"

#include "traitsBasic.h"
#include "ggmGridSearch.hpp"
#include "morganaFields.hpp"
#include "feStaticField3d.hpp"

using namespace std;


//_________________________________________________________________________________________________
// SERVO-SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Static Field 3d Global Manipulator - empty general class  */
template<typename PMAPTYPE, typename FETYPE, typename DOFTYPE, dms3d_order ORDER = dms3d_vectMajor, dms3d_mode MODE = dms3d_allMode>
class fsf3dGlobalManip;



//_________________________________________________________________________________________________
// PMAPITEM SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Static Field 3d Global Manipulator - \c pMapItem specialization */
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
class fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem                                   PMAPTYPE;
    typedef dofMapStatic3d_options                     OPTIONS;
    typedef feStaticField3d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef typename FIELD::MESH3D                     MESH3D;
    typedef typename FIELD::CONNECT3D                  CONNECT3D;
    typedef typename FIELD::DOFVECT                    DOFVECT;
    typedef typename FIELD::GEOSHAPE                   GEOSHAPE;
    
    typedef typename GEOSHAPE::GEOBSHAPE            GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE          GEOSHAPE1D;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef mesh1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef connect1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE> CONNECT1D;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fsf3dGlobalManip();
    fsf3dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fsf3dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    bool initilize(const sArray<string>      & exprString,
                   const Teuchos::RCP<FIELD> & Field);
    
    bool initilize(const sArray<string> & exprString,
                         FIELD          & Field);
    
    bool initilize(const sArray<string>      & exprString,
                   const sVect<string>       & symbs,
                         sVect<Real>           symbsVals,
                   const Teuchos::RCP<FIELD> & Field);
    
    bool initilize(const sArray<string> & exprString,
                   const sVect<string>  & symbs,
                         sVect<Real>      symbsVals,
                         FIELD          & Field);
    
    //ATTENTION: requires that the two meshes has the same gids for each geometric entity: i.e a given face should be gid-numbered in the same way
    void changeGrid(const Teuchos::RCP<MESH3D>    & newGrid,
                    const Teuchos::RCP<CONNECT3D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH3D    & newGrid,
                    CONNECT3D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;

    bool matchGrid(const Teuchos::RCP<MESH3D>    & newGrid3d,
                   const Teuchos::RCP<CONNECT3D> & newConnect3d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH3D    & newGrid3d,
                   CONNECT3D & newConnect3d,
                   OPTIONS   & options,
                   FIELD     & field) const;
    
    Real evalMax(const Teuchos::RCP<FIELD> & Field,
                 const sVect<point3d>      & evalNodes,
                 const set<UInt>           & activeIds);
    
    Real evalMax(      FIELD          & Field,
                 const sVect<point3d> & evalNodes,
                 const set<UInt>      & activeIds);
    //@}
    
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH3D>             & NewGrid,
                            const Teuchos::RCP<CONNECT3D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH3D>    & NewGrid,
                            const Teuchos::RCP<CONNECT3D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH3D>             & NewGrid,
                            const Teuchos::RCP<CONNECT3D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH3D>    & NewGrid,
                            const Teuchos::RCP<CONNECT3D> & NewConnect,
                                  FIELD                   & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fsf3dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fsf3dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fsf3dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const Teuchos::RCP<FIELD> & Field)
{
  return(initilize(exprString,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
                FIELD          & Field)
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH3D    MESH3D;
  typedef typename FETYPE::DOFCARD  DOFCARD;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  //Alloc------------------------------------------------------------
  bool logic = true;
  UInt dofId;
  Real eval;
  Real x,y,z;
  DOFTYPE V;
  point3d Y, P;
  DOFCARD dofCard;
  sVect<point3d> points;
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<MESH3D> grid3d = Field.getMesh3d();
  geoMapInterface<GEOSHAPE> geoInterface;
  
  //EXPRTK startup--------------------------------------------------
  typedef exprtk::symbol_table<Real> SYMBOLTABLE;
  typedef exprtk::expression<Real>   EXPRESSION;
  typedef exprtk::parser<Real>       PARSER;
   
  SYMBOLTABLE symbol_table;
  symbol_table.add_variable("x",x);
  symbol_table.add_variable("y",y);
  symbol_table.add_variable("z",z);
  symbol_table.add_constants();
  
  
  sArray<Teuchos::RCP<EXPRESSION> > expressions(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  sArray<Teuchos::RCP<PARSER> >     parsers(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  
  for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)
  {
    for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
    {
      expressions(I,J) = Teuchos::rcp(new EXPRESSION);
      expressions(I,J)->register_symbol_table(symbol_table);
      
      parsers(I,J) = Teuchos::rcp(new PARSER);
      logic = logic & parsers(I,J)->compile(exprString(I,J), *expressions(I,J));
      
      assert(logic);
    }
  }
  
  //Compute----------------------------------------------------------
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    assert(grid3d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid3d->getElementNodesL(i);
     
      //Loop on the local basis
      for(UInt j=1; j <= FETYPE::numBasis; ++j)
      {
        Y = FETYPE::getRefNode(j);
        P = geoInterface.getPosition(points,Y);

        //Evaluate the dof
        for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)  
        {
          for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
          {
            x = P.getX();
            y = P.getY();
            z = P.getZ();

            eval = expressions(I,J)->value();
            traitsBasic<DOFTYPE>::setIJ(I,J,eval,V);
          }
        }

        //Map the dof 
        dofCard = FETYPE::getDofCard(j);
        dofCard.setLocalElId(i);

        dofId = Field.getDofMapper().mapDofL(dofCard);
        Field.setDofL(dofId,V);
      }
      //End loop on the local basis      
    }
  }
  
  return(logic);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const sVect<string>       & symbs,
                sVect<Real>           symbsVals,
          const Teuchos::RCP<FIELD> & Field)
{
  return(initilize(exprString,symbs,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
          const sVect<string>  & symbs,
                sVect<Real>      symbsVals,
                FIELD          & Field)
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  assert(symbs.size() == symbsVals.size());
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH3D    MESH3D;
  typedef typename FETYPE::DOFCARD  DOFCARD;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  //Alloc------------------------------------------------------------
  bool logic = true;
  UInt dofId;
  Real eval;
  Real x,y,z;
  DOFTYPE V;
  point3d Y, P;
  DOFCARD dofCard;
  sVect<point3d> points;
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<MESH3D> grid3d = Field.getMesh3d();
  geoMapInterface<GEOSHAPE> geoInterface;
  
  //EXPRTK startup--------------------------------------------------
  typedef exprtk::symbol_table<Real> SYMBOLTABLE;
  typedef exprtk::expression<Real>   EXPRESSION;
  typedef exprtk::parser<Real>       PARSER;
   
  SYMBOLTABLE symbol_table;
  symbol_table.add_variable("x",x);
  symbol_table.add_variable("y",y);
  symbol_table.add_variable("z",z);
  
  for(UInt i=1; i <= symbs.size(); ++i)
  { symbol_table.add_variable(symbs(i),symbsVals(i)); }
  
  symbol_table.add_constants();
  
  
  sArray<Teuchos::RCP<EXPRESSION> > expressions(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  sArray<Teuchos::RCP<PARSER> >     parsers(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  
  for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)
  {
    for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
    {
      expressions(I,J) = Teuchos::rcp(new EXPRESSION);
      expressions(I,J)->register_symbol_table(symbol_table);
      
      parsers(I,J) = Teuchos::rcp(new PARSER);
      logic = logic & parsers(I,J)->compile(exprString(I,J), *expressions(I,J));
      
      assert(logic);
    }
  }
  
  //Compute----------------------------------------------------------
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    assert(grid3d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid3d->getElementNodesL(i);
     
      //Loop on the local basis
      for(UInt j=1; j <= FETYPE::numBasis; ++j)
      {
        Y = FETYPE::getRefNode(j);
        P = geoInterface.getPosition(points,Y);

        //Evaluate the dof
        for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)  
        {
          for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
          {
            x = P.getX();
            y = P.getY();
            z = P.getZ();

            eval = expressions(I,J)->value();
            traitsBasic<DOFTYPE>::setIJ(I,J,eval,V);
          }
        }

        //Map the dof 
        dofCard = FETYPE::getDofCard(j);
        dofCard.setLocalElId(i);

        dofId = Field.getDofMapper().mapDofL(dofCard);
        Field.setDofL(dofId,V);
      }
      //End loop on the local basis      
    }
  }
  
  return(logic);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH3D>    & newGrid,
           const Teuchos::RCP<CONNECT3D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField3d;
  newField3d.setCommunicator(commDev);
  newField3d.setGeometry(newGrid,newConnect);
  newField3d.setOptions(options);
  newField3d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField3d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField3d.setDofVect(oldDofVect);
  *Field = newField3d;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH3D    & newGrid,
           CONNECT3D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField3d;
  newField3d.setCommunicator(commDev);
  newField3d.setGeometry(newGrid,newConnect);
  newField3d.setOptions(options);
  newField3d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField3d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField3d.setDofVect(oldDofVect);
  Field = newField3d;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH3D>    & newGrid3d,
          const Teuchos::RCP<CONNECT3D> & newConnect3d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  assert(commDevLoaded);
  
  return( matchGrid(*newGrid3d,
                    *newConnect3d,
                    *options,
                    *field) );
}
    
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH3D    & newGrid3d,
          CONNECT3D & newConnect3d,
          OPTIONS   & options,
          FIELD     & field) const
{
  //Typdefs----------------------------------------------------------------------------------------
  typedef ggmGridSearch<MESH2D,MESH2D>  FACESEARCH;
  typedef ggmGridSearch<MESH1D,MESH1D>  EDGESEARCH;
  
  typedef typename FACESEARCH::OUTVECT FACEVECT;
  typedef typename EDGESEARCH::OUTVECT EDGEVECT;
  
  typedef pMap<PMAPTYPE> FACEMAP;
  typedef pMap<PMAPTYPE> EDGEMAP;
  
  typedef typename FIELD::OPTIONS FIELD_OPTIONS;
  
  //Assert-----------------------------------------------------------------------------------------
  assert(commDevLoaded);
  
  //Download data----------------------------------------------------------------------------------
  Teuchos::RCP<MESH3D>    oldGrid3d(new MESH3D(*field.getMesh3d()));
  Teuchos::RCP<CONNECT3D> oldConnect3d(new CONNECT3D(*field.getMeshConnect()));
  
  //Alloc------------------------------------------------------------------------------------------
  bool success = true;
  sVect<point3d> inLocCoords;
  
  //Map the faces----------------------------------------------------------------------------------
  Teuchos::RCP<MESH2D> oldFaceGrid2d(new MESH2D(oldGrid3d->getNodes(),
                                                oldGrid3d->getIsVertex(),
                                                oldGrid3d->getFaces()));
  
  Teuchos::RCP<MESH2D> newFaceGrid2d(new MESH2D(newGrid3d.getNodes(),
                                                newGrid3d.getIsVertex(),
                                                newGrid3d.getFaces()));
  
  Teuchos::RCP<CONNECT2D> oldFaceConnect2d(new CONNECT2D(commDev));
  oldFaceConnect2d->setMesh2d(oldFaceGrid2d);
  
  Teuchos::RCP<CONNECT2D> newFaceConnect2d(new CONNECT2D(commDev));
  newFaceConnect2d->setMesh2d(newFaceGrid2d);
  
  FACESEARCH faceMatcher(commDev);
  faceMatcher.setMesh(oldFaceGrid2d, //tgt
                      oldFaceConnect2d,
                      newFaceGrid2d, //src
                      newFaceConnect2d);
  
  FACEVECT faceVect = faceMatcher.search(inLocCoords);
  FACEMAP  faceMap(faceVect.size());
  
  for(UInt i=1; i <= faceVect.size(); ++i)
  {
    success   = success && faceVect(i).getIsNested();
    faceMap(i) = faceVect(i).getElMap();
  }
  
  if(!success) { return(success); }
  
  assert(oldGrid3d->getNumFaces() == faceMap.size());
  
  //Map the edges----------------------------------------------------------------------------------
  Teuchos::RCP<MESH1D> oldEdgeGrid1d(new MESH1D(oldGrid3d->getNodes(),
                                                oldGrid3d->getIsVertex(),
                                                oldGrid3d->getEdges()));
  
  Teuchos::RCP<MESH1D> newEdgeGrid1d(new MESH1D(newGrid3d.getNodes(),
                                                newGrid3d.getIsVertex(),
                                                newGrid3d.getEdges()));
  
  Teuchos::RCP<CONNECT1D> oldEdgeConnect1d(new CONNECT1D(commDev));
  oldEdgeConnect1d->setMesh1d(oldEdgeGrid1d);
  
  Teuchos::RCP<CONNECT1D> newEdgeConnect1d(new CONNECT1D(commDev));
  newEdgeConnect1d->setMesh1d(newEdgeGrid1d);
  
  EDGESEARCH edgeMatcher(commDev);
  edgeMatcher.setMesh(oldEdgeGrid1d, //tgt
                      oldEdgeConnect1d,
                      newEdgeGrid1d, //src
                      newEdgeConnect1d);
  
  EDGEVECT edgeVect = edgeMatcher.search(inLocCoords);
  EDGEMAP  edgeMap(edgeVect.size());
  
  for(UInt i=1; i <= edgeVect.size(); ++i)
  {
    success   = success && edgeVect(i).getIsNested();
    edgeMap(i) = edgeVect(i).getElMap();
  }
  
  if(!success) { return(success); }
  
  assert(oldGrid3d->getNumEdges() == edgeMap.size());
  
  //Change the gids--------------------------------------------------------------------------------  
  for(UInt f=1; f <= oldGrid3d->getNumFaces(); ++f)
  { oldGrid3d->getFaces().getRowMapL(f).setGid(faceMap(f).getGid()); }
  
  for(UInt d=1; d <= oldGrid3d->getNumEdges(); ++d)
  { oldGrid3d->getEdges().getRowMapL(d).setGid(edgeMap(d).getGid()); }
  
  oldConnect3d->getVertexToEdge().setColMap(oldGrid3d->getEdges().getRowMap());
  oldConnect3d->getFaceToElement().setRowMap(oldGrid3d->getFaces().getRowMap());
  oldConnect3d->getElementToEdge().setColMap(oldGrid3d->getEdges().getRowMap());
  oldConnect3d->getElementToFace().setColMap(oldGrid3d->getFaces().getRowMap());
  
  //Build old field--------------------------------------------------------------------------------
  Teuchos::RCP<FIELD_OPTIONS> oldOptions(new FIELD_OPTIONS(field.getOptions())); 
  
  Teuchos::RCP<FIELD> oldField(new FIELD);
  oldField->setCommunicator(commDev);
  oldField->setGeometry(oldGrid3d,oldConnect3d);
  oldField->setFeCards(field.getFeCards());
  oldField->setOptions(oldOptions);
  oldField->startup();
  
  for(UInt lid=1; lid <= oldField->getDofVect().size(); ++lid)
  { oldField->setDofL(lid, field.getDofL(lid)); }
  
  //Change the maps--------------------------------------------------------------------------------
  changeGrid(newGrid3d,
             newConnect3d,
             options,
            *oldField);
  
  field = *oldField;
  
  return(success);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
Real
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
evalMax(const Teuchos::RCP<FIELD> & Field,
        const sVect<point3d>      & evalNodes,
        const set<UInt>           & activeIds)
{
  return(evalMax(*Field,evalNodes,activeIds));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
Real
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
evalMax(      FIELD          & Field,
        const sVect<point3d> & evalNodes,
        const set<UInt>      & activeIds)
{
  //Assert
  assert(evalNodes.size() >= 1);
  
  //Typedefs
  typedef typename FIELD::OUTTYPE OUTTYPE;
  
  //Alloc
  Real valMax = 0.0;
  UInt geoId;
  OUTTYPE outVal;
  Teuchos::RCP<MESH3D> grid3d = Field.getMesh();
  
  //Find max
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    geoId = grid3d->getElementL(i).getGeoId();
    
    if(activeIds.count(geoId) == 1)
    {
      for(UInt j=1; j <= evalNodes.size(); ++j)
      {
        Field.evalL(i,evalNodes(j),outVal);
        valMax = std::max(traitsBasic<OUTTYPE>::norm(outVal), valMax);
      }
    }
  }
  
  //Comunication
  boost::mpi::maximum<Real> op;
  all_reduce(*commDev,valMax,valMax,op);
  
  return(valMax);
}


//_________________________________________________________________________________________________
// COMMUNICATOR MANIPULATIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH3D>             & NewGrid,
                   const Teuchos::RCP<CONNECT3D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH3D>    & NewGrid,
                   const Teuchos::RCP<CONNECT3D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,DOFTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,DOFTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Reduce comm
    manipDofVect.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH3D>             & NewGrid,
                   const Teuchos::RCP<CONNECT3D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH3D>    & NewGrid,
                   const Teuchos::RCP<CONNECT3D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,DOFTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,DOFTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Reduce comm
    manipDofVect.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}



//_________________________________________________________________________________________________
// PMAPITEMSHARE SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Static Field 3d Global Manipulator - \c pMapItemShare specialization */
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
class fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare                              PMAPTYPE;
    typedef feStaticField3d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef pMap<pMapItemSendRecv>                     SENDRECV;
    typedef dofMapStatic3d_options                     OPTIONS;
    typedef typename FIELD::MESH3D                     MESH3D;
    typedef typename FIELD::CONNECT3D                  CONNECT3D;
    typedef typename FIELD::DOFVECT                    DOFVECT;
    typedef typename FIELD::GEOSHAPE                   GEOSHAPE;
    
    typedef pVectGlobalManip<DOFTYPE,pMapItemShare>  PVGLOBMANIP;
    typedef typename PVGLOBMANIP::PVPS               PVPS;
    typedef typename PVGLOBMANIP::PVUR               PVUR;
    
    typedef typename GEOSHAPE::GEOBSHAPE            GEOSHAPE2D;
    typedef typename GEOSHAPE2D::GEOBSHAPE          GEOSHAPE1D;
    typedef mesh2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE>    MESH2D;
    typedef mesh1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE>    MESH1D;
    typedef connect2d<GEOSHAPE2D,PMAPTYPE,PMAPTYPE> CONNECT2D;
    typedef connect1d<GEOSHAPE1D,PMAPTYPE,PMAPTYPE> CONNECT1D;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fsf3dGlobalManip();
    fsf3dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fsf3dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    bool initilize(const sArray<string>      & exprString,
                   const Teuchos::RCP<FIELD> & Field);
    
    bool initilize(const sArray<string> & exprString,
                         FIELD          & Field);
    
    bool initilize(const sArray<string>      & exprString,
                   const sVect<string>       & symbs,
                   const sVect<Real>           symbsVals,
                   const Teuchos::RCP<FIELD> & Field);
    
    bool initilize(const sArray<string> & exprString,
                   const sVect<string>  & symbs,
                   const sVect<Real>      symbsVals,
                         FIELD          & Field);
    
    //ATTENTION: requires that the two meshes has the same gids for each geometric entity:
    // i.e a given face should be gid-numbered in the same way
    //ATTENTION: the feCards must be loaded after the function
    void changeGrid(const Teuchos::RCP<MESH3D>    & newGrid,
                    const Teuchos::RCP<CONNECT3D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH3D    & newGrid,
                    CONNECT3D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;

    bool matchGrid(const Teuchos::RCP<MESH3D>    & newGrid3d,
                   const Teuchos::RCP<CONNECT3D> & newConnect3d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH3D    & newGrid3d,
                   CONNECT3D & newConnect3d,
                   OPTIONS   & options,
                   FIELD     & field) const;
    
    void createSendRecvMap(const Teuchos::RCP<FIELD>    & Field,
                                 Teuchos::RCP<SENDRECV> & mapSend,
                                 Teuchos::RCP<SENDRECV> & mapRecv) const;
    
    void createSendRecvMap(const FIELD    & Field,
                                 SENDRECV & mapSend,
                                 SENDRECV & mapRecv) const;
    
    void updateData(const Teuchos::RCP<FIELD>    & Field,
                    const Teuchos::RCP<SENDRECV> & mapSend) const;
    
    void updateData(FIELD    & Field,
                    SENDRECV & mapSend) const;

    void updateData(const Teuchos::RCP<FIELD>    & Field,
                    const Teuchos::RCP<SENDRECV> & mapSend,
                    const Teuchos::RCP<SENDRECV> & mapRecv) const;
    
    void updateData(FIELD    & Field,
                    SENDRECV & mapSend,
                    SENDRECV & mapRecv) const;
    
    Real evalMax(const Teuchos::RCP<FIELD> & Field,
                 const sVect<point3d>      & evalNodes,
                 const set<UInt>           & activeIds);
    
    Real evalMax(      FIELD          & Field,
                 const sVect<point3d> & evalNodes,
                 const set<UInt>      & activeIds);
    //@}
    
    /*! @name Non-pending update functions */ //@{
  public:
    void updateDataI(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
    
    void updateDataO(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
     
    void updateDataI(FIELD          & Field,
                     SENDRECV       & mapSend,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
   
    void updateDataO(FIELD          & Field,
                     SENDRECV       & mapSend,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
     
    void updateDataI(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                     const Teuchos::RCP<SENDRECV> & mapRecv,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
    
    void updateDataO(const Teuchos::RCP<FIELD>    & Field,
                     const Teuchos::RCP<SENDRECV> & mapSend,
                     const Teuchos::RCP<SENDRECV> & mapRecv,
                           sVect<DOFVECT>         & bufSegments,
                           sVect<PVPS>            & sendPvps,
                           sVect<PVPS>            & recvPvps,
                     const UInt                   & channel = 2) const;
     
    void updateDataI(FIELD          & Field,
                     SENDRECV       & mapSend,
                     SENDRECV       & mapRecv,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
   
    void updateDataO(FIELD          & Field,
                     SENDRECV       & mapSend,
                     SENDRECV       & mapRecv,
                     sVect<DOFVECT> & bufSegments,
                     sVect<PVPS>    & sendPvps,
                     sVect<PVPS>    & recvPvps,
                     const UInt     & channel = 2) const;
    //@}
     
    /*! @name Recursive update functions */ //@{
  public:
    void updateDataRR(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
      
    void updateDataRR(      FIELD    & Field,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRI(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRI(      FIELD    & Field,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRO(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRO(      FIELD    & Field,
                      const SENDRECV & mapSend,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
     
    void updateDataRR(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
      
    void updateDataRR(      FIELD    & Field,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRI(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRI(      FIELD    & Field,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
      
    void updateDataRO(const Teuchos::RCP<FIELD>          & Field,
                      const Teuchos::RCP<const SENDRECV> & mapSend,
                      const Teuchos::RCP<const SENDRECV> & mapRecv,
                            PVUR & commBuffer,
                      const UInt & channel = 2) const;
     
    void updateDataRO(      FIELD    & Field,
                      const SENDRECV & mapSend,
                      const SENDRECV & mapRecv,
                            PVUR     & commBuffer,
                      const UInt     & channel = 2) const;
    //@}
      
    /*! @name Communicator manipulations */ //@{
  public:
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH3D>             & NewGrid,
                            const Teuchos::RCP<CONNECT3D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH3D>    & NewGrid,
                            const Teuchos::RCP<CONNECT3D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH3D>             & NewGrid,
                            const Teuchos::RCP<CONNECT3D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                     & isActive,
                            const communicator             & OldCommDev,
                            const FIELD                    & OldField,
                                  communicator             & NewCommDev,
                            const Teuchos::RCP<MESH3D>     & NewGrid,
                            const Teuchos::RCP<CONNECT3D>  & NewConnect,
                                  FIELD                    & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fsf3dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fsf3dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fsf3dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const Teuchos::RCP<FIELD> & Field)
{
  return(initilize(exprString,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
                FIELD          & Field)
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH3D    MESH3D;
  typedef typename FETYPE::DOFCARD  DOFCARD;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  //Alloc------------------------------------------------------------
  bool logic = true;
  UInt dofId;
  Real eval;
  Real x,y,z;
  DOFTYPE V;
  point3d Y, P;
  DOFCARD dofCard;
  sVect<point3d> points;
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<MESH3D> grid3d = Field.getMesh3d();
  geoMapInterface<GEOSHAPE> geoInterface;
  
  //EXPRTK startup--------------------------------------------------
  typedef exprtk::symbol_table<Real> SYMBOLTABLE;
  typedef exprtk::expression<Real>   EXPRESSION;
  typedef exprtk::parser<Real>       PARSER;
   
  SYMBOLTABLE symbol_table;
  symbol_table.add_variable("x",x);
  symbol_table.add_variable("y",y);
  symbol_table.add_variable("z",z);
  symbol_table.add_constants();
  
  
  sArray<Teuchos::RCP<EXPRESSION> > expressions(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  sArray<Teuchos::RCP<PARSER> >     parsers(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  
  for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)
  {
    for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
    {
      expressions(I,J) = Teuchos::rcp(new EXPRESSION);
      expressions(I,J)->register_symbol_table(symbol_table);
      
      parsers(I,J) = Teuchos::rcp(new PARSER);
      logic = logic & parsers(I,J)->compile(exprString(I,J), *expressions(I,J));
      
      assert(logic);
    }
  }
  
  //Compute----------------------------------------------------------
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    assert(grid3d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid3d->getElementNodesL(i);
     
      //Loop on the local basis
      for(UInt j=1; j <= FETYPE::numBasis; ++j)
      {
        Y = FETYPE::getRefNode(j);
        P = geoInterface.getPosition(points,Y);

        //Evaluate the dof
        for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)  
        {
          for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
          {
            x = P.getX();
            y = P.getY();
            z = P.getZ();

            eval = expressions(I,J)->value();
            traitsBasic<DOFTYPE>::setIJ(I,J,eval,V);
          }
        }

        //Map the dof 
        dofCard = FETYPE::getDofCard(j);
        dofCard.setLocalElId(i);

        dofId = Field.getDofMapper().mapDofL(dofCard);
        Field.setDofL(dofId,V);
      }
      //End loop on the local basis      
    }
  }
  
  return(logic);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const sVect<string>       & symbs,
                sVect<Real>           symbsVals,
          const Teuchos::RCP<FIELD> & Field)
{
  return(initilize(exprString,symbs,symbsVals,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
          const sVect<string>  & symbs,
                sVect<Real>      symbsVals,
                FIELD          & Field)
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  assert(symbs.size() == symbsVals.size());
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH3D    MESH3D;
  typedef typename FETYPE::DOFCARD  DOFCARD;
  typedef typename FETYPE::GEOSHAPE GEOSHAPE;
  
  //Alloc------------------------------------------------------------
  bool logic = true;
  UInt dofId;
  Real eval;
  Real x,y,z;
  DOFTYPE V;
  point3d Y, P;
  DOFCARD dofCard;
  sVect<point3d> points;
  
  //Alloc structures-------------------------------------------------
  Teuchos::RCP<MESH3D> grid3d = Field.getMesh3d();
  geoMapInterface<GEOSHAPE> geoInterface;
  
  //EXPRTK startup--------------------------------------------------
  typedef exprtk::symbol_table<Real> SYMBOLTABLE;
  typedef exprtk::expression<Real>   EXPRESSION;
  typedef exprtk::parser<Real>       PARSER;
   
  Real a;
  SYMBOLTABLE symbol_table;
  symbol_table.add_variable("x",x);
  symbol_table.add_variable("y",y);
  symbol_table.add_variable("z",z);
  
  for(UInt i=1; i <= symbs.size(); ++i)
  {
    a = symbsVals(i);
    symbol_table.add_variable(symbs(i),a);
  }
  
  symbol_table.add_constants();
  
  
  sArray<Teuchos::RCP<EXPRESSION> > expressions(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  sArray<Teuchos::RCP<PARSER> >     parsers(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  
  for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)
  {
    for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
    {
      expressions(I,J) = Teuchos::rcp(new EXPRESSION);
      expressions(I,J)->register_symbol_table(symbol_table);
      
      parsers(I,J) = Teuchos::rcp(new PARSER);
      logic = logic & parsers(I,J)->compile(exprString(I,J), *expressions(I,J));
      
      assert(logic);
    }
  }
  
  //Compute----------------------------------------------------------
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    assert(grid3d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid3d->getElementNodesL(i);
     
      //Loop on the local basis
      for(UInt j=1; j <= FETYPE::numBasis; ++j)
      {
        Y = FETYPE::getRefNode(j);
        P = geoInterface.getPosition(points,Y);

        //Evaluate the dof
        for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)  
        {
          for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
          {
            x = P.getX();
            y = P.getY();
            z = P.getZ();

            eval = expressions(I,J)->value();
           traitsBasic<DOFTYPE>::setIJ(I,J,eval,V);
          }
        }

        //Map the dof 
        dofCard = FETYPE::getDofCard(j);
        dofCard.setLocalElId(i);

        dofId = Field.getDofMapper().mapDofL(dofCard);
        Field.setDofL(dofId,V);
      }
      //End loop on the local basis      
    }
  }
  
  return(logic);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH3D>    & newGrid,
           const Teuchos::RCP<CONNECT3D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField3d;
  newField3d.setCommunicator(commDev);
  newField3d.setGeometry(newGrid,newConnect);
  newField3d.setOptions(options);
  newField3d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField3d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField3d.setDofVect(oldDofVect);
  *Field = newField3d;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH3D    & newGrid,
           CONNECT3D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField3d;
  newField3d.setCommunicator(commDev);
  newField3d.setGeometry(newGrid,newConnect);
  newField3d.setOptions(options);
  newField3d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField3d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField3d.setDofVect(oldDofVect);
  Field = newField3d;
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH3D>    & newGrid3d,
          const Teuchos::RCP<CONNECT3D> & newConnect3d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  assert(commDevLoaded);
  
  return( matchGrid(*newGrid3d,
                    *newConnect3d,
                    *options,
                    *field) );
}
    
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
bool
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH3D    & newGrid3d,
          CONNECT3D & newConnect3d,
          OPTIONS   & options,
          FIELD     & field) const
{
  //Typdefs----------------------------------------------------------------------------------------
  typedef ggmGridSearch<MESH2D,MESH2D>  FACESEARCH;
  typedef ggmGridSearch<MESH1D,MESH1D>  EDGESEARCH;
  
  typedef typename FACESEARCH::OUTVECT FACEVECT;
  typedef typename EDGESEARCH::OUTVECT EDGEVECT;
  
  typedef pMap<PMAPTYPE> FACEMAP;
  typedef pMap<PMAPTYPE> EDGEMAP;
  
  typedef typename FIELD::OPTIONS FIELD_OPTIONS;
  
  //Assert-----------------------------------------------------------------------------------------
  assert(commDevLoaded);
  
  //Download data----------------------------------------------------------------------------------
  Teuchos::RCP<MESH3D>    oldGrid3d(new MESH3D(*field.getMesh3d()));
  Teuchos::RCP<CONNECT3D> oldConnect3d(new CONNECT3D(*field.getMeshConnect()));
  
  //Alloc------------------------------------------------------------------------------------------
  bool success = true;
  sVect<point3d> inLocCoords;
  
  //Map the faces----------------------------------------------------------------------------------
  Teuchos::RCP<MESH2D> oldFaceGrid2d(new MESH2D(oldGrid3d->getNodes(),
                                                oldGrid3d->getIsVertex(),
                                                oldGrid3d->getFaces()));
  
  Teuchos::RCP<MESH2D> newFaceGrid2d(new MESH2D(newGrid3d.getNodes(),
                                                newGrid3d.getIsVertex(),
                                                newGrid3d.getFaces()));
  
  Teuchos::RCP<CONNECT2D> oldFaceConnect2d(new CONNECT2D(commDev));
  oldFaceConnect2d->setMesh2d(oldFaceGrid2d);
  
  Teuchos::RCP<CONNECT2D> newFaceConnect2d(new CONNECT2D(commDev));
  newFaceConnect2d->setMesh2d(newFaceGrid2d);
  
  FACESEARCH faceMatcher(commDev);
  faceMatcher.setMesh(oldFaceGrid2d, //tgt
                      oldFaceConnect2d,
                      newFaceGrid2d, //src
                      newFaceConnect2d);
  
  FACEVECT faceVect = faceMatcher.search(inLocCoords);
  FACEMAP  faceMap(faceVect.size());
  
  for(UInt i=1; i <= faceVect.size(); ++i)
  {
    success   = success && faceVect(i).getIsNested();
    faceMap(i) = faceVect(i).getElMap();
  }
  
  if(!success) { return(success); }
  
  assert(oldGrid3d->getNumFaces() == faceMap.size());
  
  //Map the edges----------------------------------------------------------------------------------
  Teuchos::RCP<MESH1D> oldEdgeGrid1d(new MESH1D(oldGrid3d->getNodes(),
                                                oldGrid3d->getIsVertex(),
                                                oldGrid3d->getEdges()));
  
  Teuchos::RCP<MESH1D> newEdgeGrid1d(new MESH1D(newGrid3d.getNodes(),
                                                newGrid3d.getIsVertex(),
                                                newGrid3d.getEdges()));
  
  Teuchos::RCP<CONNECT1D> oldEdgeConnect1d(new CONNECT1D(commDev));
  oldEdgeConnect1d->setMesh1d(oldEdgeGrid1d);
  
  Teuchos::RCP<CONNECT1D> newEdgeConnect1d(new CONNECT1D(commDev));
  newEdgeConnect1d->setMesh1d(newEdgeGrid1d);
  
  EDGESEARCH edgeMatcher(commDev);
  edgeMatcher.setMesh(oldEdgeGrid1d, //tgt
                      oldEdgeConnect1d,
                      newEdgeGrid1d, //src
                      newEdgeConnect1d);
  
  EDGEVECT edgeVect = edgeMatcher.search(inLocCoords);
  EDGEMAP  edgeMap(edgeVect.size());
  
  for(UInt i=1; i <= edgeVect.size(); ++i)
  {
    success   = success && edgeVect(i).getIsNested();
    edgeMap(i) = edgeVect(i).getElMap();
  }
  
  if(!success) { return(success); }
  
  assert(oldGrid3d->getNumEdges() == edgeMap.size());
  
  //Change the gids--------------------------------------------------------------------------------  
  for(UInt f=1; f <= oldGrid3d->getNumFaces(); ++f)
  { oldGrid3d->getFaces().getRowMapL(f).setGid(faceMap(f).getGid()); }
  
  for(UInt d=1; d <= oldGrid3d->getNumEdges(); ++d)
  { oldGrid3d->getEdges().getRowMapL(d).setGid(edgeMap(d).getGid()); }
  
  oldConnect3d->getVertexToEdge().setColMap(oldGrid3d->getEdges().getRowMap());
  oldConnect3d->getFaceToElement().setRowMap(oldGrid3d->getFaces().getRowMap());
  oldConnect3d->getElementToEdge().setColMap(oldGrid3d->getEdges().getRowMap());
  oldConnect3d->getElementToFace().setColMap(oldGrid3d->getFaces().getRowMap());
  
  //Build old field--------------------------------------------------------------------------------
  Teuchos::RCP<FIELD_OPTIONS> oldOptions(new FIELD_OPTIONS(field.getOptions())); 
  
  Teuchos::RCP<FIELD> oldField(new FIELD);
  oldField->setCommunicator(commDev);
  oldField->setGeometry(oldGrid3d,oldConnect3d);
  oldField->setFeCards(field.getFeCards());
  oldField->setOptions(oldOptions);
  oldField->startup();
  
  for(UInt lid=1; lid <= oldField->getDofVect().size(); ++lid)
  { oldField->setDofL(lid, field.getDofL(lid)); }
  
  //Change the maps--------------------------------------------------------------------------------
  changeGrid(newGrid3d,
             newConnect3d,
             options,
            *oldField);
  
  field = *oldField;
  
  return(success);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const Teuchos::RCP<FIELD>    & Field,
                        Teuchos::RCP<SENDRECV> & mapSend,
                        Teuchos::RCP<SENDRECV> & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field->getDofVect().getMapRcp(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const FIELD    & Field,
                        SENDRECV & mapSend,
                        SENDRECV & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field.getDofVect().getMapRef(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(const Teuchos::RCP<FIELD>    & Field,
           const Teuchos::RCP<SENDRECV> & mapSend) const
{
  //Checks
  assert(commDevLoaded);
    
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field->getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,*mapSend);
  
  //Reload the dofs
  Field->setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(FIELD    & Field,
           SENDRECV & mapSend) const
{
  //Checks
  assert(commDevLoaded);
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field.getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,mapSend);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(const Teuchos::RCP<FIELD>    & Field,
           const Teuchos::RCP<SENDRECV> & mapSend,
           const Teuchos::RCP<SENDRECV> & mapRecv) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field->getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,*mapSend,*mapRecv);
  
  //Reload the dofs
  Field->setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(FIELD    & Field,
           SENDRECV & mapSend,
           SENDRECV & mapRecv) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field.getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,mapSend,mapRecv);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
Real
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
evalMax(const Teuchos::RCP<FIELD> & Field,
        const sVect<point3d>      & evalNodes,
        const set<UInt>           & activeIds)
{
  return(evalMax(*Field,evalNodes,activeIds));
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
Real
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
evalMax(      FIELD          & Field,
        const sVect<point3d> & evalNodes,
        const set<UInt>      & activeIds)
{
  //Assert
  assert(evalNodes.size() >= 1);
  
  //Typedefs
  typedef typename FIELD::OUTTYPE OUTTYPE;
  
  //Alloc
  Real valMax = 0.0;
  UInt geoId;
  OUTTYPE outVal;
  Teuchos::RCP<MESH3D> grid3d = Field.getMesh();
  
  //Find max
  for(UInt i=1; i <= grid3d->getNumElements(); ++i)
  {
    geoId = grid3d->getElementL(i).getGeoId();
    
    if(activeIds.count(geoId) == 1)
    {
      for(UInt j=1; j <= evalNodes.size(); ++j)
      {
	Field.evalL(i,evalNodes(j),outVal);
	valMax = std::max(traitsBasic<OUTTYPE>::norm(outVal), valMax);
      }
    }
  }
  
  //Comunication
  boost::mpi::maximum<Real> op;
  all_reduce(*commDev,valMax,valMax,op);
  
  return(valMax);
}


//_________________________________________________________________________________________________
// NON-PENDING UPDATE FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataI(*Field,
              *mapSend,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataO(*Field,
              *mapSend,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}
   
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>     
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(FIELD          & Field,
            SENDRECV       & mapSend,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataI(dofsVect,
                          mapSend,
                          bufSegments,
                          sendPvps,
                          recvPvps,
                          channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(FIELD          & Field,
            SENDRECV       & mapSend,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataO(dofsVect,
                          mapSend,
                          bufSegments,
                          sendPvps,
                          recvPvps,
                          channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
            const Teuchos::RCP<SENDRECV> & mapRecv,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataI(*Field,
              *mapSend,
              *mapRecv,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(const Teuchos::RCP<FIELD>    & Field,
            const Teuchos::RCP<SENDRECV> & mapSend,
            const Teuchos::RCP<SENDRECV> & mapRecv,
                  sVect<DOFVECT>         & bufSegments,
                  sVect<PVPS>            & sendPvps,
                  sVect<PVPS>            & recvPvps,
            const UInt                   & channel) const
{
  assert(commDevLoaded);
  updateDataO(*Field,
              *mapSend,
              *mapRecv,
               bufSegments,
               sendPvps,
               recvPvps,
               channel);
}
   
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>     
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(FIELD          & Field,
            SENDRECV       & mapSend,
            SENDRECV       & mapRecv,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataI(dofsVect,
                          mapSend,
                          mapRecv,
                          bufSegments,
                          sendPvps,
                          recvPvps,
                          channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(FIELD          & Field,
            SENDRECV       & mapSend,
            SENDRECV       & mapRecv,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Dofs Vect
  pVect<DOFTYPE,pMapItemShare> & dofsVect = Field.dofVect;
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataO(dofsVect,
                          mapSend,
                          mapRecv,
                          bufSegments,
                          sendPvps,
                          recvPvps,
                          channel);
}


//_________________________________________________________________________________________________
// RECURSIVE UPDATE FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRR(*Field,
               *mapSend,
                commBuffer,
                channel);
}
      
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(      FIELD    & Field,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRR(Field.dofVect,
                           mapSend,
                           commBuffer,
                           channel);
}
      
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRI(*Field,
               *mapSend,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(      FIELD    & Field,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRI(Field.dofVect,
                           mapSend,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
                    PVUR & commBuffer,
              const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRO(*Field,
               *mapSend,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(      FIELD    & Field,
             const SENDRECV & mapSend,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  assert(commDevLoaded);
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRO(Field.dofVect,
                           mapSend,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRR(*Field,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRR(      FIELD    & Field,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRR(Field.dofVect,
                           mapSend,
                           mapRecv,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRI(*Field,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRI(      FIELD    & Field,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRI(Field.dofVect,
                           mapSend,
                           mapRecv,
                           commBuffer,
                           channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(const Teuchos::RCP<FIELD>          & Field,
             const Teuchos::RCP<const SENDRECV> & mapSend,
             const Teuchos::RCP<const SENDRECV> & mapRecv,
                   PVUR & commBuffer,
             const UInt & channel) const
{
  assert(commDevLoaded);
  updateDataRO(*Field,
               *mapSend,
               *mapRecv,
                commBuffer,
                channel);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>  
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataRO(      FIELD    & Field,
             const SENDRECV & mapSend,
             const SENDRECV & mapRecv,
                   PVUR     & commBuffer,
             const UInt     & channel) const
{
  //Assert
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  
  assert(commDevLoaded);
  assert(cheker.check(mapSend,mapRecv));
  
  //Update the data
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateDataRO(Field.dofVect,
                           mapSend,
                           mapRecv,
                           commBuffer,
                           channel);
}


//_________________________________________________________________________________________________
// COMMUNICATOR MANIPULATIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH3D>             & NewGrid,
                   const Teuchos::RCP<CONNECT3D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  reduceCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH3D>    & NewGrid,
                   const Teuchos::RCP<CONNECT3D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,PMAPTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,PMAPTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() < OldCommDev.size());
    
    //Reduce comm
    manipDofVect.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.reduceCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH3D>             & NewGrid,
                   const Teuchos::RCP<CONNECT3D>          & NewConnect,
                         Teuchos::RCP<FIELD>              & NewField)
{
  expandCommunicator(isActive,
                    *OldCommDev,
                    *OldField,
                    *NewCommDev,
                     NewGrid,
                     NewConnect,
                    *NewField);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
void
fsf3dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH3D>    & NewGrid,
                   const Teuchos::RCP<CONNECT3D> & NewConnect,
                         FIELD                   & NewField)
{
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::DOFVECT DOFVECT;
  typedef typename FIELD::FECARDS FECARDS;
  typedef typename FIELD::FECARD  FECARD;
  
  typedef pVectGlobalManip<DOFTYPE,DOFTYPE>  MANIP_DOFVECT;
  typedef pVectGlobalManip<FECARD,DOFTYPE>   MANIP_FECARDS;
  
  //Manipulators-----------------------------------------------------
  DOFVECT newDofVect;
  FECARDS newFeCards;

  MANIP_DOFVECT manipDofVect;
  MANIP_FECARDS manipFeCards;
  
  OPTIONS options = OldField.getOptions();
  
  //Copy to new communicator-----------------------------------------
  if(isActive)
  {
    //Assert
    assert(NewCommDev.size() > OldCommDev.size());
    
    //Reduce comm
    manipDofVect.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getDofVect(),
                                    NewCommDev,
                                    newDofVect);
    
    manipFeCards.expandCommunicator(isActive,
                                    OldCommDev,
                                    OldField.getFeCards(),
                                    NewCommDev,
                                    newFeCards);
    
    //New Field
    NewField.setCommunicator(NewCommDev);
    NewField.setGeometry(NewGrid,NewConnect);
    NewField.setOptions(options);
    NewField.setFeCards(newFeCards);
    NewField.startup();
    NewField.setDofVect(newDofVect);
  }
}



//_________________________________________________________________________________________________
// SUMMARY CLASS
//-------------------------------------------------------------------------------------------------

/*! Interface class for the manipulation of Static Field 3d */
template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER = dms3d_vectMajor, dms3d_mode MODE = dms3d_allMode>
class feStaticField3dGlobalManip : public fsf3dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>
{
     /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE PMAPTYPE;
    typedef feStaticField3d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef typename FETYPE::GEOSHAPE  GEOSHAPE;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feStaticField3dGlobalManip();
    feStaticField3dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    feStaticField3dGlobalManip(communicator & CommDev);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
feStaticField3dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField3dGlobalManip() : fsf3dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>()
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
feStaticField3dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField3dGlobalManip(const Teuchos::RCP<communicator> & CommDev) : fsf3dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dms3d_order ORDER, dms3d_mode MODE>
feStaticField3dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField3dGlobalManip(communicator & CommDev) : fsf3dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 3>::returnValue);
}

#endif
