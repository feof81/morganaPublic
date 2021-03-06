/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FESTATICFIELD1DGLOBALMANIP_HPP
#define FESTATICFIELD1DGLOBALMANIP_HPP

#include "exprtk.hpp"

#include "pMapItem.h"
#include "pMapItemShare.h"

#include "morganaFields.hpp"
#include "feStaticField1d.hpp"

using namespace std;


//_________________________________________________________________________________________________
// SERVO-SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Static Field 1d Global Manipulator - empty general class  */
template<typename PMAPTYPE, typename FETYPE, typename DOFTYPE, dms1d_order ORDER = dms1d_vectMajor, dms1d_mode MODE = dms1d_allMode>
class fsf1dGlobalManip;



//_________________________________________________________________________________________________
// PMAPITEM SPECIALIZATION
//-------------------------------------------------------------------------------------------------

/*! Static Field 2d Global Manipulator - \c pMapItem specialization */
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
class fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItem                                   PMAPTYPE;
    typedef dofMapStatic1d_options                     OPTIONS;
    typedef feStaticField1d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef typename FIELD::MESH1D                     MESH1D;
    typedef typename FIELD::CONNECT1D                  CONNECT1D;
    typedef typename FIELD::DOFVECT                    DOFVECT;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fsf1dGlobalManip();
    fsf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fsf1dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    bool initilize(const sArray<string>      & exprString,
                   const Teuchos::RCP<FIELD> & Field) const;
    
    bool initilize(const sArray<string> & exprString,
                         FIELD          & Field) const;
    
    bool initilize(const sArray<string>      & exprString,
                   const sVect<string>       & symbs,
                         sVect<Real>           symbsVals,
                   const Teuchos::RCP<FIELD> & Field) const;
    
    bool initilize(const sArray<string> & exprString,
                   const sVect<string>  & symbs,
                         sVect<Real>      symbsVals,
                         FIELD          & Field) const;
    
    void changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
                    const Teuchos::RCP<CONNECT1D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH1D    & newGrid,
                    CONNECT1D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;

    bool matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
                   const Teuchos::RCP<CONNECT1D> & newConnect1d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH1D    & newGrid1d,
                   CONNECT1D & newConnect1d,
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
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH1D>    & NewGrid,
                            const Teuchos::RCP<CONNECT1D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH1D>    & NewGrid,
                            const Teuchos::RCP<CONNECT1D> & NewConnect,
                                  FIELD                   & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fsf1dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fsf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
fsf1dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const Teuchos::RCP<FIELD> & Field) const
{
  return(initilize(exprString,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
                FIELD          & Field) const
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH1D    MESH1D;
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
  Teuchos::RCP<MESH1D> grid1d = Field.getMesh1d();
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
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    assert(grid1d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid1d->getElementNodesL(i);
     
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const sVect<string>       & symbs,
                sVect<Real>           symbsVals,
          const Teuchos::RCP<FIELD> & Field) const
{
  return(initilize(exprString,symbs,symbsVals,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
          const sVect<string>  & symbs,
                sVect<Real>      symbsVals,
                FIELD          & Field) const
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  assert(symbs.size() == symbsVals.size());
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH1D    MESH1D;
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
  Teuchos::RCP<MESH1D> grid1d = Field.getMesh1d();
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
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    assert(grid1d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid1d->getElementNodesL(i);
     
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
           const Teuchos::RCP<CONNECT1D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  *Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH1D    & newGrid,
           CONNECT1D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
          const Teuchos::RCP<CONNECT1D> & newConnect1d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  changeGrid(newGrid1d,
             newConnect1d,
             options,
             field);
  
  return(true);
}
    
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH1D    & newGrid1d,
          CONNECT1D & newConnect1d,
          OPTIONS   & options,
          FIELD     & field) const
{
  changeGrid(newGrid1d,
             newConnect1d,
             options,
             field);
  
  return(true);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
Real
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
evalMax(const Teuchos::RCP<FIELD> & Field,
        const sVect<point3d>      & evalNodes,
        const set<UInt>           & activeIds)
{
  return(evalMax(*Field,evalNodes,activeIds));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
Real
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
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
  Teuchos::RCP<MESH1D> grid1d = Field.getMesh();
  
  //Find max
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    geoId = grid1d->getElementL(i).getGeoId();
    
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
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItem,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
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

/*! Static Field 2d Global Manipulator - \c pMapItemShare specialization */
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
class fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef pMapItemShare                              PMAPTYPE;
    typedef feStaticField1d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef pMap<pMapItemSendRecv>                     SENDRECV;
    typedef dofMapStatic1d_options                     OPTIONS;
    typedef typename FIELD::MESH1D                     MESH1D;
    typedef typename FIELD::CONNECT1D                  CONNECT1D;
    typedef typename FIELD::DOFVECT                    DOFVECT;
    typedef pVectGlobalManip<DOFTYPE,pMapItemShare>    PVGLOBMANIP;
    typedef typename PVGLOBMANIP::PVPS                 PVPS;
    typedef typename PVGLOBMANIP::PVUR                 PVUR;
    //@}
  
    /*! @name Internal data and links */ //@{
  public:
    bool commDevLoaded;
    Teuchos::RCP<communicator> commDev;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    fsf1dGlobalManip();
    fsf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    fsf1dGlobalManip(communicator & CommDev);
    void setCommunicator(const Teuchos::RCP<communicator> & CommDev);
    void setCommunicator(communicator & CommDev);
    //@}
    
    /*! @name Operator functions */ //@{
  public:
    bool initilize(const sArray<string>      & exprString,
                   const Teuchos::RCP<FIELD> & Field) const;
    
    bool initilize(const sArray<string> & exprString,
                         FIELD          & Field) const;
    
    bool initilize(const sArray<string>      & exprString,
                   const sVect<string>       & symbs,
                         sVect<Real>           symbsVals,
                   const Teuchos::RCP<FIELD> & Field) const;
    
    bool initilize(const sArray<string> & exprString,
                   const sVect<string>  & symbs,
                         sVect<Real>      symbsVals,
                         FIELD          & Field) const;
    
    void changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
                    const Teuchos::RCP<CONNECT1D> & newConnect,
                    const Teuchos::RCP<OPTIONS>   & options,
                    const Teuchos::RCP<FIELD>     & Field) const;
    
    void changeGrid(MESH1D    & newGrid,
                    CONNECT1D & newConnect,
                    OPTIONS   & options,
                    FIELD     & Field) const;

    bool matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
                   const Teuchos::RCP<CONNECT1D> & newConnect1d,
                   const Teuchos::RCP<OPTIONS>   & options,
                   const Teuchos::RCP<FIELD>     & field) const;
    
    bool matchGrid(MESH1D    & newGrid1d,
                   CONNECT1D & newConnect1d,
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
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a smaller communicator, all the pids not included in the new comm will be truncated */
    void reduceCommunicator(const bool                    & isActive,
                            const communicator            & OldCommDev,
                            const FIELD                   & OldField,
                                  communicator            & NewCommDev,
                            const Teuchos::RCP<MESH1D>    & NewGrid,
                            const Teuchos::RCP<CONNECT1D> & NewConnect,
                                  FIELD                   & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                             & isActive,
                            const Teuchos::RCP<const communicator> & OldCommDev,
                            const Teuchos::RCP<const FIELD>        & OldField,
                            const Teuchos::RCP<communicator>       & NewCommDev,
                            const Teuchos::RCP<MESH1D>             & NewGrid,
                            const Teuchos::RCP<CONNECT1D>          & NewConnect,
                                  Teuchos::RCP<FIELD>              & NewField);
    
    /*! Copy the map on a bigger communicator, all the pids not included are void */
    void expandCommunicator(const bool                     & isActive,
                            const communicator             & OldCommDev,
                            const FIELD                    & OldField,
                                  communicator             & NewCommDev,
                            const Teuchos::RCP<MESH1D>     & NewGrid,
                            const Teuchos::RCP<CONNECT1D>  & NewConnect,
                                  FIELD                    & NewField);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fsf1dGlobalManip()
{
  commDevLoaded = false;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fsf1dGlobalManip(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
fsf1dGlobalManip(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(const Teuchos::RCP<communicator> & CommDev)
{
  commDevLoaded = true;
  commDev = CommDev;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
setCommunicator(communicator & CommDev)
{
  commDevLoaded = true;
  commDev = Teuchos::rcpFromRef(CommDev);
}


//_________________________________________________________________________________________________
// OPERATOR FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const Teuchos::RCP<FIELD> & Field) const
{
  return(initilize(exprString,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
                FIELD          & Field) const
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH1D    MESH1D;
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
  Teuchos::RCP<MESH1D> grid1d = Field.getMesh1d();
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
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    assert(grid1d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid1d->getElementNodesL(i);
     
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string>      & exprString,
          const sVect<string>       & symbs,
                sVect<Real>           symbsVals,
          const Teuchos::RCP<FIELD> & Field) const
{
  return(initilize(exprString,symbs,symbsVals,*Field));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
initilize(const sArray<string> & exprString,
          const sVect<string>  & symbs,
                sVect<Real>      symbsVals,
                FIELD          & Field) const
{
  assert(commDevLoaded);
  assert(FETYPE::isNodal);
  assert(FETYPE::feBaseType == scalar);
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  assert(symbs.size() == symbsVals.size());
  
  //Typedefs---------------------------------------------------------
  typedef typename FIELD::MESH1D    MESH1D;
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
  Teuchos::RCP<MESH1D> grid1d = Field.getMesh1d();
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
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    assert(grid1d->getElements().getRowMapL(i).getLid() == i);
    dofCard.setLocalElId(i);
    
    if(Field.getDofMapper().isActive(dofCard))
    {
      //Element data extraction
      points = grid1d->getElementNodesL(i);
     
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(const Teuchos::RCP<MESH1D>    & newGrid,
           const Teuchos::RCP<CONNECT1D> & newConnect,
           const Teuchos::RCP<OPTIONS>   & options,
           const Teuchos::RCP<FIELD>     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field->getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  *Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
changeGrid(MESH1D    & newGrid,
           CONNECT1D & newConnect,
           OPTIONS   & options,
           FIELD     & Field) const
{
  assert(commDevLoaded);
  
  //Create the new field
  FIELD newField1d;
  newField1d.setCommunicator(commDev);
  newField1d.setGeometry(newGrid,newConnect);
  newField1d.setOptions(options);
  newField1d.startup();
  
  //New map - old vector
  pMap<PMAPTYPE> newDofMap  = newField1d.getDofVect().getMapRef();
  DOFVECT        oldDofVect = Field.getDofVect();
  
  //Cheking
  pVectGlobalManip<DOFTYPE,PMAPTYPE> dofVectManip(commDev);
  pMapGlobalManip<PMAPTYPE>          dofMapManip(commDev);
  
  assert(dofVectManip.sizeG(oldDofVect) == dofMapManip.sizeG(newDofMap));
  
  //New dofVector
  dofVectManip.changeMap(oldDofVect,newDofMap);
  
  //Loading the new dofVector
  newField1d.setDofVect(oldDofVect);
  Field = newField1d;
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(const Teuchos::RCP<MESH1D>    & newGrid1d,
          const Teuchos::RCP<CONNECT1D> & newConnect1d,
          const Teuchos::RCP<OPTIONS>   & options,
          const Teuchos::RCP<FIELD>     & field) const
{
  changeGrid(newGrid1d,
             newConnect1d,
             options,
             field);
  
  return(true);
}
    
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
bool
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
matchGrid(MESH1D    & newGrid1d,
          CONNECT1D & newConnect1d,
          OPTIONS   & options,
          FIELD     & field) const
{
  changeGrid(newGrid1d,
             newConnect1d,
             options,
             field);
  
  return(true);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const Teuchos::RCP<FIELD>    & Field,
                        Teuchos::RCP<SENDRECV> & mapSend,
                        Teuchos::RCP<SENDRECV> & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field.getDofVect().getMapRcp(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
createSendRecvMap(const FIELD    & Field,
                        SENDRECV & mapSend,
                        SENDRECV & mapRecv) const
{
  assert(commDevLoaded);
  
  pMapGlobalManip<PMAPTYPE> mapManipulator(commDev);
  mapManipulator.createSendRecvMap(Field.getDofVect().getMapRef(), mapSend, mapRecv);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(const Teuchos::RCP<FIELD>    & Field,
           const Teuchos::RCP<SENDRECV> & mapSend,
           const Teuchos::RCP<SENDRECV> & mapRecv) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field->getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,*mapSend,*mapRecv);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateData(FIELD    & Field,
           SENDRECV & mapSend,
           SENDRECV & mapRecv) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
  //Update the data
  pVect<DOFTYPE,pMapItemShare> dofsVect = Field.getDofVect();
  
  pVectGlobalManip<DOFTYPE,pMapItemShare> manipulator(commDev);
  manipulator.updateData(dofsVect,mapSend,mapRecv);
  
  //Reload the dofs
  Field.setDofVect(dofsVect);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
Real
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
evalMax(const Teuchos::RCP<FIELD> & Field,
        const sVect<point3d>      & evalNodes,
        const set<UInt>           & activeIds)
{
  return(evalMax(*Field,evalNodes,activeIds));
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
Real
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
  Teuchos::RCP<MESH1D> grid1d = Field.getMesh();
  
  //Find max
  for(UInt i=1; i <= grid1d->getNumElements(); ++i)
  {
    geoId = grid1d->getElementL(i).getGeoId();
    
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
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
   
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>     
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
   
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>     
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataI(FIELD          & Field,
            SENDRECV       & mapSend,
            SENDRECV       & mapRecv,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
updateDataO(FIELD          & Field,
            SENDRECV       & mapSend,
            SENDRECV       & mapRecv,
            sVect<DOFVECT> & bufSegments,
            sVect<PVPS>    & sendPvps,
            sVect<PVPS>    & recvPvps,
            const UInt     & channel) const
{
  //Checks
  assert(commDevLoaded);
  pMapGlobalManip<pMapItemSendRecv> cheker(commDev);
  cheker.check(mapSend,mapRecv);
  
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
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
      
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
      
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
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
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
reduceCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                             & isActive,
                   const Teuchos::RCP<const communicator> & OldCommDev,
                   const Teuchos::RCP<const FIELD>        & OldField,
                   const Teuchos::RCP<communicator>       & NewCommDev,
                   const Teuchos::RCP<MESH1D>             & NewGrid,
                   const Teuchos::RCP<CONNECT1D>          & NewConnect,
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

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
void
fsf1dGlobalManip<pMapItemShare,FETYPE,DOFTYPE,ORDER,MODE>::
expandCommunicator(const bool                    & isActive,
                   const communicator            & OldCommDev,
                   const FIELD                   & OldField,
                         communicator            & NewCommDev,
                   const Teuchos::RCP<MESH1D>    & NewGrid,
                   const Teuchos::RCP<CONNECT1D> & NewConnect,
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

/*! Interface class for the manipulation of Static Field 2d */
template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER = dms1d_vectMajor, dms1d_mode MODE = dms1d_allMode>
class feStaticField1dGlobalManip : public fsf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>
{
     /*! @name Typedefs */ //@{
  public:
    typedef typename FETYPE::FETYPE_PMAPTYPE PMAPTYPE;
    typedef feStaticField1d<FETYPE,DOFTYPE,ORDER,MODE> FIELD;
    typedef typename FETYPE::GEOSHAPE  GEOSHAPE;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    feStaticField1dGlobalManip();
    feStaticField1dGlobalManip(const Teuchos::RCP<communicator> & CommDev);
    feStaticField1dGlobalManip(communicator & CommDev);
    //@}
};


template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
feStaticField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField1dGlobalManip() : fsf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>()
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
feStaticField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField1dGlobalManip(const Teuchos::RCP<communicator> & CommDev) : fsf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

template<typename FETYPE, typename DOFTYPE, dms1d_order ORDER, dms1d_mode MODE>
feStaticField1dGlobalManip<FETYPE,DOFTYPE,ORDER,MODE>::
feStaticField1dGlobalManip(communicator & CommDev) : fsf1dGlobalManip<typename FETYPE::FETYPE_PMAPTYPE,FETYPE,DOFTYPE,ORDER,MODE>(CommDev)
{
  assert(staticAssert<GEOSHAPE::nDim == 1>::returnValue);
}

#endif
