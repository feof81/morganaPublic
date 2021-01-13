/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef ANALYTICFIELD_HPP
#define ANALYTICFIELD_HPP

#include "typesInterface.hpp"
#include "traitsBasic.h"
#include "exprtk.hpp"


/*! Analytic definition of a function */
template<typename DOFTYPE>
class analyticField
{
    /*! @name Typedefs */ //@{
  public:
    typedef exprtk::symbol_table<Real> SYMBOLTABLE;
    typedef exprtk::expression<Real>   EXPRESSION;
    typedef exprtk::parser<Real>       PARSER;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    bool startupOk;
    Real x,y,z;
    sArray<string>                    exprString;
    Teuchos::RCP<SYMBOLTABLE>         symbol_table;
    sArray<Teuchos::RCP<EXPRESSION> > expressions;
    sArray<Teuchos::RCP<PARSER> >     parsers;
    //@}
    
    /*! @name Functions */ //@{
  public:
    analyticField();
    analyticField(const sArray<string> & ExprString);
    analyticField(const analyticField & A);
    void set(const sArray<string> & ExprString);
    void startup();
    DOFTYPE eval(const point3d & Pg);
    //@}
};


//_________________________________________________________________________________________________
// FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename DOFTYPE>
analyticField<DOFTYPE>::
analyticField()
{
  startupOk = false;
}

template<typename DOFTYPE>
analyticField<DOFTYPE>::
analyticField(const sArray<string> & ExprString)
{
  assert(ExprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(ExprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  
  exprString = ExprString;
  startup();
}

template<typename DOFTYPE>
analyticField<DOFTYPE>::
analyticField(const analyticField & A)
{
  exprString = A.exprString;
  startup();
}

template<typename DOFTYPE>
void
analyticField<DOFTYPE>::
set(const sArray<string> & ExprString)
{
  assert(ExprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(ExprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  
  exprString = ExprString;
  startup();
}

template<typename DOFTYPE>
void
analyticField<DOFTYPE>::
startup()
{
  //Assert
  assert(exprString.nrows() == traitsBasic<DOFTYPE>::numI);
  assert(exprString.ncols() == traitsBasic<DOFTYPE>::numJ);
  
  //Logics
  startupOk = true;
  
  //Symbol table
  symbol_table = Teuchos::rcp(new SYMBOLTABLE);
  symbol_table->add_variable("x",x);
  symbol_table->add_variable("y",y);
  symbol_table->add_variable("z",z);
  symbol_table->add_constants();
  
  //Startup
  expressions.reshape(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  parsers.reshape(traitsBasic<DOFTYPE>::numI, traitsBasic<DOFTYPE>::numJ);
  
  bool logic;
  
  for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)
  {
    for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
    {
      expressions(I,J) = Teuchos::rcp(new EXPRESSION);
      expressions(I,J)->register_symbol_table(*symbol_table);
      
      parsers(I,J) = Teuchos::rcp(new PARSER);
      logic = parsers(I,J)->compile(exprString(I,J), *expressions(I,J));
      
      assert(logic);
    }
  }
}

template<typename DOFTYPE>
DOFTYPE
analyticField<DOFTYPE>::
eval(const point3d & Pg)
{
  //Assert
  assert(startupOk);
  
  //Evaluation
  DOFTYPE V;
  Real temp;
  
  for(UInt I=1; I <= traitsBasic<DOFTYPE>::numI; ++I)  
  {
    for(UInt J=1; J <= traitsBasic<DOFTYPE>::numJ; ++J)
    {
      x = Pg.getX();
      y = Pg.getY();
      z = Pg.getZ();
      
      temp = expressions(I,J)->value();
      traitsBasic<DOFTYPE>::setIJ(I,J,temp,V);
    }
  }
  
  return(V);
}


#endif
