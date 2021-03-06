/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef ELCARD1D_HPP
#define ELCARD1D_HPP

#include "typesInterface.hpp"

#include "geoShapes.h"
#include "morganaGeometry.hpp"
#include "mesh2d.hpp"


/*! Contains some topological information of a 1d element */
template<typename GEOSHAPE, typename PMAPTYPE>
class elCard1d
{
  /*! @name Typedefs */ //@{
  public:   
    typedef sVect<point3d>  NODES;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    NODES nodes;
    //@}
    
    /*! @name Constructors and operators */ //@{
  public:
    elCard1d();
    elCard1d(const elCard1d & inCard);
    elCard1d operator=(const elCard1d & inCard);
    //@}
    
    /*! @name Get - Set functions */ //@{
  public:
    NODES         & getNodes();
    point3d       & getNode(const UInt & i);
    const point3d & getNode(const UInt & i) const;
    const NODES   & getNodes() const;
    void setNodes(const NODES & Nodes);
    void setNode(const UInt & i, const point3d & P);
    //@}
    
    /*! @name Printout */ //@{
  public:
    template<typename G, typename T>
    friend ostream & operator<<(ostream & f, const elCard1d<G,T> & V);
    //@}
};



//_________________________________________________________________________________________________
// CONSTRUCTORS AND OPERATORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
elCard1d<GEOSHAPE,PMAPTYPE>::
elCard1d()
{
  nodes.resize(GEOSHAPE::numVertices);
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard1d<GEOSHAPE,PMAPTYPE>::
elCard1d(const elCard1d & inCard)
{ 
  nodes = inCard.nodes;
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard1d<GEOSHAPE,PMAPTYPE>
elCard1d<GEOSHAPE,PMAPTYPE>::
operator=(const elCard1d & inCard)
{
  nodes = inCard.nodes;
  
  return(*this);
}



//_________________________________________________________________________________________________
// GET FUNCTIONS - REFERENCE
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
typename elCard1d<GEOSHAPE,PMAPTYPE>::NODES &
elCard1d<GEOSHAPE,PMAPTYPE>::
getNodes()
{
  return(nodes);
}


//_________________________________________________________________________________________________
// GET FUNCTIONS - REFERENCE/ITEM
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
point3d &
elCard1d<GEOSHAPE,PMAPTYPE>::
getNode(const UInt & i)
{
  assert(i <= nodes.size());
  return(nodes(i));
}



//_________________________________________________________________________________________________
// GET FUNCTIONS - REFERENCE/ITEM
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
const point3d &
elCard1d<GEOSHAPE,PMAPTYPE>::
getNode(const UInt & i) const
{
  assert(i <= nodes.size());
  return(nodes(i));
}



//_________________________________________________________________________________________________
// GET FUNCTIONS - CONST REFERENCE
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
const typename  elCard1d<GEOSHAPE,PMAPTYPE>::NODES &
elCard1d<GEOSHAPE,PMAPTYPE>::
getNodes() const
{
  return(nodes);
}



//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
void
elCard1d<GEOSHAPE,PMAPTYPE>::
setNodes(const NODES & Nodes)
{
  nodes = Nodes;
}



//_______________________________________________________________________________________________________
// SET FUNCTIONS PER ITEM
//-------------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
void
elCard1d<GEOSHAPE,PMAPTYPE>::
setNode(const UInt & i, const point3d & P)
{
  assert(i <= nodes.size());
  nodes(i) = P;
}



//_______________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------------
template<typename G, typename T>
ostream & operator<<(ostream & f, const elCard1d<G,T> & V)
{
  f << "Nodes--------------------------------" << endl;
  f << V.nodes << endl << endl;
   
  return(f);
}


#endif
