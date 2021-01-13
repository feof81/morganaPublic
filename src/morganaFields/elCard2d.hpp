/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef ELCARD2D_HPP
#define ELCARD2D_HPP

#include "typesInterface.hpp"

#include "geoShapes.h"
#include "morganaGeometry.hpp"
#include "mesh2d.hpp"


/*! Contains some topological information of a 2d element */
template<typename GEOSHAPE, typename PMAPTYPE>
class elCard2d
{
  /*! @name Typedefs */ //@{
  public:   
    typedef sVect<point3d> NODES;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    NODES nodes;
    //@}
    
    /*! @name Constructors and operators */ //@{
  public:
    elCard2d();
    elCard2d(const elCard2d & inCard);
    elCard2d operator=(const elCard2d & inCard);
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
    friend ostream & operator<<(ostream & f, const elCard2d<G,T> & V);
    //@}
};


//_________________________________________________________________________________________________
// CONSTRUCTORS AND OPERATORS
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
elCard2d<GEOSHAPE,PMAPTYPE>::
elCard2d()
{
  nodes.resize(GEOSHAPE::numVertices);
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard2d<GEOSHAPE,PMAPTYPE>::
elCard2d(const elCard2d & inCard)
{ 
  nodes = inCard.nodes;
}

template<typename GEOSHAPE, typename PMAPTYPE>
elCard2d<GEOSHAPE,PMAPTYPE>
elCard2d<GEOSHAPE,PMAPTYPE>::
operator=(const elCard2d & inCard)
{
  nodes = inCard.nodes;  
  return(*this);
}



//_________________________________________________________________________________________________
// GET FUNCTIONS - REFERENCE
//-------------------------------------------------------------------------------------------------
template<typename GEOSHAPE, typename PMAPTYPE>
typename elCard2d<GEOSHAPE,PMAPTYPE>::NODES &
elCard2d<GEOSHAPE,PMAPTYPE>::
getNodes()
{
  return(nodes);
}


template<typename GEOSHAPE, typename PMAPTYPE>
point3d &
elCard2d<GEOSHAPE,PMAPTYPE>::
getNode(const UInt & i)
{
  assert(i <= nodes.size());
  return(nodes(i));
}


template<typename GEOSHAPE, typename PMAPTYPE>
const point3d &
elCard2d<GEOSHAPE,PMAPTYPE>::
getNode(const UInt & i) const
{
  assert(i <= nodes.size());
  return(nodes(i));
}

template<typename GEOSHAPE, typename PMAPTYPE>
const typename  elCard2d<GEOSHAPE,PMAPTYPE>::NODES &
elCard2d<GEOSHAPE,PMAPTYPE>::
getNodes() const
{
  return(nodes);
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
elCard2d<GEOSHAPE,PMAPTYPE>::
setNodes(const NODES & Nodes)
{
  nodes = Nodes;
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
elCard2d<GEOSHAPE,PMAPTYPE>::
setNode(const UInt & i, const point3d & P)
{
  assert(i <= nodes.size());
  nodes(i) = P;
}


//_______________________________________________________________________________________________________
// PRINTOUT
//-------------------------------------------------------------------------------------------------------
template<typename G, typename T>
ostream & operator<<(ostream & f, const elCard2d<G,T> & V)
{
  f << "Nodes--------------------------------" << endl;
  f << V.nodes << endl << endl;  
  
  return(f);
}



#endif
