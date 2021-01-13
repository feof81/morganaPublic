/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEP1BP13D_HPP
#define FEP1BP13D_HPP

#include "morganaTypes.hpp"

#include "geoMapInterface.hpp"

#include "morganaFields.hpp"
#include "elCard3d.hpp"
#include "feStaticDofCard3d.h"
#include "feStaticEvalIterators.hpp"


//_________________________________________________________________________________________________
// FE CARD
//-------------------------------------------------------------------------------------------------
class feP1BP13d_card
{
  public:
    feP1BP13d_card() { };
};


//Forward declaration______________________________________________________________________________
enum feP1B13d_piece {feP1B13d_full, feP1B13d_regular, feP1B13d_bubble};

template<typename PMAPTYPE, feP1B13d_piece FEPIECE = feP1B13d_full> class feP1BP13d;


//Full description_________________________________________________________________________________

/*! Linear tetrahedral finite element with a linear bubble enrichment */
template<typename PMAPTYPE>
class feP1BP13d< PMAPTYPE,feP1B13d_full>
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feP1BP13d_card              FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_P1B1_3d;
    static const FEBaseLabel feBaseLabel = BL_P1B1_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 8;
    static const UInt numBasis = 5;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 1;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 0;
    static const UInt dofPerVolume = 1;
    static const bool isNodal = false;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool elCardLoaded;
    bool feCardLoaded;
    FECARD FeCard;
    ELCARD ElCard;
    //@}
    
    /*! @name Constructors and setting functions */ //@{
  public:
    feP1BP13d();
    feP1BP13d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    bool baseIsActive(const UInt & i, const point3d & Y) const;
    //@}
    
    /*! @name Eval functions */ //@{
  public:
    static DOFCARD getDofCard(const UInt & i);
    static point3d getRefNode(const UInt & i);
    void localEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    void globalEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    //@}
};



template<typename PMAPTYPE>
feP1BP13d< PMAPTYPE,feP1B13d_full>::
feP1BP13d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feP1BP13d< PMAPTYPE,feP1B13d_full>::
feP1BP13d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feP1BP13d< PMAPTYPE,feP1B13d_full>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
bool
feP1BP13d< PMAPTYPE,feP1B13d_full>::
baseIsActive(const UInt & i, const point3d & Y) const
{
   Real dist[4];
   
   static const point3d N(1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0));
   static const point3d X(1.0,0.0,0.0);
    
   dist[0] = fabs( point3d::dot(Y-X, N) );
   dist[1] = Y.getX();
   dist[2] = Y.getY();
   dist[3] = Y.getZ();
    
   return( min( min(dist[0],dist[1]), min(dist[2],dist[3])) == dist[(i-1)] );
}

template<typename PMAPTYPE>
typename feP1BP13d< PMAPTYPE,feP1B13d_full>::DOFCARD
feP1BP13d< PMAPTYPE,feP1B13d_full>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 5);
  
  static DOFCARD list[5] =
  {
    DOFCARD(VERTEX,1,1),
    DOFCARD(VERTEX,1,2),
    DOFCARD(VERTEX,1,3),
    DOFCARD(VERTEX,1,4),
    DOFCARD(VOLUME,1,1) };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feP1BP13d< PMAPTYPE,feP1B13d_full>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 5);
  
  static point3d list[5] =
  {
    point3d(0.0,  0.0,  0.0),
    point3d(1.0,  0.0,  0.0),
    point3d(0.0,  1.0,  0.0),
    point3d(0.0,  0.0,  1.0),
    point3d(0.25, 0.25, 0.25)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feP1BP13d< PMAPTYPE,feP1B13d_full>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numBasis);
  
  //Polynomial evaluation
  sVect<Real> polyVal(numPoly);
  feStaticRealEval<feLabel,numPoly>::eval(polyVal,Y);
  
  //Active terms
  val(1) = polyVal(1);
  val(2) = polyVal(2);
  val(3) = polyVal(3);
  val(4) = polyVal(4);
  val(5) = polyVal(5) * baseIsActive(1,Y) +
           polyVal(6) * baseIsActive(2,Y) +
           polyVal(7) * baseIsActive(3,Y) +
           polyVal(8) * baseIsActive(4,Y);
	   	   
  assert( (baseIsActive(1,Y) +
           baseIsActive(2,Y) +
	   baseIsActive(3,Y) +
	   baseIsActive(4,Y) ) == 1 );
}

template<typename PMAPTYPE>
void
feP1BP13d< PMAPTYPE,feP1B13d_full>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Polynomial evaluation
  sVect<Real> polyGradX(numPoly);
  sVect<Real> polyGradY(numPoly);
  sVect<Real> polyGradZ(numPoly);
  
  feStaticRealDx<feLabel,numPoly>::eval(polyGradX,Y);
  feStaticRealDy<feLabel,numPoly>::eval(polyGradY,Y);
  feStaticRealDz<feLabel,numPoly>::eval(polyGradZ,Y);
  
  //Active terms X
  gradX(1) = polyGradX(1);
  gradX(2) = polyGradX(2);
  gradX(3) = polyGradX(3);
  gradX(4) = polyGradX(4);
  gradX(5) = polyGradX(5) * baseIsActive(1,Y) +
             polyGradX(6) * baseIsActive(2,Y) +
             polyGradX(7) * baseIsActive(3,Y) +
             polyGradX(8) * baseIsActive(4,Y);
	     
  //Active terms Y	     
  gradY(1) = polyGradY(1);
  gradY(2) = polyGradY(2);
  gradY(3) = polyGradY(3);
  gradY(4) = polyGradY(4);
  gradY(5) = polyGradY(5) * baseIsActive(1,Y) +
             polyGradY(6) * baseIsActive(2,Y) +
             polyGradY(7) * baseIsActive(3,Y) +
             polyGradY(8) * baseIsActive(4,Y);
	          
  //Active terms Y	     
  gradZ(1) = polyGradZ(1);
  gradZ(2) = polyGradZ(2);
  gradZ(3) = polyGradZ(3);
  gradZ(4) = polyGradZ(4);
  gradZ(5) = polyGradZ(5) * baseIsActive(1,Y) +
             polyGradZ(6) * baseIsActive(2,Y) +
             polyGradZ(7) * baseIsActive(3,Y) +
             polyGradZ(8) * baseIsActive(4,Y);
	     
	     
  assert( (baseIsActive(1,Y) +
           baseIsActive(2,Y) +
	   baseIsActive(3,Y) +
	   baseIsActive(4,Y) ) == 1 );
}

template<typename PMAPTYPE>
void
feP1BP13d< PMAPTYPE,feP1B13d_full>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numBasis);
  
  //Polynomial evaluation
  localEval(Y,val) ;
}

template<typename PMAPTYPE>
void
feP1BP13d< PMAPTYPE,feP1B13d_full>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
  T.computeInverse();
  
  //Local grad
  sVect<BASETYPE> locX(numBasis), locY(numBasis), locZ(numBasis);
  localGrad(Y,locX,locY,locZ);
  
  //Compute gradient
  for(UInt i=1; i <= numBasis; ++i)
  {
    gradX(i) = locX(i) * T(1,1) + locY(i) * T(2,1) + locZ(i) * T(3,1);
    gradY(i) = locX(i) * T(1,2) + locY(i) * T(2,2) + locZ(i) * T(3,2);
    gradZ(i) = locX(i) * T(1,3) + locY(i) * T(2,3) + locZ(i) * T(3,3);
  }
}




//Description without bubble_______________________________________________________________________

/*! Linear tetrahedral finite element with a linear bubble enrichment -> without the bubble */
template<typename PMAPTYPE>
class feP1BP13d<PMAPTYPE, feP1B13d_regular>
{
  /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feP1BP13d_card              FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass    feClass    = feStatic;
    static const FEGridType feGridType = primal;
    static const FEBaseType feBaseType = scalar;
    static const FELabel    feLabel    = FE_P1B1_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 4;
    static const UInt numBasis = 4;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 1;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 0;
    static const UInt dofPerVolume = 0;
    static const bool isNodal = true;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool elCardLoaded;
    bool feCardLoaded;
    FECARD FeCard;
    ELCARD ElCard;
    //@}
    
    /*! @name Constructors and setting functions */ //@{
  public:
    feP1BP13d();
    feP1BP13d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    //@}
    
    /*! @name Eval functions */ //@{
  public:
    static DOFCARD getDofCard(const UInt & i);
    static point3d getRefNode(const UInt & i);
    void localEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    void globalEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    //@}
};

template<typename PMAPTYPE>
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
feP1BP13d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
feP1BP13d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}
    
template<typename PMAPTYPE>
typename feP1BP13d<PMAPTYPE, feP1B13d_regular>::DOFCARD
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  static DOFCARD list[4] =
  {
    DOFCARD(VERTEX,1,1),
    DOFCARD(VERTEX,1,2),
    DOFCARD(VERTEX,1,3),
    DOFCARD(VERTEX,1,4) };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  static point3d list[4] =
  {
    point3d(0.0,  0.0,  0.0),
    point3d(1.0,  0.0,  0.0),
    point3d(0.0,  1.0,  0.0),
    point3d(0.0,  0.0,  1.0)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation
  feStaticRealEval<feLabel,numPoly>::eval(val,Y);
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  //Polynomial evaluation
  feStaticRealDx<feLabel,numPoly>::eval(gradX,Y);
  feStaticRealDy<feLabel,numPoly>::eval(gradY,Y);
  feStaticRealDz<feLabel,numPoly>::eval(gradZ,Y);
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numBasis);
  
  //Polynomial evaluation
  feStaticRealEval<feLabel,numBasis>::eval(val,Y);
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_regular>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
  T.computeInverse();
  
  //Local grad
  sVect<BASETYPE> locX(numBasis), locY(numBasis), locZ(numBasis);
  localGrad(Y,locX,locY,locZ);
  
  //Compute gradient
  for(UInt i=1; i <= numBasis; ++i)
  {
    gradX(i) = locX(i) * T(1,1) + locY(i) * T(2,1) + locZ(i) * T(3,1);
    gradY(i) = locX(i) * T(1,2) + locY(i) * T(2,2) + locZ(i) * T(3,2);
    gradZ(i) = locX(i) * T(1,3) + locY(i) * T(2,3) + locZ(i) * T(3,3);
  }
}


//Description without bubble_______________________________________________________________________

/*! Linear tetrahedral finite element with a linear bubble enrichment -> bubble only */
template<typename PMAPTYPE>
class feP1BP13d<PMAPTYPE, feP1B13d_bubble>
{
  /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feP1BP13d_card              FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass    feClass    = feStatic;
    static const FEGridType feGridType = primal;
    static const FEBaseType feBaseType = scalar;
    static const FELabel    feLabel    = FE_P1B1_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 4;
    static const UInt numBasis = 1;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 0;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 0;
    static const UInt dofPerVolume = 1;
    static const bool isNodal = true;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool elCardLoaded;
    bool feCardLoaded;
    FECARD FeCard;
    ELCARD ElCard;
    //@}
    
    /*! @name Constructors and setting functions */ //@{
  public:
    feP1BP13d();
    feP1BP13d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    bool baseIsActive(const UInt & i, const point3d & Y) const;
    //@}
    
    /*! @name Eval functions */ //@{
  public:
    static DOFCARD getDofCard(const UInt & i);
    static point3d getRefNode(const UInt & i);
    void localEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    void globalEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    //@}
};



template<typename PMAPTYPE>
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
feP1BP13d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
feP1BP13d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
bool
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
baseIsActive(const UInt & i, const point3d & Y) const
{
   Real dist[4];
   
   point3d N(1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0));
   point3d X(1.0,0.0,0.0);
    
   dist[0] = fabs( point3d::dot(Y-X, N) );
   dist[1] = Y.getX();
   dist[2] = Y.getY();
   dist[3] = Y.getZ();
    
   return( min( min(dist[0],dist[1]), min(dist[2],dist[3])) == dist[(i-1)] );
}

template<typename PMAPTYPE>
typename feP1BP13d<PMAPTYPE, feP1B13d_bubble>::DOFCARD
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
getDofCard(const UInt & i)
{
  assert(i == 1);
  
  static DOFCARD list[1] =
  {
    DOFCARD(VOLUME,1,1)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 1);
  
  static point3d list[1] =
  {
    point3d(0.25, 0.25, 0.25)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numBasis);
  
  //Polynomial evaluation
  sVect<Real> polyVal(numPoly);
  feStaticRealEval<feLabel,numPoly>::eval(polyVal,Y);
  
  //Active terms
  val(1) = polyVal(1) * baseIsActive(1,Y) +
           polyVal(2) * baseIsActive(2,Y) +
           polyVal(3) * baseIsActive(3,Y) +
           polyVal(4) * baseIsActive(4,Y);
	   	   
  assert( (baseIsActive(1,Y) +
           baseIsActive(2,Y) +
	   baseIsActive(3,Y) +
	   baseIsActive(4,Y) ) == 1 );
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Polynomial evaluation
  sVect<Real> polyGradX(numPoly);
  sVect<Real> polyGradY(numPoly);
  sVect<Real> polyGradZ(numPoly);
  
  feStaticRealDx<feLabel,numPoly>::eval(polyGradX,Y);
  feStaticRealDy<feLabel,numPoly>::eval(polyGradY,Y);
  feStaticRealDz<feLabel,numPoly>::eval(polyGradZ,Y);
  
  //Active terms X
  gradX(1) = polyGradX(1) * baseIsActive(1,Y) +
             polyGradX(2) * baseIsActive(2,Y) +
             polyGradX(3) * baseIsActive(3,Y) +
             polyGradX(4) * baseIsActive(4,Y);
	     
  //Active terms Y	     
  gradY(1) = polyGradY(1) * baseIsActive(1,Y) +
             polyGradY(2) * baseIsActive(2,Y) +
             polyGradY(3) * baseIsActive(3,Y) +
             polyGradY(4) * baseIsActive(4,Y);
	          
  //Active terms Y	     
  gradZ(1) = polyGradZ(1) * baseIsActive(1,Y) +
             polyGradZ(2) * baseIsActive(2,Y) +
             polyGradZ(3) * baseIsActive(3,Y) +
             polyGradZ(4) * baseIsActive(4,Y);
	     
	     
  assert( (baseIsActive(1,Y) +
           baseIsActive(2,Y) +
	   baseIsActive(3,Y) +
	   baseIsActive(4,Y) ) == 1 );
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numBasis);
  
  //Polynomial evaluation
  localEval(Y,val) ;
}

template<typename PMAPTYPE>
void
feP1BP13d<PMAPTYPE, feP1B13d_bubble>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
  T.computeInverse();
  
  //Local grad
  sVect<BASETYPE> locX(numBasis), locY(numBasis), locZ(numBasis);
  localGrad(Y,locX,locY,locZ);
  
  //Compute gradient
  for(UInt i=1; i <= numBasis; ++i)
  {
    gradX(i) = locX(i) * T(1,1) + locY(i) * T(2,1) + locZ(i) * T(3,1);
    gradY(i) = locX(i) * T(1,2) + locY(i) * T(2,2) + locZ(i) * T(3,2);
    gradZ(i) = locX(i) * T(1,3) + locY(i) * T(2,3) + locZ(i) * T(3,3);
  }
}
  
#endif
