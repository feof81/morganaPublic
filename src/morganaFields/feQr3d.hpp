/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEQR3D_HPP
#define FEQR3D_HPP

#include "morganaTypes.hpp"

#include "geoMapInterface.hpp"

#include "morganaFields.hpp"
#include "elCard3d.hpp"
#include "feStaticDofCard3d.h"
#include "feStaticEvalIterators.hpp"



//_________________________________________________________________________________________________
// FE CARD
//-------------------------------------------------------------------------------------------------

/*! Dummy card */
class feQr3d_card
{
  public:
    feQr3d_card() { };
};



//Forward declaration______________________________________________________________________________
template<Int R, typename PMAPTYPE> class feQr3d;



//_________________________________________________________________________________________________
//   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0   Q0  
//-------------------------------------------------------------------------------------------------

/*! Hexahedral nodal finite element - Q0 */
template<typename PMAPTYPE>
class feQr3d<0, PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearHexa                  GEOSHAPE;
    typedef feQr3d_card                 FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
     
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_Q0_3d;
    static const FEBaseLabel feBaseLabel = BL_Qr_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 1;
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
    feQr3d();
    feQr3d(const FECARD & FECard, const ELCARD & ELCard);
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
feQr3d<0, PMAPTYPE>::
feQr3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feQr3d<0, PMAPTYPE>::
feQr3d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}
    
template<typename PMAPTYPE>
void
feQr3d<0, PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}
    
template<typename PMAPTYPE>
typename feQr3d<0, PMAPTYPE>::DOFCARD
feQr3d<0, PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i == 1);  
  return( DOFCARD(VOLUME,1,1) );
}

template<typename PMAPTYPE>
point3d
feQr3d<0, PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i == 1);  
  return(point3d(0.5, 0.5, 0.5));
}

template<typename PMAPTYPE>
void
feQr3d<0, PMAPTYPE>::
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
feQr3d<0, PMAPTYPE>::
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
feQr3d<0, PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation
  feStaticRealEval<feLabel,numPoly>::eval(val,Y);
}

template<typename PMAPTYPE>
void
feQr3d<0, PMAPTYPE>::
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



//_________________________________________________________________________________________________
//   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1   Q1  
//-------------------------------------------------------------------------------------------------

/*! Hexahedral nodal finite element - Q1 */
template<typename PMAPTYPE>
class feQr3d<1, PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearHexa                  GEOSHAPE;
    typedef feQr3d_card                 FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
     
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_Q1_3d;
    static const FEBaseLabel feBaseLabel = BL_Qr_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 8;
    static const UInt numBasis = 8;
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
    feQr3d();
    feQr3d(const FECARD & FECard, const ELCARD & ELCard);
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
feQr3d<1, PMAPTYPE>::
feQr3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feQr3d<1, PMAPTYPE>::
feQr3d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}
    
template<typename PMAPTYPE>
void
feQr3d<1, PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feQr3d<1, PMAPTYPE>::DOFCARD
feQr3d<1, PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 8);
  
  static DOFCARD list[8] =
  {
    DOFCARD(VERTEX,1,1),
    DOFCARD(VERTEX,1,2),
    DOFCARD(VERTEX,1,3),
    DOFCARD(VERTEX,1,4),
    DOFCARD(VERTEX,1,5),
    DOFCARD(VERTEX,1,6),
    DOFCARD(VERTEX,1,7),
    DOFCARD(VERTEX,1,8) };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feQr3d<1, PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 8);
  
  static point3d list[8] =
  {
    point3d(0.0,  0.0,  0.0),
    point3d(1.0,  0.0,  0.0),
    point3d(1.0,  1.0,  0.0),
    point3d(0.0,  1.0,  0.0),
    
    point3d(0.0,  0.0,  1.0),
    point3d(1.0,  0.0,  1.0),
    point3d(1.0,  1.0,  1.0),
    point3d(0.0,  1.0,  1.0)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feQr3d<1, PMAPTYPE>::
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
feQr3d<1, PMAPTYPE>::
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
feQr3d<1, PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation
  feStaticRealEval<feLabel,numPoly>::eval(val,Y);
}

template<typename PMAPTYPE>
void
feQr3d<1, PMAPTYPE>::
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



//_________________________________________________________________________________________________
//   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2   Q2  
//-------------------------------------------------------------------------------------------------

/*! Hexahedral nodal finite element - Q2 */
template<typename PMAPTYPE>
class feQr3d<2, PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearHexa                  GEOSHAPE;
    typedef feQr3d_card                 FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
     
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_Q2_3d;
    static const FEBaseLabel feBaseLabel = BL_Qr_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 27;
    static const UInt numBasis = 27;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 1;
    static const UInt dofPerEdge   = 1;
    static const UInt dofPerFace   = 1;
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
    feQr3d();
    feQr3d(const FECARD & FECard, const ELCARD & ELCard);
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
feQr3d<2, PMAPTYPE>::
feQr3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feQr3d<2, PMAPTYPE>::
feQr3d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}
    
template<typename PMAPTYPE>
void
feQr3d<2, PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feQr3d<2, PMAPTYPE>::DOFCARD
feQr3d<2, PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 27);
  
  static DOFCARD list[27] =
  {
    DOFCARD(VERTEX,1,1),
    DOFCARD(VERTEX,1,2),
    DOFCARD(VERTEX,1,3),
    DOFCARD(VERTEX,1,4),
    DOFCARD(EDGE,1,1),
    DOFCARD(EDGE,1,2),
    DOFCARD(EDGE,1,3),
    DOFCARD(EDGE,1,4),
    DOFCARD(FACE,1,1),
    
    DOFCARD(EDGE,1,5),
    DOFCARD(EDGE,1,6),
    DOFCARD(EDGE,1,7),
    DOFCARD(EDGE,1,8),
    DOFCARD(FACE,1,2),
    DOFCARD(FACE,1,3),
    DOFCARD(FACE,1,4),
    DOFCARD(FACE,1,5),
    DOFCARD(VOLUME,1,1),
    
    DOFCARD(VERTEX,1,5),
    DOFCARD(VERTEX,1,6),
    DOFCARD(VERTEX,1,7),
    DOFCARD(VERTEX,1,8),
    DOFCARD(EDGE,1,9),
    DOFCARD(EDGE,1,10),
    DOFCARD(EDGE,1,11),
    DOFCARD(EDGE,1,12),
    DOFCARD(FACE,1,6) };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feQr3d<2, PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 27);
  
  static point3d list[27] =
  {
    point3d(0.0,  0.0,  0.0),
    point3d(1.0,  0.0,  0.0),
    point3d(1.0,  1.0,  0.0),
    point3d(0.0,  1.0,  0.0),
    point3d(0.5,  0.0,  0.0),
    point3d(1.0,  0.5,  0.0),
    point3d(0.5,  1.0,  0.0),
    point3d(0.0,  0.5,  0.0),
    point3d(0.5,  0.5,  0.0),
    
    point3d(0.0,  0.0,  0.5),
    point3d(1.0,  0.0,  0.5),
    point3d(1.0,  1.0,  0.5),
    point3d(0.0,  1.0,  0.5),
    point3d(0.5,  0.0,  0.5),
    point3d(1.0,  0.5,  0.5),
    point3d(0.5,  1.0,  0.5),
    point3d(0.0,  0.5,  0.5),
    point3d(0.5,  0.5,  0.5),
    
    point3d(0.0,  0.0,  1.0),
    point3d(1.0,  0.0,  1.0),
    point3d(1.0,  1.0,  1.0),
    point3d(0.0,  1.0,  1.0),
    point3d(0.5,  0.0,  1.0),
    point3d(1.0,  0.5,  1.0),
    point3d(0.5,  1.0,  1.0),
    point3d(0.0,  0.5,  1.0),
    point3d(0.5,  0.5,  1.0)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feQr3d<2, PMAPTYPE>::
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
feQr3d<2, PMAPTYPE>::
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
feQr3d<2, PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation
  feStaticRealEval<feLabel,numPoly>::eval(val,Y);
}

template<typename PMAPTYPE>
void
feQr3d<2, PMAPTYPE>::
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
