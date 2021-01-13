/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEMAPPIOLA2D_HPP
#define FEMAPPIOLA2D_HPP

#include "elCard2d.hpp"
#include "geoMapInterface.hpp"


/*! The two dimensional Piola-Kirkoff map */
template<typename GEOSHAPE, typename PMAPTYPE>
class feMapPiola2d
{
    /*! @name Typedefs */ //@{
  public:
    typedef geoMapInterface<GEOSHAPE>   GEOMAP;
    typedef elCard2d<GEOSHAPE,PMAPTYPE> ELCARD;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    GEOMAP geoMap;
    ELCARD ElCard;
    bool cardLoaded;
    //@}
    
    /*! @name Constructors and set functions */ //@{
  public:
    feMapPiola2d();
    feMapPiola2d(const ELCARD & ELCard);
    void setCard(const ELCARD & ELCard);
    //@}
    
    /*! @name Internal computing functions */ //@{
  public:
    /*! Compute the gradient of the jacobian */
    tensor3d computeFgradX(const point3d & Y) const;
    tensor3d computeFgradY(const point3d & Y) const;
    tensor3d computeFgradZ(const point3d & Y) const;
    
    /*! Compute the gradient of the determinant of the jacobian */
    Real computeFgradDet(const tensor3d & T, const tensor3d & TG) const;
    //@}
    
    /*! @name Eval functions */ //@{
  public:   
    /*! Eval F(V) with F the Piola transformation 
    \param V input vector 
    \param Y local coordinate */
    point3d eval(const point3d & V, const point3d & Y) const;
    
    /*! Gradient of F(V) with F the Piola transformation - \c V vector \c W gradient vector */
    point3d evalGradX(const point3d & V, const point3d & W, const point3d & Y) const;
    point3d evalGradY(const point3d & V, const point3d & W, const point3d & Y) const;
    point3d evalGradZ(const point3d & V, const point3d & W, const point3d & Y) const;
    //@}
};


template<typename GEOSHAPE, typename PMAPTYPE>
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
feMapPiola2d()
{
  cardLoaded = false;
}

template<typename GEOSHAPE, typename PMAPTYPE>
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
feMapPiola2d(const ELCARD & ELCard)
{
  cardLoaded = true;
  ElCard     = ELCard;
}

template<typename GEOSHAPE, typename PMAPTYPE>
void
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
setCard(const ELCARD & ELCard)
{
  cardLoaded = true;
  ElCard     = ELCard;
}

template<typename GEOSHAPE, typename PMAPTYPE>
tensor3d
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
computeFgradX(const point3d & Y) const
{
  assert(cardLoaded);
  
  tensor3d TX;
  
  TX.setCol(1, geoMap.template getDerivative <2,0,0> (ElCard.getNodes(),Y) );
  TX.setCol(2, geoMap.template getDerivative <1,1,0> (ElCard.getNodes(),Y) );
  TX.setCol(3, geoMap.template getDerivative <1,0,1> (ElCard.getNodes(),Y) );
  
  if( point3d::norm2(TX.getCol(1) ^ TX.getCol(2)) != 0.0 )
  { TX.completeThirdColoumn(); }
    
  return(TX);
}

template<typename GEOSHAPE, typename PMAPTYPE>
tensor3d
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
computeFgradY(const point3d & Y) const
{
  assert(cardLoaded);
  
  tensor3d TY;
  
  TY.setCol(1, geoMap.template getDerivative <1,1,0> (ElCard.getNodes(),Y) );
  TY.setCol(2, geoMap.template getDerivative <0,2,0> (ElCard.getNodes(),Y) );
  TY.setCol(3, geoMap.template getDerivative <0,1,1> (ElCard.getNodes(),Y) );
  
  if( point3d::norm2(TY.getCol(1) ^ TY.getCol(2)) != 0.0 )
  { TY.completeThirdColoumn(); }
  
  return(TY);
}

template<typename GEOSHAPE, typename PMAPTYPE>
tensor3d
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
computeFgradZ(const point3d & Y) const
{
  assert(cardLoaded);
  
  tensor3d TZ;
  
  TZ.setCol(1, geoMap.template getDerivative <1,0,1> (ElCard.getNodes(),Y) );
  TZ.setCol(2, geoMap.template getDerivative <0,1,1> (ElCard.getNodes(),Y) );
  TZ.setCol(3, geoMap.template getDerivative <0,0,2> (ElCard.getNodes(),Y) );
  
  if( point3d::norm2(TZ.getCol(1) ^ TZ.getCol(2)) != 0.0 )
  { TZ.completeThirdColoumn(); }
  
  return(TZ);
}

template<typename GEOSHAPE, typename PMAPTYPE>
Real
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
computeFgradDet(const tensor3d & T, const tensor3d & TG) const
{
  assert(cardLoaded);  
  tensor3d G1(T), G2(T), G3(T);
  
  G1.setCol(1,TG.getCol(1));
  G2.setCol(2,TG.getCol(2));
  G3.setCol(3,TG.getCol(3));
  
  return( G1.getThirdInvariant() + G2.getThirdInvariant() + G3.getThirdInvariant() );
}

template<typename GEOSHAPE, typename PMAPTYPE>
point3d
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
eval(const point3d & V, const point3d & Y) const
{
  assert(cardLoaded);
  
  tensor3d T = geoMap.getGradient(ElCard.getNodes(),Y);
  T.completeThirdColoumn();
  
  return(T.secondIndexSaturation(V) / T.getThirdInvariant());
}

template<typename GEOSHAPE, typename PMAPTYPE>
point3d
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
evalGradX(const point3d & V, const point3d & W, const point3d & Y) const
{
  assert(cardLoaded);
  
  tensor3d  T  = geoMap.getGradient(ElCard.getNodes(),Y);
  T.completeThirdColoumn();
  
  tensor3d  G  = computeFgradX(Y);
  Real detGrad = computeFgradDet(T,G);
  Real    det  = T.getThirdInvariant();
  
  G = T * (- (detGrad / pow(det,2.0))) + G / det;  
  return(G.secondIndexSaturation(V) + (T.secondIndexSaturation(W) / det));
}

template<typename GEOSHAPE, typename PMAPTYPE>
point3d
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
evalGradY(const point3d & V, const point3d & W, const point3d & Y) const
{
  assert(cardLoaded);
  
  tensor3d  T  = geoMap.getGradient(ElCard.getNodes(),Y);
  T.completeThirdColoumn();
  
  tensor3d  G  = computeFgradY(Y);
  Real detGrad = computeFgradDet(T,G);
  Real    det  = T.getThirdInvariant();
  
  G = T * (- (detGrad / pow(det,2.0))) + G / det;
  return(G.secondIndexSaturation(V) + (T.secondIndexSaturation(W) / det));
}

template<typename GEOSHAPE, typename PMAPTYPE>
point3d
feMapPiola2d<GEOSHAPE,PMAPTYPE>::
evalGradZ(const point3d & V, const point3d & W, const point3d & Y) const
{
  assert(cardLoaded);
  
  tensor3d  T  = geoMap.getGradient(ElCard.getNodes(),Y);
  T.completeThirdColoumn();
  
  tensor3d  G  = computeFgradZ(Y);
  Real detGrad = computeFgradDet(T,G);
  Real    det  = T.getThirdInvariant();
  
  G = T * (- (detGrad / pow(det,2.0))) + G / det;
  return(G.secondIndexSaturation(V) + (T.secondIndexSaturation(W) / det));
}

#endif
