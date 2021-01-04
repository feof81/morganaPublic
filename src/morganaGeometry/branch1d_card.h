/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef BRANCH1D_CARD_H
#define BRANCH1D_CARD_H

#include "typesInterface.hpp"
#include "simpleFormats.hpp"


class branch1d_card
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt    nodeGid;
    UInt    intElementGid;
    point3d intNormal;
    UInt           extNodeGid;
    sVect<UInt>    extElementsGid;
    sVect<point3d> extNormals;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    branch1d_card();
    branch1d_card(const branch1d_card & card);
    branch1d_card & operator=(const branch1d_card & card);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setNodeGid(const UInt & NodeGid);
    void setIntElementGid(const UInt & IntElementGid);
    void setIntNormal(const point3d & IntNormal);
    void setExtNodeGid(const UInt & ExtNodeGid);
    void setExtElementsGid(const sVect<UInt> & ExtElementsGid);
    void setExtNormals(const sVect<point3d> & ExtNormals);
    //@}
    
    /*! @name Add functions */ //@{
  public:
    void addExtElementsGid(const UInt & gid);
    void addExtNormals(const point3d & P);
    //@}
    
    /*! @name Get end nodes functions - const */ //@{
  public:
    const UInt           & getNodeGid() const;
    const UInt           & getIntElementGid() const;
    const point3d        & getIntNormal() const;
    const UInt           & getExtNodeGid() const;
          UInt             getExtElementsNum() const;
    const UInt           & getExtElementsGid(const UInt & i) const;
    const sVect<UInt>    & getExtElementsGid() const;
          UInt             getExtNormalsNum() const;
    const point3d        & getExtNormals(const UInt & i) const;
    const sVect<point3d> & getExtNormals() const;
    //@}
    
    /*! @name Get end nodes functions */ //@{
  public:
    UInt           & getNodeGid();
    UInt           & getIntElementGid();
    point3d        & getIntNormal();
    UInt           & getExtNodeGid();
    UInt           & getExtElementsGid(const UInt & i);
    sVect<UInt>    & getExtElementsGid();
    point3d        & getExtNormals(const UInt & i);
    sVect<point3d> & getExtNormals();
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    friend ostream & operator<<(ostream & f, const branch1d_card & card);
    //@}
};


//_________________________________________________________________________________________________
// TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
branch1d_card::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & nodeGid;
  ar & intElementGid;
  ar & intNormal;
  ar & extNodeGid;
  extElementsGid.serialize(ar,version);
  extNormals.serialize(ar,version);
}

    
#endif
