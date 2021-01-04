/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef BRANCH1D_SEARCHCARD_H
#define BRANCH1D_SEARCHCARD_H

#include "typesInterface.hpp"
#include "simpleFormats.hpp"


class branch1d_searchCard
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt    tgtNodeGid;
    point3d tgtP;
    UInt    tgtElGid;
    point3d tgtNormal;
    
    UInt srcNodeGid;
    sVect<UInt>    srcElementsGids;
    sVect<point3d> srcNormals;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    branch1d_searchCard();
    branch1d_searchCard(const branch1d_searchCard & card);
    branch1d_searchCard & operator=(const branch1d_searchCard & card);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setTgtNodeGid(const UInt & TgtNodeGid);
    void setTgtNode(const point3d & Q);
    void setTgtElGid(const UInt & TgtElGid);
    void setTgtNormal(const point3d & TgtNormal);    
    void setSrcNodeGid(const UInt & SrcNodeGid);
    void setSrcElementsGids(const sVect<UInt> & SrcElementsGids);
    void setSrcNormals(const sVect<point3d> & SrcNormals);
    //@}
    
    /*! @name Add functions */ //@{
  public:
    void addSrcElementsGids(const UInt & gid);
    void addSrcNormals(const point3d & P);
    //@}
    
    /*! @name Get functions - const */ //@{
  public:
    const UInt    & getTgtNodeGid() const;
    const point3d & getTgtNode() const;
    const UInt    & getTgtElGid() const;
    const point3d & getTgtNormal() const;
    
    const UInt           & getSrcNodeGid() const;
          UInt             getSrcElementsNum() const;
    const UInt           & getSrcElementsGids(const UInt & i) const;
    const sVect<UInt>    & getSrcElementsGids() const;
          UInt             getSrcNormalsNum() const;
    const point3d        & getSrcNormals(const UInt & i) const;
    const sVect<point3d> & getSrcNormals() const;
    //@}
    
    /*! @name Get functions */ //@{
  public:
    UInt    & getTgtNodeGid();
    point3d & getTgtNode();
    UInt    & getTgtElGid();
    point3d & getTgtNormal();
    
    UInt           & getSrcNodeGid();
    sVect<UInt>    & getSrcElementsGids();
    UInt           & getSrcElementsGids(const UInt & i);
    sVect<point3d> & getSrcNormals();
    point3d        & getSrcNormals(const UInt & i);
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    friend ostream & operator<<(ostream & f, const branch1d_searchCard & card);
    //@}
};


//_________________________________________________________________________________________________
// TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
branch1d_searchCard::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & tgtNodeGid;
  ar & tgtP;
  ar & tgtElGid;
  ar & tgtNormal;
  ar & srcNodeGid;
  srcElementsGids.serialize(ar,version);
  srcNormals.serialize(ar,version);
}


#endif
