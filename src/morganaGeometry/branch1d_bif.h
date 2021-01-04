/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef BRANCH1D_BIF_H
#define BRANCH1D_BIF_H

#include "typesInterface.hpp"
#include "simpleFormats.hpp"


class branch1d_bif
{
    /*! @name Parallel support */ //@{
  public:
    friend class boost::serialization::access;
    
    template<class ARK>
    void serialize(ARK & ar, const unsigned int version);
    //@}
    
    /*! @name Internal data */ //@{
  public:
    UInt           nodeGid;
    sVect<UInt>    intElementsGid;
    sVect<point3d> intNormals;
    
    UInt    extNodeGid;
    UInt    extElementGid;
    point3d extNormal;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    branch1d_bif();
    branch1d_bif(const branch1d_bif & card);
    branch1d_bif & operator=(const branch1d_bif & card);
    //@}
    
    /*! @name Set functions */ //@{
  public:
    void setNodeGid(const UInt & NodeGid);
    void setIntElementsGid(const sVect<UInt> & IntElementsGid);
    void setIntNormal(const sVect<point3d> & IntNormal);
    void setExtNodeGid(const UInt & ExtNodeGid);
    void setExtElementGid(const UInt & ExtElementGid);
    void setExtNormal(const point3d & ExtNormal);
    //@}
    
    /*! @name Add functions */ //@{
  public:
    void addIntElementsGid(const UInt & gid);
    void addIntNormals(const point3d & P);
    //@}
    
    /*! @name Get end nodes functions - const */ //@{
  public:
    const UInt           & getNodeGid() const;
          UInt             getIntElementsGidNum() const;  
    const sVect<UInt>    & getIntElementsGid() const;
          UInt             getIntNormalsNum() const;
    const sVect<point3d> & getIntNormals() const;
    const UInt           & getExtNodeGid() const;
    const UInt           & getExtElementGid() const;
    const point3d        & getExtNormal() const;
    //@}
    
    /*! @name Get end nodes functions - const */ //@{
  public:
    UInt           & getNodeGid(); 
    sVect<UInt>    & getIntElementsGid();
    sVect<point3d> & getIntNormals();
    UInt           & getExtNodeGid();
    UInt           & getExtElementGid();
    point3d        & getExtNormal();
    //@}
    
    /*! @name Outstream operators */ //@{
  public:
    friend ostream & operator<<(ostream & f, const branch1d_bif & card);
    //@}
};


//_________________________________________________________________________________________________
// TEMPLATED FUNCTIONS
//-------------------------------------------------------------------------------------------------
template<class ARK>
void
branch1d_bif::
serialize(ARK & ar, const unsigned int version)
{
  assert(version == version);
  
  ar & nodeGid;
  intElementsGid.serialize(ar,version);
  intNormals.serialize(ar,version);
  ar & extNodeGid;
  ar & extElementGid;
  ar & extNormal;
}


#endif
