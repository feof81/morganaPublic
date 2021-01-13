/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef LOAD_H
#define LOAD_H

#include "typesInterface.hpp"


/*! Class used to read and extract data from formatted files */
class load
{
    /*! @name Variables */ //@{
  public:
    string strProjects, strGeometries, strResumePath;
    vector<string> VS;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    load();
    //@}
    
    /*! @name Read commands */ //@{
  public:
    void controllo(ifstream & file, string & file_name);
    UInt leggiriga(ifstream & file);
    //@}
    
    /*! @name Extract information */ //@{
  public:
    void assegnavalore(const UInt &i, UInt & val);
    void assegnavalore(const UInt &i, Real & val);
    void assegnavalore(const UInt &i, Real * val);
    void assegnavalore(const UInt &i, string & val);
    void assegnavalore(const UInt &i, char   & val);
    //@}
    
    /*! @name Single line info */ //@{
  public:
    string readString(ifstream & file);
    string readFirstString(ifstream & file);
    UInt   readFirstUInt(ifstream & file);
    int    readFirstInt(ifstream & file);
    Real   readFirstReal(ifstream & file);
    //@}
    
    /*! @name Service functions */ //@{
  public:
    void caseinsensitive();
    void stoupper(string & s);
    void stolower(string & s);
    //@}
    
    /*! @name Segmentation functions */ //@{
  public:
    UInt segmentationUpper(const UInt & pid, const UInt & npid, const UInt & n) const;
    UInt segmentationLower(const UInt & pid, const UInt & npid, const UInt & n) const;
    //@}
};

#endif
