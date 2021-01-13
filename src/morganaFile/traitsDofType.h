/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef TRAITSDOFTYPE_H
#define TRAITSDOFTYPE_H

#include "typesInterface.hpp"
#include "printMesh.hpp"
#include "../morganaDofs/komplex.h"

#include "EpetraExt_DistArray.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "EpetraExt_HDF5.h"


using namespace std;


//-------------------------------------------------------------------------------------------------
/*! Generic traits - empty */
template<typename TYPE>
class traitsDofType;


//-------------------------------------------------------------------------------------------------
/*! Real trait */
template<> class traitsDofType<Real>
{
    /*! @name Typedefs */ //@{
  public:
    typedef EpetraExt::HDF5 HDF5;
    //@}  
  
    /*! @name Constructor */ //@{
  public:
    traitsDofType();
    //@}
    
    /*! @name Print to UCD */ //@{
  public:
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const Real & a);
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    template<typename OUTVECT>
    void printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap);
    void printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows);
    //@}
};


template<typename OUTVECT>
void
traitsDofType<Real>::
printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap)
{
  string s = "field";
  EpetraExt::DistArray<double> dataEpetra(epetraMap,1);
  
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  {
    dataEpetra(i-1,0) = dataMorgana(i);
  }
    
  hdf5.Write(s,dataEpetra);
}



//-------------------------------------------------------------------------------------------------
/*! Point2d trait */
template<> class traitsDofType<point2d>
{
    /*! @name Typedefs */ //@{
  public:
    typedef EpetraExt::HDF5 HDF5;
    //@}  
  
    /*! @name Constructor */ //@{
  public:
    traitsDofType();
    //@}
    
    /*! @name Print to UCD */ //@{
  public:
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const point2d & a);
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    template<typename OUTVECT>
    void printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap);
    void printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows);
    //@}
};


template<typename OUTVECT>
void
traitsDofType<point2d>::
printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap)
{
  string s = "field";
  EpetraExt::DistArray<double> dataEpetra(epetraMap,3);
  
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  {
    dataEpetra(i-1,0) = dataMorgana(i).getX();
    dataEpetra(i-1,1) = dataMorgana(i).getY();
    dataEpetra(i-1,2) = 0.0;
  }
    
  hdf5.Write(s,dataEpetra); 
}



//-------------------------------------------------------------------------------------------------
/*! Tensor2d trait */
template<> class traitsDofType<tensor2d>
{
    /*! @name Typedefs */ //@{
  public:
    typedef EpetraExt::HDF5 HDF5;
    //@}  
  
    /*! @name Constructor */ //@{
  public:
    traitsDofType();
    //@}
    
    /*! @name Print to UCD */ //@{
  public:
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const tensor2d & a);
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    template<typename OUTVECT>
    void printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap);
    void printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows);
    //@}
};


template<typename OUTVECT>
void
traitsDofType<tensor2d>::
printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap)
{
  string s;
  EpetraExt::DistArray<double> dataEpetra(epetraMap,1);

  // XX
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(1,1); }
  s = "fieldXX";
  hdf5.Write(s,dataEpetra);
  
  // XY
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(1,2); }
  s = "fieldXY";
  hdf5.Write(s,dataEpetra);
  
  
  // YX
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(2,1); }
  s = "fieldYX";
  hdf5.Write(s,dataEpetra);
  
  // YY
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(2,2); }
  s = "fieldYY";
  
  
  hdf5.Write(s,dataEpetra);
}



//-------------------------------------------------------------------------------------------------
/*! Point3d trait */
template<> class traitsDofType<point3d>
{
    /*! @name Typedefs */ //@{
  public:
    typedef EpetraExt::HDF5 HDF5;
    //@}
  
    /*! @name Constructor */ //@{
  public:
    traitsDofType();
    //@}
     
    /*! @name Print to UCD */ //@{
  public:
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const point3d & a);
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    template<typename OUTVECT>
    void printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap);
    void printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows);
    //@}
};


template<typename OUTVECT>
void
traitsDofType<point3d>::
printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap)
{
  string s = "field";
  EpetraExt::DistArray<double> dataEpetra(epetraMap,3);
  
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  {
    dataEpetra(i-1,0) = dataMorgana(i).getX();
    dataEpetra(i-1,1) = dataMorgana(i).getY();
    dataEpetra(i-1,2) = dataMorgana(i).getZ();
  }
    
  hdf5.Write(s,dataEpetra); 
}



//-------------------------------------------------------------------------------------------------
/*! Tensor3d trait */
template<> class traitsDofType<tensor3d>
{
    /*! @name Typedefs */ //@{
  public:
    typedef EpetraExt::HDF5 HDF5;
    //@}
  
    /*! @name Constructor */ //@{
  public:
    traitsDofType();
    //@}
    
    /*! @name Print to UCD */ //@{
  public:
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const tensor3d & a);
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    template<typename OUTVECT>
    void printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap);
    void printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows);
    //@}
};


template<typename OUTVECT>
void
traitsDofType<tensor3d>::
printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap)
{
  string s;
  EpetraExt::DistArray<double> dataEpetra(epetraMap,1);

  // XX
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(1,1); }
  s = "fieldXX";
  hdf5.Write(s,dataEpetra);
  
  // XY
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(1,2); }
  s = "fieldXY";
  hdf5.Write(s,dataEpetra);
  
  // XZ
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(1,3); }
  s = "fieldXZ";
  hdf5.Write(s,dataEpetra);
  
  
  // YX
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(2,1); }
  s = "fieldYX";
  hdf5.Write(s,dataEpetra);
  
  // YY
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(2,2); }
  s = "fieldYY";
  hdf5.Write(s,dataEpetra);
  
  // YZ
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(2,3); }
  s = "fieldYZ";
  hdf5.Write(s,dataEpetra);
  
  
  // ZX
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(3,1); }
  s = "fieldZX";
  hdf5.Write(s,dataEpetra);
  
  // ZY
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(3,2); }
  s = "fieldZY";
  hdf5.Write(s,dataEpetra);
  
  // ZZ
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getIJ(3,3); }
  s = "fieldZZ";
  hdf5.Write(s,dataEpetra);
}



//-------------------------------------------------------------------------------------------------
/*! Static vector trait */
template<size_t N> class traitsDofType<staticVector<N> >
{
    /*! @name Typedefs */ //@{
  public:
    typedef EpetraExt::HDF5 HDF5;
    //@}
  
    /*! @name Constructor */ //@{
  public:
    traitsDofType();
    //@}
   
    /*! @name Print to UCD */ //@{
  public:
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const staticVector<N> & a);
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    template<typename OUTVECT>
    void printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap);
    void printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows);
    //@}
};

template<size_t N>
traitsDofType<staticVector<N> >::
traitsDofType()
{ }

template<size_t N>
void
traitsDofType<staticVector<N> >::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << N << " ";
  
  for(UInt i=1; i <= N; ++i)
  {
    out << " " << 1;
  }
  out << endl;
  
  for(UInt i=1; i <= N; ++i)
  {
    out << value << i << ", " << unit << endl;
  }  
}

template<size_t N>
void
traitsDofType<staticVector<N> >::
printDof(ofstream & out, const staticVector<N> & a)
{
  for(UInt i=1; i <= N; ++i)
  {
    out << a(i) << " ";
  }
  out << endl;
}

template<size_t N>
template<typename OUTVECT>
void
traitsDofType<staticVector<N> >::
printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap)
{
  string s;
  EpetraExt::DistArray<double> dataEpetra(epetraMap,1);
  
  for(UInt k=1; k <= N; ++k)
  {
    for(UInt i=1; i <= dataMorgana.size(); ++i)
    {
      dataEpetra(i-1,0) = dataMorgana(i)(k);
    }
    
    s = "field" + num2str(k);
    hdf5.Write(s,dataEpetra);
  }  
}

template<size_t N>
void
traitsDofType<staticVector<N> >::
printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows)
{
  string s;
  
  for(UInt k=1; k <= N; ++k)
  {
    s = "field" + num2str(k);
    
    xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
    xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
    xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
    xdmfFile << "        </DataItem>"                                                                               << endl;
    xdmfFile << "      </Attribute>"                                                                                << endl;
  }    
}



//-------------------------------------------------------------------------------------------------
/*! Static komplex trait */
template<> class traitsDofType<komplex>
{
    /*! @name Typedefs */ //@{
  public:
    typedef EpetraExt::HDF5 HDF5;
    //@}
  
    /*! @name Constructor */ //@{
  public:
    traitsDofType();
    //@}
   
    /*! @name Print to UCD */ //@{
  public:
    void paraviewString(ofstream & out, const string & value, const string & unit);
    void printDof(ofstream & out, const komplex & a);
    //@}
    
    /*! @name Print HDF5 */ //@{
  public:
    template<typename OUTVECT>
    void printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap);
    void printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows);
    //@}
};


template<typename OUTVECT>
void
traitsDofType<komplex>::
printHDF5(HDF5 & hdf5, const OUTVECT & dataMorgana, const Epetra_Map & epetraMap)
{
  string s;
  EpetraExt::DistArray<double> dataEpetra(epetraMap,1);
  
  //Real component
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getReal(); }
  s = "real";
  hdf5.Write(s,dataEpetra);
  
  //Imag component
  for(UInt i=1; i <= dataMorgana.size(); ++i)
  { dataEpetra(i-1,0) = dataMorgana(i).getImag(); }
  s = "imag";
  hdf5.Write(s,dataEpetra);
}


#endif
