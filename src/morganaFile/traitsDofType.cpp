/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "traitsDofType.h"


//-------------------------------------------------------------------------------------------------
/*! Real trait */
traitsDofType<Real>::
traitsDofType()
{ }

void
traitsDofType<Real>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "1 1" << endl;
  out << value << ", " << unit << endl;
}

void
traitsDofType<Real>::
printDof(ofstream & out, const Real & a)
{
  out << a << endl;
}

void
traitsDofType<Real>::
printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows)
{
  string s = "field";

  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;   
}



//-------------------------------------------------------------------------------------------------
/*! Point2d trait */
traitsDofType<point2d>::
traitsDofType()
{ }

void
traitsDofType<point2d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "1 3" << endl;
  out << value << ", " << unit << endl;
}

void
traitsDofType<point2d>::
printDof(ofstream & out, const point2d & a)
{
  out << a.getX() << " " << a.getY() << " " << 0 << endl;
}

void
traitsDofType<point2d>::
printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows)
{
  string s = "field";

  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Vector\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 3 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
}



//-------------------------------------------------------------------------------------------------
/*! Tensor2d trait */
traitsDofType<tensor2d>::
traitsDofType()
{ }

void
traitsDofType<tensor2d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "4 1 1 1 1"  << endl;
  out << value << "xx" << ", " << unit << endl;
  out << value << "xy" << ", " << unit << endl;
  out << value << "yx" << ", " << unit << endl;
  out << value << "yy" << ", " << unit << endl;
}

void
traitsDofType<tensor2d>::
printDof(ofstream & out, const tensor2d & a)
{
  out << a.getIJ(1,1) << " " << a.getIJ(1,2) << " " << a.getIJ(2,1) << " " << a.getIJ(2,2) << endl;
}

void
traitsDofType<tensor2d>::
printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows)
{
  string s;
  
  s = "fieldXX";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
 
  s = "fieldXY";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
  
  s = "fieldYX";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
 
  s = "fieldYY";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
}



//-------------------------------------------------------------------------------------------------
/*! Point3d trait */
traitsDofType<point3d>::
traitsDofType()
{ }

void
traitsDofType<point3d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "1 3" << endl;
  out << value << ", " << unit << endl;
}

void
traitsDofType<point3d>::
printDof(ofstream & out, const point3d & a)
{
  out << a.getX() << " " << a.getY() << " " << a.getZ() << endl;
}

void
traitsDofType<point3d>::
printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows)
{
  string s = "field";

  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Vector\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 3 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
}



//-------------------------------------------------------------------------------------------------
/*! Tensor3d trait */
traitsDofType<tensor3d>::
traitsDofType()
{ }

void
traitsDofType<tensor3d>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << "9 1 1 1 1 1 1 1 1 1"  << endl;
  out << value << "xx" << ", " << unit << endl;
  out << value << "xy" << ", " << unit << endl;
  out << value << "xz" << ", " << unit << endl;
  out << value << "yx" << ", " << unit << endl;
  out << value << "yy" << ", " << unit << endl;
  out << value << "yz" << ", " << unit << endl;
  out << value << "zx" << ", " << unit << endl;
  out << value << "zy" << ", " << unit << endl;
  out << value << "zz" << ", " << unit << endl;
}

void
traitsDofType<tensor3d>::
printDof(ofstream & out, const tensor3d & a)
{
  out << a.getIJ(1,1) << " " << a.getIJ(1,2) << " " << a.getIJ(1,3)
      << a.getIJ(2,1) << " " << a.getIJ(2,2) << " " << a.getIJ(2,3)
      << a.getIJ(3,1) << " " << a.getIJ(3,2) << " " << a.getIJ(3,3) << endl;
}

void
traitsDofType<tensor3d>::
printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows)
{
  string s;
  
  s = "fieldXX";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
 
  s = "fieldXY";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
  
  s = "fieldXZ";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
  
  s = "fieldYX";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
 
  s = "fieldYY";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
  
  s = "fieldYZ";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
  
  s = "fieldZX";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
 
  s = "fieldZY";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
  
  s = "fieldZZ";  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
}



//-------------------------------------------------------------------------------------------------
/*! Komplex */
traitsDofType<komplex>::
traitsDofType()
{ }

void
traitsDofType<komplex>::
paraviewString(ofstream & out, const string & value, const string & unit)
{
  out << 2 << " ";
  
  out << " " << 1;
  out << " " << 1;
  out << endl;
  
  out << value << 1 << ", " << unit << endl;
  out << value << 2 << ", " << unit << endl;  
}

void
traitsDofType<komplex>::
printDof(ofstream & out, const komplex & a)
{
  out << a.getReal() << " ";
  out << a.getImag() << " ";
  out << endl;
}

void
traitsDofType<komplex>::
printXDMF(ofstream & xdmfFile, const string & center, const string & fileName, const UInt & numRows)
{
  string s;
  
  //Real
  s = "real";
  
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;
  
  //Complex
  s = "imag";
    
  xdmfFile << "      <Attribute Name=\"" << s << "\" AttributeType=\"Scalar\"" << " Center=\"" + center + "\">"   << endl;
  xdmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Dimensions=\"" << numRows << " " << 1 << "\">" << endl;
  xdmfFile << "          " << fileName << ".h5:/" << s << "/Values"                                               << endl;
  xdmfFile << "        </DataItem>"                                                                               << endl;
  xdmfFile << "      </Attribute>"                                                                                << endl;  
}
