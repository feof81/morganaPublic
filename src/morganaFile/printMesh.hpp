/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef PRINTMESH_HPP
#define PRINTMESH_HPP

#include "load.h"
#include "traitsMapItemFixer.hpp"

#include "pVect.hpp"
#include "pGraph.hpp"

#include "geoShapes.h"
#include "mesh1d.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"

using namespace std;


//_______________________________________________________________________________________________________
// A USEFULL FUNCTION
//-------------------------------------------------------------------------------------------------------
inline string num2str(const int &i)
{
  stringstream ss;
  string strNum;
  
  ss << i;
  ss >> strNum;
  
  // MASSIMO NUMERO DI ITERAZIONI CONSENTITE: 100000
  if (i < 10)
    strNum = "0000" + strNum;
  else if (i < 100)
    strNum = "000"  + strNum;
  else if (i < 1000)
    strNum = "00"   + strNum;
  else if (i < 10000)
    strNum = "0"    + strNum;
  else if (i < 100000)
    ;
  else 
    cout << "ERRORE: Massimo numero di iterazioni per plottaggio 100000" << endl;
  
  return(strNum);
}


//_______________________________________________________________________________________________________
// NOT SPECIALIZED
//-------------------------------------------------------------------------------------------------------
/*! Not specialized */
template<typename GEOSHAPE, typename ELMAP, typename NODEMAP>
class printMesh
{
};



//_______________________________________________________________________________________________________
// LINEAR LINE SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Print Mesh - Line specialization */
template<typename ELMAP, typename NODEMAP>
class printMesh<linearLine,ELMAP,NODEMAP>
{
  /*! @name Typedefs */ //@{
  public:
    typedef mesh1d<linearLine,ELMAP,NODEMAP>  MESH1D;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    static const Real paraFix;
    printMesh();
    //@}
    
    /*! @name Print functions */ //@{
  public:
    void paraviewSerial(const string & s, const MESH1D & Mesh1d);
    void paraviewLocal(const string & s, const UInt & pid, const MESH1D & Mesh1d);
    //@}
};

template<typename ELMAP, typename NODEMAP>
const Real printMesh<linearLine,ELMAP,NODEMAP>::paraFix = 1.000000001;

template<typename ELMAP, typename NODEMAP>
printMesh<linearLine,ELMAP,NODEMAP>::
printMesh()
{
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearLine,ELMAP,NODEMAP>::
paraviewSerial(const string & s, const MESH1D & Mesh1d)
{
  ofstream out(s.c_str());

  out << Mesh1d.getNumVertices() << " " << Mesh1d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh1d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh1d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh1d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh1d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  {
    out << i << " ";
    out << Mesh1d.getElementL(i).getGeoId() << "  ";
    out << "line ";
    
    assert(Mesh1d.getElementL(i).size() == 2);
    
    for(UInt k=1; k <= Mesh1d.getElementL(i).size(); ++k)
    {
      out << Mesh1d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearLine,ELMAP,NODEMAP>::
paraviewLocal(const string & s, const UInt & pid, const MESH1D & Mesh1d)
{
  ofstream out(s.c_str());

  out << Mesh1d.getNumVertices() << " " << Mesh1d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh1d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh1d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh1d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh1d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh1d.getNumElements(); ++i)
  {
    out << i << " ";
    out << pid << "  ";
    out << "line ";
    
    assert(Mesh1d.getElementL(i).size() == 2);
    
    for(UInt k=1; k <= Mesh1d.getElementL(i).size(); ++k)
    {
      out << Mesh1d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}



//_______________________________________________________________________________________________________
// LINEAR TRIANGLE SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Print Mesh - Triangle specialization */
template<typename ELMAP, typename NODEMAP>
class printMesh<linearTriangle,ELMAP,NODEMAP>
{
  /*! @name Typedefs */ //@{
  public:
    typedef mesh2d<linearTriangle,ELMAP,NODEMAP> MESH2D;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    static const Real paraFix;
    printMesh();
    //@}
    
    /*! @name Print functions */ //@{
  public:
    void paraviewSerial(const string & s, const MESH2D & Mesh2d);
    void paraviewLocal(const string & s, const UInt & pid, const MESH2D & Mesh2d);
    //@}
};

template<typename ELMAP, typename NODEMAP>
const Real printMesh<linearTriangle,ELMAP,NODEMAP>::paraFix = 1.000000001;

template<typename ELMAP, typename NODEMAP>
printMesh<linearTriangle,ELMAP,NODEMAP>::
printMesh()
{
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearTriangle,ELMAP,NODEMAP>::
paraviewSerial(const string & s, const MESH2D & Mesh2d)
{
  ofstream out(s.c_str());

  out << Mesh2d.getNumVertices() << " " << Mesh2d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh2d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh2d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    out << i << " ";
    out << Mesh2d.getElementL(i).getGeoId() << "  ";
    out << "tri ";
    
    assert(Mesh2d.getElementL(i).size() == 3);
    
    for(UInt k=1; k <= Mesh2d.getElementL(i).size(); ++k)
    {
      out << Mesh2d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearTriangle,ELMAP,NODEMAP>::
paraviewLocal(const string & s, const UInt & pid, const MESH2D & Mesh2d)
{
  ofstream out(s.c_str());

  out << Mesh2d.getNumVertices() << " " << Mesh2d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh2d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh2d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    out << i << " ";
    out << pid << "  ";
    out << "tri ";
    
    assert(Mesh2d.getElementL(i).size() == 3);
    
    for(UInt k=1; k <= Mesh2d.getElementL(i).size(); ++k)
    {
      out << Mesh2d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}



//_______________________________________________________________________________________________________
// LINEAR QUAD SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Print Mesh - Quad specialization */
template<typename ELMAP, typename NODEMAP>
class printMesh<linearQuad,ELMAP,NODEMAP>
{
  /*! @name Typedefs */ //@{
  public:
    typedef mesh2d<linearQuad,ELMAP,NODEMAP> MESH2D;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    static const Real paraFix;
    printMesh();
    //@}
    
    /*! @name Print functions */ //@{
  public:
    void paraviewSerial(const string & s, const MESH2D & Mesh2d);
    void paraviewLocal(const string & s, const UInt & pid, const MESH2D & Mesh2d);
    //@}
};

template<typename ELMAP, typename NODEMAP>
const Real printMesh<linearQuad,ELMAP,NODEMAP>::paraFix = 1.000000001;

template<typename ELMAP, typename NODEMAP>
printMesh<linearQuad,ELMAP,NODEMAP>::
printMesh()
{
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearQuad,ELMAP,NODEMAP>::
paraviewSerial(const string & s, const MESH2D & Mesh2d)
{
  ofstream out(s.c_str());

  out << Mesh2d.getNumVertices() << " " << Mesh2d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh2d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh2d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    out << i << " ";
    out << Mesh2d.getElementL(i).getGeoId() << "  ";
    out << "quad ";
    
    assert(Mesh2d.getElementL(i).size() == 4);
    
    for(UInt k=1; k <= Mesh2d.getElementL(i).size(); ++k)
    {
      out << Mesh2d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearQuad,ELMAP,NODEMAP>::
paraviewLocal(const string & s, const UInt & pid, const MESH2D & Mesh2d)
{
  ofstream out(s.c_str());

  out << Mesh2d.getNumVertices() << " " << Mesh2d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh2d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh2d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh2d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh2d.getNumElements(); ++i)
  {
    out << i << " ";
    out << pid << "  ";
    out << "quad ";
    
    assert(Mesh2d.getElementL(i).size() == 4);
    
    for(UInt k=1; k <= Mesh2d.getElementL(i).size(); ++k)
    {
      out << Mesh2d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}



//_______________________________________________________________________________________________________
// LINEAR TETRA SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Print Mesh - Tetra specialization */
template<typename ELMAP, typename NODEMAP>
class printMesh<linearTetra,ELMAP,NODEMAP>
{
  /*! @name Typedefs */ //@{
  public:
    typedef mesh3d<linearTetra,ELMAP,NODEMAP> MESH3D;
    typedef typename MESH3D::NODESVECT        NODESVECT;
    typedef typename MESH3D::GRAPH3D          GRAPH3D;
    typedef typename MESH3D::GRAPH2D          GRAPH2D;
    //@}
    
    /*! @name Constructors and values */ //@{
  public:
    static const Real paraFix;
    printMesh();
    //@}
    
    /*! @name Print functions */ //@{
  public:
    void paraviewSerial(const string & s, const MESH3D & Mesh3d);
    void paraviewLocal(const string & s, const UInt & pid, const MESH3D & Mesh3d);
    //@}
};

template<typename ELMAP, typename NODEMAP>
const Real printMesh<linearTetra,ELMAP,NODEMAP>::paraFix = 1.000000001;

template<typename ELMAP, typename NODEMAP>
printMesh<linearTetra,ELMAP,NODEMAP>::
printMesh()
{
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearTetra,ELMAP,NODEMAP>::
paraviewSerial(const string & s, const MESH3D & Mesh3d)
{
  ofstream out(s.c_str());

  out << Mesh3d.getNumVertices() << " " << Mesh3d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh3d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh3d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    out << i << " ";
    out << Mesh3d.getElementL(i).getGeoId() << "  ";
    out << "tet ";
    
    assert(Mesh3d.getElementL(i).size() == 4);
    
    for(UInt k=1; k <= Mesh3d.getElementL(i).size(); ++k)
    {
      out << Mesh3d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearTetra,ELMAP,NODEMAP>::
paraviewLocal(const string & s, const UInt & pid, const MESH3D & Mesh3d)
{
  ofstream out(s.c_str());

  out << Mesh3d.getNumVertices() << " " << Mesh3d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh3d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh3d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getZ() << endl;
  }
  
  // passo alle facce
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    out << i << " ";
    out << pid << "  ";
    out << "tet ";
    
    assert(Mesh3d.getElementL(i).size() == 4);
    
    for(UInt k=1; k <= Mesh3d.getElementL(i).size(); ++k)
    {
      out << Mesh3d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}



//_______________________________________________________________________________________________________
// LINEAR HEXAHEDRA SPECIALIZATION
//-------------------------------------------------------------------------------------------------------
/*! Print Mesh - Hexa specialization */
template<typename ELMAP, typename NODEMAP>
class printMesh<linearHexa,ELMAP,NODEMAP>
{
  /*! @name Typedefs */ //@{
  public:
    typedef mesh3d<linearHexa,ELMAP,NODEMAP>  MESH3D;
    typedef typename MESH3D::NODESVECT        NODESVECT;
    typedef typename MESH3D::GRAPH3D          GRAPH3D;
    typedef typename MESH3D::GRAPH2D          GRAPH2D;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    static const Real paraFix;
    printMesh();
    //@}
    
    /*! @name Print functions */ //@{
  public:
    void paraviewSerial(const string & s, const MESH3D & Mesh3d);
    void paraviewLocal(const string & s, const UInt & pid, const MESH3D & Mesh3d);
    //@}
};

template<typename ELMAP, typename NODEMAP>
const Real printMesh<linearHexa,ELMAP,NODEMAP>::paraFix = 1.000000001;

template<typename ELMAP, typename NODEMAP>
printMesh<linearHexa,ELMAP,NODEMAP>::
printMesh()
{
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearHexa,ELMAP,NODEMAP>::
paraviewSerial(const string & s, const MESH3D & Mesh3d)
{
  ofstream out(s.c_str());

  out << Mesh3d.getNumVertices() << " " << Mesh3d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh3d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh3d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    out << i << " ";
    out << Mesh3d.getElementL(i).getGeoId() << "  ";
    out << "hex ";
    
    assert(Mesh3d.getElementL(i).size() == 8);
    
    for(UInt k=1; k <= Mesh3d.getElementL(i).size(); ++k)
    {
      out << Mesh3d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}

template<typename ELMAP, typename NODEMAP>
void
printMesh<linearHexa,ELMAP,NODEMAP>::
paraviewLocal(const string & s, const UInt & pid, const MESH3D & Mesh3d)
{
  ofstream out(s.c_str());

  out << Mesh3d.getNumVertices() << " " << Mesh3d.getNumElements() << " 0 0 0" << endl;
  
  for(UInt i=1; i<= Mesh3d.getNumVertices(); ++i)
  {
    out << i << " ";
    out << paraFix * Mesh3d.getNodeL(i).getX() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getY() << " ";
    out << paraFix * Mesh3d.getNodeL(i).getZ() << endl;
  }
  
  for(UInt i=1; i <= Mesh3d.getNumElements(); ++i)
  {
    out << i << " ";
    out << pid << "  ";
    out << "hex ";
    
    assert(Mesh3d.getElementL(i).size() == 8);
    
    for(UInt k=1; k <= Mesh3d.getElementL(i).size(); ++k)
    {
      out << Mesh3d.getElementL(i).getCid(k) << "  ";
    }
    out << endl;
  }
}

#endif
