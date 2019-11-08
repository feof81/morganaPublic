#include "tensor3d.h"

//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
tensor3d::
tensor3d() : id(1)
{
  T[0][0] = 0.0;
  T[0][1] = 0.0;
  T[0][2] = 0.0;
  
  T[1][0] = 0.0;
  T[1][1] = 0.0;
  T[1][2] = 0.0;
  
  T[2][0] = 0.0;
  T[2][1] = 0.0;
  T[2][2] = 0.0;
}

tensor3d::
tensor3d(const Real & XX, const Real & XY, const Real & XZ,
	     const Real & YX, const Real & YY, const Real & YZ,
	     const Real & ZX, const Real & ZY, const Real & ZZ) : id(1)
{
  T[0][0] = XX;
  T[0][1] = XY;
  T[0][2] = XZ;
  
  T[1][0] = YX;
  T[1][1] = YY;
  T[1][2] = YZ;
  
  T[2][0] = ZX;
  T[2][1] = ZY;
  T[2][2] = ZZ;
}

tensor3d::
tensor3d(const point3d & P1, const point3d & P2, const point3d & P3) : id(1)
{
  T[0][0] = P1.X[0];
  T[0][1] = P2.X[0];
  T[0][2] = P3.X[0];
  
  T[1][0] = P1.X[1];
  T[1][1] = P2.X[1];
  T[1][2] = P3.X[1];
  
  T[2][0] = P1.X[2];
  T[2][1] = P2.X[2];
  T[2][2] = P3.X[2];
}

tensor3d::
tensor3d(const tensor3d &v) : id(v.id)
{
  T[0][0] = v.T[0][0];
  T[0][1] = v.T[0][1];
  T[0][2] = v.T[0][2];
  
  T[1][0] = v.T[1][0];
  T[1][1] = v.T[1][1];
  T[1][2] = v.T[1][2];
  
  T[2][0] = v.T[2][0];
  T[2][1] = v.T[2][1];
  T[2][2] = v.T[2][2];
}

tensor3d::
~tensor3d()
{ }



//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
tensor3d::
setIJ(const UInt & i, const UInt & j, const Real & value)
{
  assert(i>=1); assert(j>=1);
  assert(i<=3); assert(j<=3);
  
  T[i-1][j-1] = value;
}

void
tensor3d::
sumIJ(const UInt & i, const UInt & j, const Real & value)
{
  assert(i>=1); assert(j>=1);
  assert(i<=3); assert(j<=3);
  
  T[i-1][j-1] += value;
}

void
tensor3d::
setRow(const UInt & i, const point3d & value)
{
  assert(i >= 1);
  assert(i <= 3);
  
  T[i-1][0] = value.getX();
  T[i-1][1] = value.getY();
  T[i-1][2] = value.getZ();
}
      
void
tensor3d::
setCol(const UInt & j, const point3d & value)
{
  assert(j >= 1);
  assert(j <= 3);
  
  T[0][j-1] = value.getX();
  T[1][j-1] = value.getY();
  T[2][j-1] = value.getZ();
}
    
void
tensor3d::
sumRow(const UInt & i, const point3d & value)
{
  assert(i >= 1);
  assert(i <= 3);
  
  T[i-1][0] += value.getX();
  T[i-1][1] += value.getY();
  T[i-1][2] += value.getZ();
}
    
void
tensor3d::
sumCol(const UInt & j, const point3d & value)
{
  assert(j >= 1);
  assert(j <= 3);
  
  T[0][j-1] += value.getX();
  T[1][j-1] += value.getY();
  T[2][j-1] += value.getZ();
}
      
void
tensor3d::
setId(const UInt & Id)
{
  id = Id;
}



//_________________________________________________________________________________________________
// GET FUNCTIONS
//------------------------------------------------------------------------------------------------- 
point3d 
tensor3d::
getRow(const UInt & i) const
{
  assert(i >= 1);
  assert(i <= 3);
 
  return(point3d(T[i-1][0], T[i-1][1], T[i-1][2]));
}
        
point3d 
tensor3d::
getCol(const UInt & j) const
{
  assert(j >= 1);
  assert(j <= 3);
  
  return(point3d(T[0][j-1], T[1][j-1], T[2][j-1]));
}
         
UInt
tensor3d::
getId() const
{
  return(id);
}



//_________________________________________________________________________________________________
// OPERATORS
//-------------------------------------------------------------------------------------------------
Real &
tensor3d::
operator()(const UInt & i, const UInt & j)
{
  assert((i>=1) & (j>=1));
  assert((i<=3) & (j<=3));
  
  return(T[i-1][j-1]);
}
    
const Real &
tensor3d::
operator()(const UInt & i, const UInt & j) const
{
  assert((i>=1) & (j>=1));
  assert((i<=3) & (j<=3));
  
  return(T[i-1][j-1]);
}
    
tensor3d &
tensor3d::
operator=(const tensor3d & v)
{
  T[0][0] = v.T[0][0];
  T[0][1] = v.T[0][1];
  T[0][2] = v.T[0][2];
  
  T[1][0] = v.T[1][0];
  T[1][1] = v.T[1][1];
  T[1][2] = v.T[1][2];
  
  T[2][0] = v.T[2][0];
  T[2][1] = v.T[2][1];
  T[2][2] = v.T[2][2];
  
  id = v.id;
  
  return *this;
}
    
void
tensor3d::
operator+=(const tensor3d & V)
{
  T[0][0] += V.T[0][0];
  T[0][1] += V.T[0][1];
  T[0][2] += V.T[0][2];
  
  T[1][0] += V.T[1][0];
  T[1][1] += V.T[1][1];
  T[1][2] += V.T[1][2];
  
  T[2][0] += V.T[2][0];
  T[2][1] += V.T[2][1];
  T[2][2] += V.T[2][2];
}

void
tensor3d::
operator-=(const tensor3d & V)
{
  T[0][0] += V.T[0][0];
  T[0][1] += V.T[0][1];
  T[0][2] += V.T[0][2];
  
  T[1][0] += V.T[1][0];
  T[1][1] += V.T[1][1];
  T[1][2] += V.T[1][2];
  
  T[2][0] += V.T[2][0];
  T[2][1] += V.T[2][1];
  T[2][2] += V.T[2][2];
}
    
void
tensor3d::
operator*=(const Real & a)
{
  T[0][0] *= a;
  T[0][1] *= a;
  T[0][2] *= a;
  
  T[1][0] *= a;
  T[1][1] *= a;
  T[1][2] *= a;
  
  T[2][0] *= a;
  T[2][1] *= a;
  T[2][2] *= a;
}
  
void
tensor3d::
operator/=(const Real & a)
{
  T[0][0] /= a;
  T[0][1] /= a;
  T[0][2] /= a;
  
  T[1][0] /= a;
  T[1][1] /= a;
  T[1][2] /= a;
  
  T[2][0] /= a;
  T[2][1] /= a;
  T[2][2] /= a;
}



//_________________________________________________________________________________________________
// INTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
tensor3d::
completeSecondColoumn()
{
  point3d X(1.0,0.0,0.0);
  point3d Y(0.0,1.0,0.0);
  
  point3d N(this->getCol(1));
  N /= N.norm2();
  
  point3d Vx = X.combinationGS(N);
  point3d Vy = Y.combinationGS(N);
  
  if(Vx.norm2() > Vy.norm2())
  {
    Vx /= Vx.norm2();
    assert(Vx.norm2() != 0.0);
    
    this->setCol(2,Vx);
  }
  else
  {
    Vy /= Vy.norm2();
    assert(Vy.norm2() != 0.0);
    
    this->setCol(2,Vy);
  }
}

void
tensor3d::
completeSecondRow()
{
  point3d X(1.0,0.0,0.0);
  point3d Y(0.0,1.0,0.0);
  
  point3d N(this->getRow(1));
  N /= N.norm2();
  
  point3d Vx = X.combinationGS(N);
  point3d Vy = Y.combinationGS(N);
  
  if(Vx.norm2() > Vy.norm2())
  {
    Vx /= Vx.norm2();
    assert(Vx.norm2() != 0.0);
    
    this->setRow(2,Vx);
  }
  else
  {
    Vy /= Vy.norm2();
    assert(Vy.norm2() != 0.0);
    
    this->setRow(2,Vy);
  }
}

void
tensor3d::
completeThirdColoumn()
{
  point3d V1 = this->getCol(1);
  point3d V2 = this->getCol(2);
  point3d V3 = V1^V2;
  
  V3 /= V3.norm2();
  assert(V3.norm2() != 0.0);
  
  this->setCol(3,V3);
}

void
tensor3d::
completeThirdRow()
{
  point3d V1 = this->getRow(1);
  point3d V2 = this->getRow(2);
  point3d V3 = V1^V2;
  
  V3 /= V3.norm2();
  assert(V3.norm2() != 0.0);
  
  this->setRow(3,V3);
}



//_________________________________________________________________________________________________
// EXTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
void 
tensor3d::
computeInverse()
{
  Real determinante;
  Real inv[3][3];
  
  determinante =  T[0][0] * T[1][1] * T[2][2] - T[0][0] * T[1][2] * T[2][1] - T[1][0] * T[0][1] * T[2][2] 
  + T[1][0] * T[0][2] * T[2][1] + T[2][0] * T[0][1] * T[1][2] - T[2][0] * T[0][2] * T[1][1];
  
  inv[0][0] =  (T[1][1] * T[2][2] - T[1][2] * T[2][1]) / determinante;
  inv[0][1] = -(T[0][1] * T[2][2] - T[0][2] * T[2][1]) / determinante;
  inv[0][2] =  (T[0][1] * T[1][2] - T[0][2] * T[1][1]) / determinante;
  inv[1][0] = -(T[1][0] * T[2][2] - T[1][2] * T[2][0]) / determinante;
  inv[1][1] =  (T[0][0] * T[2][2] - T[0][2] * T[2][0]) / determinante;
  inv[1][2] = -(T[0][0] * T[1][2] - T[0][2] * T[1][0]) / determinante;
  inv[2][0] =  (T[1][0] * T[2][1] - T[1][1] * T[2][0]) / determinante;
  inv[2][1] = -(T[0][0] * T[2][1] - T[0][1] * T[2][0]) / determinante;
  inv[2][2] =  (T[0][0] * T[1][1] - T[0][1] * T[1][0]) / determinante;
  
  
  T[0][0] = inv[0][0];
  T[0][1] = inv[0][1];
  T[0][2] = inv[0][2];
  
  T[1][0] = inv[1][0];
  T[1][1] = inv[1][1];
  T[1][2] = inv[1][2];
  
  T[2][0] = inv[2][0];
  T[2][1] = inv[2][1];
  T[2][2] = inv[2][2];
}

void 
tensor3d::
transpose()
{
  Real transp[3][3];
  
  //Copia
  for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      transp[i][j] = T[i][j];
    }
  }
  
  //Scrittura trasposta
  for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      T[i][j] = transp[j][i];
    }
  }
}



//_________________________________________________________________________________________________
// EXTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
double 
tensor3d::
getFirstInvariant() const
{
  return(T[0][0] + T[1][1] + T[2][2]);
}

double
tensor3d::
getSecondInvariant() const
{
  return( T[0][0]*(T[1][1] + T[2][2]) + T[1][1]*T[2][2] - (T[0][1]*T[1][0] + T[0][2]*T[2][0] + T[1][2]*T[2][1]) );
}

double
tensor3d::
getThirdInvariant() const
{
  return(
    T[0][0] * T[1][1] * T[2][2] - T[0][0] * T[1][2] * T[2][1] 
  - T[1][0] * T[0][1] * T[2][2] + T[1][0] * T[0][2] * T[2][1] 
  + T[2][0] * T[0][1] * T[1][2] - T[2][0] * T[0][2] * T[1][1] );
}

double
tensor3d::
getFrobenius() const
{
  Real tot = 0.0;
  
  for(UInt i=0;i<3;++i)
  {
    for(UInt j=0;j<3;++j)
    {
      tot += pow(T[i][j],2.0);
    }
  }
  
  return(pow(tot,0.5));
}

point3d
tensor3d::
firstIndexSaturation(const point3d & v) const
{
  Real X = T[0][0] * v.getX() + T[1][0] * v.getY() + T[2][0] * v.getZ();
  Real Y = T[0][1] * v.getX() + T[1][1] * v.getY() + T[2][1] * v.getZ();
  Real Z = T[0][2] * v.getX() + T[1][2] * v.getY() + T[2][2] * v.getZ();
  
  return(point3d(X,Y,Z));
}

point3d
tensor3d::
secondIndexSaturation(const point3d & v) const
{
  Real X = T[0][0] * v.getX() + T[0][1] * v.getY() + T[0][2] * v.getZ();
  Real Y = T[1][0] * v.getX() + T[1][1] * v.getY() + T[1][2] * v.getZ();
  Real Z = T[2][0] * v.getX() + T[2][1] * v.getY() + T[2][2] * v.getZ();
  
  return(point3d(X,Y,Z));
}

void 
tensor3d::
directProduct(const point3d & a, const point3d & b)
{
  T[0][0] = a.getX() * b.getX();
  T[0][1] = a.getX() * b.getY();
  T[0][2] = a.getX() * b.getZ();
  
  T[1][0] = a.getY() * b.getX();
  T[1][1] = a.getY() * b.getY();
  T[1][2] = a.getY() * b.getZ();
  
  T[2][0] = a.getZ() * b.getX();
  T[2][1] = a.getZ() * b.getY();
  T[2][2] = a.getZ() * b.getZ();
}

tensor3d
tensor3d::
matrixProduct(const tensor3d & A) const
{
  tensor3d out;
  
  for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      for(UInt k=0; k<3; ++k)
      {
	out.sumIJ(i+1, j+1, T[i][k] * A.T[k][j] );
      }
    }
  }
  
  return out;
}

Real
tensor3d::
scalarProduct(const tensor3d & A) const
{
  Real out = 0.0;
  
  for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      out += T[i][j] * A.T[i][j];
    }
  }
  
  return(out);
}

void
tensor3d::
ricciThirdSaturation(const point3d & v)
{
  T[0][1] =  v.getZ();
  T[0][2] = -v.getY();
  T[1][0] = -v.getZ();
  T[1][2] =  v.getX();
  T[2][0] =  v.getY();
  T[2][1] = -v.getX();
}

point3d
tensor3d::
antiSymmetric(const tensor3d & A) const
{
  return(point3d(
  A.getIJ(2,3)-A.getIJ(3,2),
  A.getIJ(3,1)-A.getIJ(1,3),
  A.getIJ(1,2)-A.getIJ(2,1)
  ));
}



//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool
tensor3d::
operator<(const tensor3d & V) const
{
  /*for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      if( T[i][j] < (V.T[i][j] - geoToll) )
      {
	return(true);
      }
      
      if( T[i][j] > (V.T[i][j] + geoToll) )
      {
	return(false);
      }
    }
  }
  
  return(false);*/
  
  for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      if( (geoToll * round(T[i][j] / geoToll)) < (geoToll * round(V.T[i][j] / geoToll)) ) { return(true); }
      if( (geoToll * round(T[i][j] / geoToll)) > (geoToll * round(V.T[i][j] / geoToll)) ) { return(false); }
    }
  }
  
  return(false);
}

bool 
tensor3d::
operator!=(const tensor3d & V) const
{
  /*for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      if( ( T[i][j] > (V.T[i][j] + geoToll) ) || ( T[i][j] < (V.T[i][j] - geoToll) ) )
      {
	return(true);
      }
    }
  }
  
  return(false);*/
  
  for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      if( (geoToll * round(T[i][j] / geoToll)) != (geoToll * round(V.T[i][j] / geoToll)) ) { return(true); }
    }
  }
  
  return(false);
}



//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
ostream &
operator<<(ostream & f, const tensor3d & A )
{
  f << "Tensor components" << std::endl;
  
  for(UInt i=0; i<3; ++i)
  {
    for(UInt j=0; j<3; ++j)
    {
      f << A.T[i][j] << " ";
    }
    
    f << std::endl;
  }
  
  f << endl;
  
  return f;
}

