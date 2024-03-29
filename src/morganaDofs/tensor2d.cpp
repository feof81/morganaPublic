#include "tensor2d.h"

//_________________________________________________________________________________________________
// CONSTRUCTORS
//-------------------------------------------------------------------------------------------------
tensor2d::
tensor2d() : id(1)
{
  T[0][0] = 0.0;
  T[0][1] = 0.0;
  
  T[1][0] = 0.0;
  T[1][1] = 0.0;
}

tensor2d::
tensor2d(const tensor2d &v) : id(v.id)
{
  T[0][0] = v.T[0][0];
  T[0][1] = v.T[0][1];
  
  T[1][0] = v.T[1][0];
  T[1][1] = v.T[1][1];
}

tensor2d::
~tensor2d()
{ }



//_________________________________________________________________________________________________
// SET FUNCTIONS
//-------------------------------------------------------------------------------------------------
void
tensor2d::
setIJ(const UInt & i, const UInt & j, const Real & value)
{
  assert(i>=1); assert(j>=1);
  assert(i<=2); assert(j<=2);
  
  T[i-1][j-1] = value;
}

void
tensor2d::
sumIJ(const UInt & i, const UInt & j, const Real & value)
{
  assert(i>=1); assert(j>=1);
  assert(i<=2); assert(j<=2);
  
  T[i-1][j-1] += value;
}

void
tensor2d::
setRow(const UInt & i, const point2d & value)
{
  assert(i >= 1);
  assert(i <= 2);
  
  T[i-1][0] = value.getX();
  T[i-1][1] = value.getY();
}
      
void
tensor2d::
setCol(const UInt & j, const point2d & value)
{
  assert(j >= 1);
  assert(j <= 2);
  
  T[0][j-1] = value.getX();
  T[1][j-1] = value.getY();
}
    
void
tensor2d::
sumRow(const UInt & i, const point2d & value)
{
  assert(i >= 1);
  assert(i <= 2);
  
  T[i-1][0] += value.getX();
  T[i-1][1] += value.getY();
}
    
void
tensor2d::
sumCol(const UInt & j, const point2d & value)
{
  assert(j >= 1);
  assert(j <= 2);
  
  T[0][j-1] += value.getX();
  T[1][j-1] += value.getY();
}
      
void
tensor2d::
setId(const UInt & Id)
{
  id = Id;
}



//_________________________________________________________________________________________________
// GET FUNCTIONS
//------------------------------------------------------------------------------------------------- 
point2d 
tensor2d::
getRow(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 2);
 
  return(point2d(T[i-1][0], T[i-1][1]));
}
        
point2d 
tensor2d::
getCol(const UInt & j)
{
  assert(j >= 1);
  assert(j <= 2);
  
  return(point2d(T[0][j-1], T[1][j-1]));
}
         
UInt
tensor2d::
getId() const
{
  return(id);
}


//_________________________________________________________________________________________________
// OPERATORS
//-------------------------------------------------------------------------------------------------
Real
tensor2d::
operator()(const UInt & i, const UInt & j)
{
  assert((i>=1) & (j>=1));
  assert((i<=2) & (j<=2));
  
  return(T[i-1][j-1]);
}
    
const Real &
tensor2d::
operator()(const UInt & i, const UInt & j) const
{
  assert((i>=1) & (j>=1));
  assert((i<=2) & (j<=2));
  
  return(T[i-1][j-1]);
}
    
tensor2d &
tensor2d::
operator=(const tensor2d & v)
{
  T[0][0] = v.T[0][0];
  T[0][1] = v.T[0][1];
  
  T[1][0] = v.T[1][0];
  T[1][1] = v.T[1][1];
  
  id = v.id;
  
  return *this;
}
    
void
tensor2d::
operator+=(const tensor2d & V)
{
  T[0][0] += V.T[0][0];
  T[0][1] += V.T[0][1];
  
  T[1][0] += V.T[1][0];
  T[1][1] += V.T[1][1];
}

void
tensor2d::
operator-=(const tensor2d & V)
{
  T[0][0] += V.T[0][0];
  T[0][1] += V.T[0][1];
  
  T[1][0] += V.T[1][0];
  T[1][1] += V.T[1][1];
}
    
void
tensor2d::
operator*=(const Real & a)
{
  T[0][0] *= a;
  T[0][1] *= a;
  
  T[1][0] *= a;
  T[1][1] *= a;
}
  
void
tensor2d::
operator/=(const Real & a)
{
  T[0][0] /= a;
  T[0][1] /= a;
  
  T[1][0] /= a;
  T[1][1] /= a;
}



//_________________________________________________________________________________________________
// EXTERNAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
void 
tensor2d::
computeInverse()
{
  Real determinante;
  Real inv[2][2];
  
  determinante = T[0][0] * T[1][1] - T[0][1] * T[1][0];
  
  inv[0][0] =  T[1][1] / determinante;
  inv[0][1] = -T[0][1] / determinante;
  inv[1][0] = -T[1][0] / determinante;
  inv[1][1] =  T[0][0] / determinante;
  
  
  T[0][0] = inv[0][0];
  T[0][1] = inv[0][1];
  
  T[1][0] = inv[1][0];
  T[1][1] = inv[1][1];
}

void 
tensor2d::
transpose()
{
  Real transp[2][2];
  
  //Copia
  for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    {
      transp[i][j] = T[i][j];
    }
  }
  
  //Scrittura trasposta
  for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    {
      T[i][j] = transp[j][i];
    }
  }
}



//_________________________________________________________________________________________________
// EXTERNAL FUNCTIONS 
//-------------------------------------------------------------------------------------------------
Real
tensor2d::
getFrobenius() const
{
  Real tot = 0.0;
  
  for(UInt i=0;i<2;++i)
  {
    for(UInt j=0;j<2;++j)
    {
      tot += pow(T[i][j],2.0);
    }
  }
  
  return(pow(tot,0.5));
}
      
point2d
tensor2d::
firstIndexSaturation(const point2d & v) const
{
  Real X = T[0][0] * v.getX() + T[1][0] * v.getY();
  Real Y = T[0][1] * v.getX() + T[1][1] * v.getY();
  
  return(point2d(X,Y));
}
      
point2d
tensor2d::
secondIndexSaturation(const point2d & v) const
{
  Real X = T[0][0] * v.getX() + T[0][1] * v.getY();
  Real Y = T[1][0] * v.getX() + T[1][1] * v.getY();
  
  return(point2d(X,Y));
}
      
void
tensor2d::
directProduct(const point2d & a, const point2d & b)
{
  T[0][0] = a.getX() * b.getX();
  T[0][1] = a.getX() * b.getY();
  
  T[1][0] = a.getY() * b.getX();
  T[1][1] = a.getY() * b.getY();
}
      
tensor2d
tensor2d::
matrixProduct(const tensor2d & A) const
{
  tensor2d out;
  
  for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    {
      for(UInt k=0; k<2; ++k)
      {
	out.sumIJ(i+1, j+1, T[i][k] * A.T[k][j] );
      }
    }
  }
  
  return out;
}
      
Real
tensor2d::
scalarProduct(const tensor2d & A) const
{
  Real out = 0.0;
  
  for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    {
      out += T[i][j] * A.T[i][j];
    }
  }
  
  return(out);
}



//_________________________________________________________________________________________________
// ORDINAL FUNCTIONS
//-------------------------------------------------------------------------------------------------
bool
tensor2d::
operator<(const tensor2d & V) const
{
  /*for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
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
  
  for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    {
      if( (geoToll * round(T[i][j] / geoToll)) < (geoToll * round(V.T[i][j] / geoToll)) ) { return(true); }
      if( (geoToll * round(T[i][j] / geoToll)) > (geoToll * round(V.T[i][j] / geoToll)) ) { return(false); }
    }
  }
  
  return(false);
}

bool 
tensor2d::
operator!=(const tensor2d & V) const
{
  /*for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    {
      if( ( T[i][j] > (V.T[i][j] + geoToll) ) || ( T[i][j] < (V.T[i][j] - geoToll) ) )
      {
	return(true);
      }
    }
  }
  
  return(false);*/
  
  for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    { if( (geoToll * round(T[i][j] / geoToll)) != (geoToll * round(V.T[i][j] / geoToll)) ) { return(true); } }
  }
  
  return(false);
}


//_________________________________________________________________________________________________
// OTHER FUNCTIONS
//-------------------------------------------------------------------------------------------------
ostream &
operator<<(ostream & f, const tensor2d & A )
{
  f << "Tensor components" << std::endl;
  
  for(UInt i=0; i<2; ++i)
  {
    for(UInt j=0; j<2; ++j)
    {
      f << A.T[i][j] << " ";
    }
    
    f << std::endl;
  }
  
  f << endl;
  
  return f;
}

size_t
tensor2d::
memSize() const
{
  return(sizeof(tensor2d));
}
