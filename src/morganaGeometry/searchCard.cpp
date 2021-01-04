#include "searchCard.h"

searchCard::
searchCard()
{
}

void
searchCard::
setPmin(const point3d & P)
{
  Pmin = P;
}

void
searchCard::
setPmax(const point3d & P)
{
  Pmax = P;
}

void
searchCard::
setConnected(const UInt & num)
{
  numConnected = num;
}

const point3d &
searchCard::
getPmin() const
{
  return(Pmin);
}

const point3d &
searchCard::
getPmax() const
{
  return(Pmax);
}

point3d &
searchCard::
getPmin()
{
  return(Pmin);
}

point3d &
searchCard::
getPmax()
{
  return(Pmax);
}

const UInt &
searchCard::
getConnected() const
{
  return(numConnected);
}

ostream &
operator<<(ostream & f, const searchCard & C)
{
  return f << "Pmin : " << C.Pmin << "Pmax :" << C.Pmax << "NumElements :" << C.numConnected;
}