/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef POLYCARDS_H
#define POLYCARDS_H

#include "typesInterface.hpp"

typedef const UInt SInt;
typedef const Real SReal;




/*! -------------------------------------------------------------------------------------------- */
/*!                                         POLYCARDS 1D                                         */
/*! -------------------------------------------------------------------------------------------- */

/*! Polynomial P0 - 1D */
template<Int N> struct trunk_P0_1d_A;
template<>      struct trunk_P0_1d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct P0_1d_A {static const UInt S = 1; typedef trunk_P0_1d_A<1> ROOT;};


/*! Polynomial P1 - 1D */
template<Int N> struct trunk_P1_1d_A;
template<>      struct trunk_P1_1d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1_1d_A<2> {typedef trunk_P1_1d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct P1_1d_A {static const UInt S = 2; typedef trunk_P1_1d_A<2> ROOT; };

template<Int N> struct trunk_P1_1d_B;
template<>      struct trunk_P1_1d_B<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct P1_1d_B {static const UInt S = 1; typedef trunk_P1_1d_B<1> ROOT; };


/*! Polynomial P2 - 1D */
template<Int N> struct trunk_P2_1d_A;
template<>      struct trunk_P2_1d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_1d_A<2> {typedef trunk_P2_1d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_1d_A<3> {typedef trunk_P2_1d_A<2> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct P2_1d_A {static const UInt S = 3; typedef trunk_P2_1d_A<3> ROOT; };

template<Int N> struct trunk_P2_1d_B;
template<>      struct trunk_P2_1d_B<1> {                              static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_1d_B<2> {typedef trunk_P2_1d_B<1> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct P2_1d_B {static const UInt S = 2; typedef trunk_P2_1d_B<2> ROOT; };

template<Int N> struct trunk_P2_1d_C;
template<>      struct trunk_P2_1d_C<1> {                              static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_1d_C<2> {typedef trunk_P2_1d_C<1> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct P2_1d_C {static const UInt S = 2; typedef trunk_P2_1d_C<2> ROOT; };


/*! Polynomial D0 - 1D */
template<Int N> struct trunk_D0_1d_A;
template<>      struct trunk_D0_1d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct D0_1d_A {static const UInt S = 1; typedef trunk_D0_1d_A<1> ROOT; };


/*! Polynomial DG0 - 1D */
template<Int N> struct trunk_DG0_1d_A;
template<>      struct trunk_DG0_1d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct DG0_1d_A {static const UInt S = 1; typedef trunk_DG0_1d_A<1> ROOT; };


/*! Polynomial DG1 - 1D */
template<Int N> struct trunk_DG1_1d_A;
template<>      struct trunk_DG1_1d_A<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct DG1_1d_A {static const UInt S = 1; typedef trunk_DG1_1d_A<1> ROOT; };

template<Int N> struct trunk_DG1_1d_B;
template<>      struct trunk_DG1_1d_B<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_DG1_1d_B<2> {typedef trunk_DG1_1d_B<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct DG1_1d_B {static const UInt S = 2; typedef trunk_DG1_1d_B<2> ROOT; };




/*! -------------------------------------------------------------------------------------------- */
/*!                                         POLYCARDS 2D                                         */
/*! -------------------------------------------------------------------------------------------- */

/*! Polynomial P0 - 2D */
template<Int N> struct trunk_P0_2d_A;
template<>      struct trunk_P0_2d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P0_2d_A    {static const UInt S = 1; typedef trunk_P0_2d_A<1> ROOT;};      


/*! Polynomial P1 - 2D */
template<Int N> struct trunk_P1_2d_A;
template<>      struct trunk_P1_2d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1_2d_A<2> {typedef trunk_P1_2d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1_2d_A<3> {typedef trunk_P1_2d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P1_2d_A    {static const UInt S = 3; typedef trunk_P1_2d_A<3> ROOT;};

template<Int N> struct trunk_P1_2d_B;
template<>      struct trunk_P1_2d_B<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P1_2d_B    {static const UInt S = 1; typedef trunk_P1_2d_B<1> ROOT;};

template<Int N> struct trunk_P1_2d_C;
template<>      struct trunk_P1_2d_C<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P1_2d_C    {static const UInt S = 1; typedef trunk_P1_2d_C<1> ROOT;};


/*! Polynomial P2 - 2D */
template<Int N> struct trunk_P2_2d_A;
template<>      struct trunk_P2_2d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_A<2> {typedef trunk_P2_2d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_A<3> {typedef trunk_P2_2d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_A<4> {typedef trunk_P2_2d_A<3> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_A<5> {typedef trunk_P2_2d_A<4> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_A<6> {typedef trunk_P2_2d_A<5> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P2_2d_A    {static const UInt S = 6; typedef trunk_P2_2d_A<6> ROOT;};

template<Int N> struct trunk_P2_2d_B;
template<>      struct trunk_P2_2d_B<1> {                              static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_B<2> {typedef trunk_P2_2d_B<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P2_2d_B    {static const UInt S = 2; typedef trunk_P2_2d_B<2> ROOT;};

template<Int N> struct trunk_P2_2d_C;
template<>      struct trunk_P2_2d_C<1> {                              static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_C<2> {typedef trunk_P2_2d_C<1> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P2_2d_C    {static const UInt S = 2; typedef trunk_P2_2d_C<2> ROOT;};      

template<Int N> struct trunk_P2_2d_D;
template<>      struct trunk_P2_2d_D<1> {                              static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_D<2> {typedef trunk_P2_2d_D<1> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_D<3> {typedef trunk_P2_2d_D<2> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P2_2d_D    {static const UInt S = 3; typedef trunk_P2_2d_D<3> ROOT;};       

template<Int N> struct trunk_P2_2d_E;
template<>      struct trunk_P2_2d_E<1> {static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P2_2d_E    {static const UInt S = 1; typedef trunk_P2_2d_E<1> ROOT;};      

template<Int N> struct trunk_P2_2d_F;
template<>      struct trunk_P2_2d_F<1> {                              static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_F<2> {typedef trunk_P2_2d_F<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_2d_F<3> {typedef trunk_P2_2d_F<2> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
struct                       P2_2d_F    {static const UInt S = 3; typedef trunk_P2_2d_F<3> ROOT;};        


/*! Polynomial Q0 - 2D */
template<Int N> struct trunk_Q0_2d_A;
template<>      struct trunk_Q0_2d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       Q0_2d_A    {static const UInt S = 1; typedef trunk_Q0_2d_A<1> ROOT;};


/*! Polynomial Q1 - 2D */
template<Int N> struct trunk_Q1_2d_A;
template<>      struct trunk_Q1_2d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_2d_A<2> {typedef trunk_Q1_2d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_2d_A<3> {typedef trunk_Q1_2d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_2d_A<4> {typedef trunk_Q1_2d_A<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q1_2d_A    {static const UInt S = 4; typedef trunk_Q1_2d_A<4> ROOT;};

template<Int N> struct trunk_Q1_2d_B;
template<>      struct trunk_Q1_2d_B<1> {                               static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_2d_B<2> {typedef  trunk_Q1_2d_B<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q1_2d_B    {static const UInt S = 2; typedef trunk_Q1_2d_B<2> ROOT;};

template<Int N> struct trunk_Q1_2d_C;
template<>      struct trunk_Q1_2d_C<1> {static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q1_2d_C    {static const UInt S = 1; typedef trunk_Q1_2d_C<1> ROOT;};

template<Int N> struct trunk_Q1_2d_D;
template<>      struct trunk_Q1_2d_D<1> {                               static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_2d_D<2> {typedef  trunk_Q1_2d_D<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q1_2d_D    {static const UInt S = 2; typedef trunk_Q1_2d_D<2> ROOT;};


/*! Polynomial Q2 - 2D */
template<Int N> struct trunk_Q2_2d_A;
template<>      struct trunk_Q2_2d_A<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<2> {typedef trunk_Q2_2d_A<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<3> {typedef trunk_Q2_2d_A<2> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<4> {typedef trunk_Q2_2d_A<3> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<5> {typedef trunk_Q2_2d_A<4> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<6> {typedef trunk_Q2_2d_A<5> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<7> {typedef trunk_Q2_2d_A<6> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<8> {typedef trunk_Q2_2d_A<7> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_A<9> {typedef trunk_Q2_2d_A<8> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_A    {static const UInt S = 9; typedef trunk_Q2_2d_A<9> ROOT;};

template<Int N> struct trunk_Q2_2d_B;
template<>      struct trunk_Q2_2d_B<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_B<2> {typedef trunk_Q2_2d_B<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_B<3> {typedef trunk_Q2_2d_B<2> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_B<4> {typedef trunk_Q2_2d_B<3> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_B<5> {typedef trunk_Q2_2d_B<4> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_B<6> {typedef trunk_Q2_2d_B<5> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_B    {static const UInt S = 6; typedef trunk_Q2_2d_B<6> ROOT;};

template<Int N> struct trunk_Q2_2d_C;
template<>      struct trunk_Q2_2d_C<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_C<2> {typedef trunk_Q2_2d_C<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_C<3> {typedef trunk_Q2_2d_C<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_C<4> {typedef trunk_Q2_2d_C<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_C    {static const UInt S = 4; typedef trunk_Q2_2d_C<4> ROOT;};

template<Int N> struct trunk_Q2_2d_D;
template<>      struct trunk_Q2_2d_D<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_D<2> {typedef trunk_Q2_2d_D<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_D<3> {typedef trunk_Q2_2d_D<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_D<4> {typedef trunk_Q2_2d_D<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_D<5> {typedef trunk_Q2_2d_D<4> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_D<6> {typedef trunk_Q2_2d_D<5> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_D    {static const UInt S = 6; typedef trunk_Q2_2d_D<6> ROOT;};

template<Int N> struct trunk_Q2_2d_E;
template<>      struct trunk_Q2_2d_E<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_E<2> {typedef trunk_Q2_2d_E<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_E<3> {typedef trunk_Q2_2d_E<2> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_E<4> {typedef trunk_Q2_2d_E<3> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_E<5> {typedef trunk_Q2_2d_E<4> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_E<6> {typedef trunk_Q2_2d_E<5> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_E    {static const UInt S = 6; typedef trunk_Q2_2d_E<6> ROOT;};

template<Int N> struct trunk_Q2_2d_F;
template<>      struct trunk_Q2_2d_F<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_F<2> {typedef trunk_Q2_2d_F<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_F<3> {typedef trunk_Q2_2d_F<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_F<4> {typedef trunk_Q2_2d_F<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_F    {static const UInt S = 4; typedef trunk_Q2_2d_F<4> ROOT;};

template<Int N> struct trunk_Q2_2d_G;
template<>      struct trunk_Q2_2d_G<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_G<2> {typedef trunk_Q2_2d_G<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_G<3> {typedef trunk_Q2_2d_G<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_G<4> {typedef trunk_Q2_2d_G<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_G    {static const UInt S = 4; typedef trunk_Q2_2d_G<4> ROOT;};

template<Int N> struct trunk_Q2_2d_H;
template<>      struct trunk_Q2_2d_H<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_H<2> {typedef trunk_Q2_2d_H<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_H<3> {typedef trunk_Q2_2d_H<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_H<4> {typedef trunk_Q2_2d_H<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_H<5> {typedef trunk_Q2_2d_H<4> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_H<6> {typedef trunk_Q2_2d_H<5> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_H    {static const UInt S = 6; typedef trunk_Q2_2d_H<6> ROOT;};

template<Int N> struct trunk_Q2_2d_I;
template<>      struct trunk_Q2_2d_I<1> {                              static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_I<2> {typedef trunk_Q2_2d_I<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_I<3> {typedef trunk_Q2_2d_I<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_2d_I<4> {typedef trunk_Q2_2d_I<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       Q2_2d_I    {static const UInt S = 4; typedef trunk_Q2_2d_I<4> ROOT;};


/*! Polynomial D0 - 2D */
template<Int N> struct trunk_D0_2d_A;
template<>      struct trunk_D0_2d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       D0_2d_A    {static const UInt S = 1; typedef trunk_D0_2d_A<1> ROOT;};      

template<Int N> struct trunk_D0_2d_B;
template<>      struct trunk_D0_2d_B<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       D0_2d_B    {static const UInt S = 1; typedef trunk_D0_2d_B<1> ROOT;}; 
      
template<Int N> struct trunk_D0_2d_C;
template<>      struct trunk_D0_2d_C<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       D0_2d_C    {static const UInt S = 1; typedef trunk_D0_2d_C<1> ROOT;};       
      

/*! Polynomial RT0 - 2D */
template<Int N> struct trunk_RT0_2d_Ax;
template<>      struct trunk_RT0_2d_Ax<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Ax    {static const UInt S = 1; typedef trunk_RT0_2d_Ax<1> ROOT;}; 

template<Int N> struct trunk_RT0_2d_Ay;
template<>      struct trunk_RT0_2d_Ay<1> {                                 static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_RT0_2d_Ay<2> {typedef  trunk_RT0_2d_Ay<1> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Ay    {static const UInt S = 2; typedef trunk_RT0_2d_Ay<2> ROOT;}; 

template<Int N> struct trunk_RT0_2d_Az;
template<>      struct trunk_RT0_2d_Az<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Az    {static const UInt S = 1; typedef trunk_RT0_2d_Az<1> ROOT;}; 

template<Int N> struct trunk_RT0_2d_Bx;
template<>      struct trunk_RT0_2d_Bx<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Bx    {static const UInt S = 1; typedef trunk_RT0_2d_Bx<1> ROOT;};

template<Int N> struct trunk_RT0_2d_By;
template<>      struct trunk_RT0_2d_By<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_By    {static const UInt S = 1; typedef trunk_RT0_2d_By<1> ROOT;};

template<Int N> struct trunk_RT0_2d_Bz;
template<>      struct trunk_RT0_2d_Bz<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Bz    {static const UInt S = 1; typedef trunk_RT0_2d_Bz<1> ROOT;}; 

template<Int N> struct trunk_RT0_2d_Cx;
template<>      struct trunk_RT0_2d_Cx<1> {                                static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_RT0_2d_Cx<2> {typedef trunk_RT0_2d_Cx<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Cx    {static const UInt S = 2; typedef trunk_RT0_2d_Cx<2> ROOT;};
      
template<Int N> struct trunk_RT0_2d_Cy;
template<>      struct trunk_RT0_2d_Cy<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Cy    {static const UInt S = 1; typedef trunk_RT0_2d_Cy<1> ROOT;};    

template<Int N> struct trunk_RT0_2d_Cz;
template<>      struct trunk_RT0_2d_Cz<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_2d_Cz    {static const UInt S = 1; typedef trunk_RT0_2d_Cz<1> ROOT;};       


/*! Polynomial DG0 - 2D */
template<Int N> struct trunk_DG0_2d_A;
template<>      struct trunk_DG0_2d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       DG0_2d_A    {static const UInt S = 1; typedef trunk_DG0_2d_A<1> ROOT;};


/*! Polynomial DG1 - 2D */
template<Int N> struct trunk_DG1_2d_A;
template<>      struct trunk_DG1_2d_A<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_DG1_2d_A<2> {typedef trunk_DG1_2d_A<1> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       DG1_2d_A    {static const UInt S = 2; typedef trunk_DG1_2d_A<2> ROOT;};

template<Int N> struct trunk_DG1_2d_B;
template<>      struct trunk_DG1_2d_B<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_DG1_2d_B<2> {typedef trunk_DG1_2d_B<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_DG1_2d_B<3> {typedef trunk_DG1_2d_B<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       DG1_2d_B    {static const UInt S = 3; typedef trunk_DG1_2d_B<3> ROOT;};

template<Int N> struct trunk_DG1_2d_C;
template<>      struct trunk_DG1_2d_C<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_DG1_2d_C<2> {typedef trunk_DG1_2d_C<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       DG1_2d_C    {static const UInt S = 2; typedef trunk_DG1_2d_C<2> ROOT;};




/*! -------------------------------------------------------------------------------------------- */
/*!                                         POLYCARDS 3D                                         */
/*! -------------------------------------------------------------------------------------------- */

/*! Polynomial P0 - 3D */
template<Int N> struct trunk_P0_3d_A;
template<>      struct trunk_P0_3d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P0_3d_A    {static const UInt S = 1; typedef trunk_P0_3d_A<1> ROOT;};

/*! Polynomial P1 - 3D */
template<Int N> struct trunk_P1_3d_A;
template<>      struct trunk_P1_3d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1_3d_A<2> {typedef trunk_P1_3d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1_3d_A<3> {typedef trunk_P1_3d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1_3d_A<4> {typedef trunk_P1_3d_A<3> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P1_3d_A    {static const UInt S = 4; typedef trunk_P1_3d_A<4> ROOT;};

template<Int N> struct trunk_P1_3d_B;
template<>      struct trunk_P1_3d_B<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P1_3d_B    {static const UInt S = 1; typedef trunk_P1_3d_B<1> ROOT;};

template<Int N> struct trunk_P1_3d_C;
template<>      struct trunk_P1_3d_C<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P1_3d_C    {static const UInt S = 1; typedef trunk_P1_3d_C<1> ROOT;};

template<Int N> struct trunk_P1_3d_D;
template<>      struct trunk_P1_3d_D<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P1_3d_D    {static const UInt S = 1; typedef trunk_P1_3d_D<1> ROOT;};     

      
/*! Polynomial P2 - 3D */
template<Int N> struct trunk_P2_3d_A;
template<>      struct trunk_P2_3d_A<1>  {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_A<2>  {typedef trunk_P2_3d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_A<3>  {typedef trunk_P2_3d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_A<4>  {typedef trunk_P2_3d_A<3> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_P2_3d_A<5>  {typedef trunk_P2_3d_A<4> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_A<6>  {typedef trunk_P2_3d_A<5> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_A<7>  {typedef trunk_P2_3d_A<6> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_P2_3d_A<8>  {typedef trunk_P2_3d_A<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_A<9>  {typedef trunk_P2_3d_A<8> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_P2_3d_A<10> {typedef trunk_P2_3d_A<9> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       P2_3d_A     {static const UInt S = 10; typedef trunk_P2_3d_A<10> ROOT;};

template<Int N> struct trunk_P2_3d_B;
template<>      struct trunk_P2_3d_B<1> {                              static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_B<2> {typedef trunk_P2_3d_B<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P2_3d_B    {static const UInt S = 2; typedef trunk_P2_3d_B<2> ROOT;};

template<Int N> struct trunk_P2_3d_C;
template<>      struct trunk_P2_3d_C<1> {                              static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_C<2> {typedef trunk_P2_3d_C<1> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P2_3d_C    {static const UInt S = 2; typedef trunk_P2_3d_C<2> ROOT;};      

template<Int N> struct trunk_P2_3d_D;
template<>      struct trunk_P2_3d_D<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_P2_3d_D<2> {typedef trunk_P2_3d_D<1> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P2_3d_D    {static const UInt S = 2; typedef trunk_P2_3d_D<2> ROOT;};       

template<Int N> struct trunk_P2_3d_E;
template<>      struct trunk_P2_3d_E<1> {                              static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_E<2> {typedef trunk_P2_3d_E<1> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_E<3> {typedef trunk_P2_3d_E<2> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_E<4> {typedef trunk_P2_3d_E<3> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P2_3d_E    {static const UInt S = 4; typedef trunk_P2_3d_E<4> ROOT;};      

template<Int N> struct trunk_P2_3d_F;
template<>      struct trunk_P2_3d_F<1> {static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P2_3d_F    {static const UInt S = 1; typedef trunk_P2_3d_F<1> ROOT;}; 
     
template<Int N> struct trunk_P2_3d_G;
template<>      struct trunk_P2_3d_G<1> {                              static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_G<2> {typedef trunk_P2_3d_G<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_G<3> {typedef trunk_P2_3d_G<2> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P2_3d_G<4> {typedef trunk_P2_3d_G<3> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       P2_3d_G    {static const UInt S = 4; typedef trunk_P2_3d_G<4> ROOT;}; 
      
template<Int N> struct trunk_P2_3d_H;
template<>      struct trunk_P2_3d_H<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_P2_3d_H<2> {typedef trunk_P2_3d_H<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_P2_3d_H<3> {typedef trunk_P2_3d_H<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_P2_3d_H<4> {typedef trunk_P2_3d_H<3> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
struct                       P2_3d_H    {static const UInt S = 4; typedef trunk_P2_3d_H<4> ROOT;};
      
template<Int N> struct trunk_P2_3d_I;
template<>      struct trunk_P2_3d_I<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P2_3d_I    {static const UInt S = 1; typedef trunk_P2_3d_I<1> ROOT;}; 

template<Int N> struct trunk_P2_3d_L;
template<>      struct trunk_P2_3d_L<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       P2_3d_L    {static const UInt S = 1; typedef trunk_P2_3d_L<1> ROOT;}; 


/*! Polynomial Q0 - 3D */
template<Int N> struct trunk_Q0_3d_A;
template<>      struct trunk_Q0_3d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       Q0_3d_A    {static const UInt S = 1; typedef trunk_Q0_3d_A<1> ROOT;}; 


/*! Polynomial Q1 - 3D */
template<Int N> struct trunk_Q1_3d_A;
template<>      struct trunk_Q1_3d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_A<2> {typedef trunk_Q1_3d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_A<3> {typedef trunk_Q1_3d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_A<4> {typedef trunk_Q1_3d_A<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_A<5> {typedef trunk_Q1_3d_A<4> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_A<6> {typedef trunk_Q1_3d_A<5> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_A<7> {typedef trunk_Q1_3d_A<6> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_A<8> {typedef trunk_Q1_3d_A<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_A    {static const UInt S = 8; typedef trunk_Q1_3d_A<8> ROOT;}; 

template<Int N> struct trunk_Q1_3d_B;
template<>      struct trunk_Q1_3d_B<1> {                              static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_B<2> {typedef trunk_Q1_3d_B<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_B<3> {typedef trunk_Q1_3d_B<2> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_B<4> {typedef trunk_Q1_3d_B<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_B    {static const UInt S = 4; typedef trunk_Q1_3d_B<4> ROOT;};

template<Int N> struct trunk_Q1_3d_C;
template<>      struct trunk_Q1_3d_C<1> {                              static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_C<2> {typedef trunk_Q1_3d_C<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_C    {static const UInt S = 2; typedef trunk_Q1_3d_C<2> ROOT;};

template<Int N> struct trunk_Q1_3d_D;
template<>      struct trunk_Q1_3d_D<1> {                              static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_D<2> {typedef trunk_Q1_3d_D<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q1_3d_D<3> {typedef trunk_Q1_3d_D<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_D<4> {typedef trunk_Q1_3d_D<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_D    {static const UInt S = 4; typedef trunk_Q1_3d_D<4> ROOT;};

template<Int N> struct trunk_Q1_3d_E;
template<>      struct trunk_Q1_3d_E<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_E<2> {typedef trunk_Q1_3d_E<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_E<3> {typedef trunk_Q1_3d_E<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_E<4> {typedef trunk_Q1_3d_E<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_E    {static const UInt S = 4; typedef trunk_Q1_3d_E<4> ROOT;};

template<Int N> struct trunk_Q1_3d_F;
template<>      struct trunk_Q1_3d_F<1> {                              static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_F<2> {typedef trunk_Q1_3d_F<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_F    {static const UInt S = 2; typedef trunk_Q1_3d_F<2> ROOT;};

template<Int N> struct trunk_Q1_3d_G;
template<>      struct trunk_Q1_3d_G<1> {static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_G    {static const UInt S = 1; typedef trunk_Q1_3d_G<1> ROOT;};

template<Int N> struct trunk_Q1_3d_H;
template<>      struct trunk_Q1_3d_H<1> {                              static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q1_3d_H<2> {typedef trunk_Q1_3d_H<1> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q1_3d_H    {static const UInt S = 2; typedef trunk_Q1_3d_H<2> ROOT;};


/*! Polynomial Q2 - 3D */
template<Int N> struct trunk_Q2_3d_AA;
template<>      struct trunk_Q2_3d_AA<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<2>  {typedef trunk_Q2_3d_AA<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<3>  {typedef trunk_Q2_3d_AA<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<4>  {typedef trunk_Q2_3d_AA<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<5>  {typedef trunk_Q2_3d_AA<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<6>  {typedef trunk_Q2_3d_AA<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<7>  {typedef trunk_Q2_3d_AA<6>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<8>  {typedef trunk_Q2_3d_AA<7>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<9>  {typedef trunk_Q2_3d_AA<8>  SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AA<10> {typedef trunk_Q2_3d_AA<9>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<11> {typedef trunk_Q2_3d_AA<10> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<12> {typedef trunk_Q2_3d_AA<11> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<13> {typedef trunk_Q2_3d_AA<12> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<14> {typedef trunk_Q2_3d_AA<13> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<15> {typedef trunk_Q2_3d_AA<14> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<16> {typedef trunk_Q2_3d_AA<15> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<17> {typedef trunk_Q2_3d_AA<16> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<18> {typedef trunk_Q2_3d_AA<17> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AA<19> {typedef trunk_Q2_3d_AA<18> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<20> {typedef trunk_Q2_3d_AA<19> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<21> {typedef trunk_Q2_3d_AA<20> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<22> {typedef trunk_Q2_3d_AA<21> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<23> {typedef trunk_Q2_3d_AA<22> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<24> {typedef trunk_Q2_3d_AA<23> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<25> {typedef trunk_Q2_3d_AA<24> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<26> {typedef trunk_Q2_3d_AA<25> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AA<27> {typedef trunk_Q2_3d_AA<26> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AA    {static const UInt S = 27; typedef trunk_Q2_3d_AA<27> ROOT;}; 

template<Int N> struct trunk_Q2_3d_AB;
template<>      struct trunk_Q2_3d_AB<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AB<2>  {typedef trunk_Q2_3d_AB<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AB<3>  {typedef trunk_Q2_3d_AB<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AB<4>  {typedef trunk_Q2_3d_AB<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AB<5>  {typedef trunk_Q2_3d_AB<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AB<6>  {typedef trunk_Q2_3d_AB<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AB<7>  {typedef trunk_Q2_3d_AB<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AB<8>  {typedef trunk_Q2_3d_AB<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AB<9>  {typedef trunk_Q2_3d_AB<8>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AB<10> {typedef trunk_Q2_3d_AB<9>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AB<11> {typedef trunk_Q2_3d_AB<10> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AB<12> {typedef trunk_Q2_3d_AB<11> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AB<13> {typedef trunk_Q2_3d_AB<12> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AB<14> {typedef trunk_Q2_3d_AB<13> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AB<15> {typedef trunk_Q2_3d_AB<14> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AB<16> {typedef trunk_Q2_3d_AB<15> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AB<17> {typedef trunk_Q2_3d_AB<16> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AB<18> {typedef trunk_Q2_3d_AB<17> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AB    {static const UInt S = 18; typedef trunk_Q2_3d_AB<18> ROOT;};

template<Int N> struct trunk_Q2_3d_AC;
template<>      struct trunk_Q2_3d_AC<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AC<2>  {typedef trunk_Q2_3d_AC<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AC<3>  {typedef trunk_Q2_3d_AC<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AC<4>  {typedef trunk_Q2_3d_AC<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AC<5>  {typedef trunk_Q2_3d_AC<4>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AC<6>  {typedef trunk_Q2_3d_AC<5>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AC<7>  {typedef trunk_Q2_3d_AC<6>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AC<8>  {typedef trunk_Q2_3d_AC<7>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AC<9>  {typedef trunk_Q2_3d_AC<8>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AC<10> {typedef trunk_Q2_3d_AC<9>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AC<11> {typedef trunk_Q2_3d_AC<10> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AC<12> {typedef trunk_Q2_3d_AC<11> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AC     {static const UInt S = 12; typedef trunk_Q2_3d_AC<12> ROOT;};

template<Int N> struct trunk_Q2_3d_AD;
template<>      struct trunk_Q2_3d_AD<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AD<2>  {typedef trunk_Q2_3d_AD<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AD<3>  {typedef trunk_Q2_3d_AD<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AD<4>  {typedef trunk_Q2_3d_AD<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AD<5>  {typedef trunk_Q2_3d_AD<4>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AD<6>  {typedef trunk_Q2_3d_AD<5>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AD<7>  {typedef trunk_Q2_3d_AD<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AD<8>  {typedef trunk_Q2_3d_AD<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AD<9>  {typedef trunk_Q2_3d_AD<8>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AD<10> {typedef trunk_Q2_3d_AD<9>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AD<11> {typedef trunk_Q2_3d_AD<10> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AD<12> {typedef trunk_Q2_3d_AD<11> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AD<13> {typedef trunk_Q2_3d_AD<12> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AD<14> {typedef trunk_Q2_3d_AD<13> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AD<15> {typedef trunk_Q2_3d_AD<14> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AD<16> {typedef trunk_Q2_3d_AD<15> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AD<17> {typedef trunk_Q2_3d_AD<16> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AD<18> {typedef trunk_Q2_3d_AD<17> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AD    {static const UInt S = 18; typedef trunk_Q2_3d_AD<18> ROOT;};


template<Int N> struct trunk_Q2_3d_AE;
template<>      struct trunk_Q2_3d_AE<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AE<2>  {typedef trunk_Q2_3d_AE<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AE<3>  {typedef trunk_Q2_3d_AE<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AE<4>  {typedef trunk_Q2_3d_AE<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AE<5>  {typedef trunk_Q2_3d_AE<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AE<6>  {typedef trunk_Q2_3d_AE<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AE<7>  {typedef trunk_Q2_3d_AE<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AE<8>  {typedef trunk_Q2_3d_AE<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AE<9>  {typedef trunk_Q2_3d_AE<8>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AE<10> {typedef trunk_Q2_3d_AE<9>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AE<11> {typedef trunk_Q2_3d_AE<10> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AE<12> {typedef trunk_Q2_3d_AE<11> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AE<13> {typedef trunk_Q2_3d_AE<12> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AE<14> {typedef trunk_Q2_3d_AE<13> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AE<15> {typedef trunk_Q2_3d_AE<14> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AE<16> {typedef trunk_Q2_3d_AE<15> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AE<17> {typedef trunk_Q2_3d_AE<16> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AE<18> {typedef trunk_Q2_3d_AE<17> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AE     {static const UInt S = 18; typedef trunk_Q2_3d_AE<18> ROOT;};

template<Int N> struct trunk_Q2_3d_AF;
template<>      struct trunk_Q2_3d_AF<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AF<2>  {typedef trunk_Q2_3d_AF<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AF<3>  {typedef trunk_Q2_3d_AF<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AF<4>  {typedef trunk_Q2_3d_AF<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AF<5>  {typedef trunk_Q2_3d_AF<4>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AF<6>  {typedef trunk_Q2_3d_AF<5>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AF<7>  {typedef trunk_Q2_3d_AF<6>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AF<8>  {typedef trunk_Q2_3d_AF<7>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AF<9>  {typedef trunk_Q2_3d_AF<8>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AF<10> {typedef trunk_Q2_3d_AF<9>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AF<11> {typedef trunk_Q2_3d_AF<10> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AF<12> {typedef trunk_Q2_3d_AF<11> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AF     {static const UInt S = 12; typedef trunk_Q2_3d_AF<12> ROOT;};

template<Int N> struct trunk_Q2_3d_AG;
template<>      struct trunk_Q2_3d_AG<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AG<2>  {typedef trunk_Q2_3d_AG<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AG<3>  {typedef trunk_Q2_3d_AG<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AG<4>  {typedef trunk_Q2_3d_AG<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AG<5>  {typedef trunk_Q2_3d_AG<4>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AG<6>  {typedef trunk_Q2_3d_AG<5>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AG<7>  {typedef trunk_Q2_3d_AG<6>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AG<8>  {typedef trunk_Q2_3d_AG<7>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AG<9>  {typedef trunk_Q2_3d_AG<8>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AG<10> {typedef trunk_Q2_3d_AG<9>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AG<11> {typedef trunk_Q2_3d_AG<10> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AG<12> {typedef trunk_Q2_3d_AG<11> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AG     {static const UInt S = 12; typedef trunk_Q2_3d_AG<12> ROOT;};
 
template<Int N> struct trunk_Q2_3d_AH;
template<>      struct trunk_Q2_3d_AH<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AH<2>  {typedef trunk_Q2_3d_AH<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AH<3>  {typedef trunk_Q2_3d_AH<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AH<4>  {typedef trunk_Q2_3d_AH<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AH<5>  {typedef trunk_Q2_3d_AH<4>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AH<6>  {typedef trunk_Q2_3d_AH<5>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AH<7>  {typedef trunk_Q2_3d_AH<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AH<8>  {typedef trunk_Q2_3d_AH<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AH<9>  {typedef trunk_Q2_3d_AH<8>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AH<10> {typedef trunk_Q2_3d_AH<9>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AH<11> {typedef trunk_Q2_3d_AH<10> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AH<12> {typedef trunk_Q2_3d_AH<11> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AH<13> {typedef trunk_Q2_3d_AH<12> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AH<14> {typedef trunk_Q2_3d_AH<13> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AH<15> {typedef trunk_Q2_3d_AH<14> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AH<16> {typedef trunk_Q2_3d_AH<15> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AH<17> {typedef trunk_Q2_3d_AH<16> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AH<18> {typedef trunk_Q2_3d_AH<17> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AH     {static const UInt S = 18; typedef trunk_Q2_3d_AH<18> ROOT;};

template<Int N> struct trunk_Q2_3d_AI;
template<>      struct trunk_Q2_3d_AI<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AI<2>  {typedef trunk_Q2_3d_AI<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AI<3>  {typedef trunk_Q2_3d_AI<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AI<4>  {typedef trunk_Q2_3d_AI<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_Q2_3d_AI<5>  {typedef trunk_Q2_3d_AI<4>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AI<6>  {typedef trunk_Q2_3d_AI<5>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AI<7>  {typedef trunk_Q2_3d_AI<6>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AI<8>  {typedef trunk_Q2_3d_AI<7>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_AI<9>  {typedef trunk_Q2_3d_AI<8>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AI<10> {typedef trunk_Q2_3d_AI<9>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AI<11> {typedef trunk_Q2_3d_AI<10> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_AI<12> {typedef trunk_Q2_3d_AI<11> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_AI     {static const UInt S = 12; typedef trunk_Q2_3d_AI<12> ROOT;};

template<Int N> struct trunk_Q2_3d_BA;
template<>      struct trunk_Q2_3d_BA<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<2>  {typedef trunk_Q2_3d_BA<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<3>  {typedef trunk_Q2_3d_BA<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<4>  {typedef trunk_Q2_3d_BA<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<5>  {typedef trunk_Q2_3d_BA<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<6>  {typedef trunk_Q2_3d_BA<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<7>  {typedef trunk_Q2_3d_BA<6>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<8>  {typedef trunk_Q2_3d_BA<7>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<9>  {typedef trunk_Q2_3d_BA<8>  SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BA<10> {typedef trunk_Q2_3d_BA<9>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<11> {typedef trunk_Q2_3d_BA<10> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<12> {typedef trunk_Q2_3d_BA<11> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<13> {typedef trunk_Q2_3d_BA<12> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<14> {typedef trunk_Q2_3d_BA<13> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<15> {typedef trunk_Q2_3d_BA<14> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<16> {typedef trunk_Q2_3d_BA<15> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<17> {typedef trunk_Q2_3d_BA<16> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BA<18> {typedef trunk_Q2_3d_BA<17> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BA     {static const UInt S = 18; typedef trunk_Q2_3d_BA<18> ROOT;};

template<Int N> struct trunk_Q2_3d_BB;
template<>      struct trunk_Q2_3d_BB<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BB<2>  {typedef trunk_Q2_3d_BB<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BB<3>  {typedef trunk_Q2_3d_BB<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BB<4>  {typedef trunk_Q2_3d_BB<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BB<5>  {typedef trunk_Q2_3d_BB<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BB<6>  {typedef trunk_Q2_3d_BB<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BB<7>  {typedef trunk_Q2_3d_BB<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BB<8>  {typedef trunk_Q2_3d_BB<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BB<9>  {typedef trunk_Q2_3d_BB<8>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BB<10> {typedef trunk_Q2_3d_BB<9>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BB<11> {typedef trunk_Q2_3d_BB<10> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BB<12> {typedef trunk_Q2_3d_BB<11> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BB     {static const UInt S = 12; typedef trunk_Q2_3d_BB<12> ROOT;};

template<Int N> struct trunk_Q2_3d_BC;
template<>      struct trunk_Q2_3d_BC<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BC<2> {typedef trunk_Q2_3d_BC<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BC<3> {typedef trunk_Q2_3d_BC<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BC<4> {typedef trunk_Q2_3d_BC<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BC<5> {typedef trunk_Q2_3d_BC<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BC<6> {typedef trunk_Q2_3d_BC<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BC<7> {typedef trunk_Q2_3d_BC<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BC<8> {typedef trunk_Q2_3d_BC<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BC    {static const UInt S = 8; typedef trunk_Q2_3d_BC<8> ROOT;};

template<Int N> struct trunk_Q2_3d_BD;
template<>      struct trunk_Q2_3d_BD<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BD<2>  {typedef trunk_Q2_3d_BD<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BD<3>  {typedef trunk_Q2_3d_BD<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BD<4>  {typedef trunk_Q2_3d_BD<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BD<5>  {typedef trunk_Q2_3d_BD<4>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BD<6>  {typedef trunk_Q2_3d_BD<5>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BD<7>  {typedef trunk_Q2_3d_BD<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BD<8>  {typedef trunk_Q2_3d_BD<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BD<9>  {typedef trunk_Q2_3d_BD<8>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BD<10> {typedef trunk_Q2_3d_BD<9>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BD<11> {typedef trunk_Q2_3d_BD<10> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BD<12> {typedef trunk_Q2_3d_BD<11> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BD     {static const UInt S = 12; typedef trunk_Q2_3d_BD<12> ROOT;};

template<Int N> struct trunk_Q2_3d_BE;
template<>      struct trunk_Q2_3d_BE<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BE<2>  {typedef trunk_Q2_3d_BE<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BE<3>  {typedef trunk_Q2_3d_BE<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BE<4>  {typedef trunk_Q2_3d_BE<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BE<5>  {typedef trunk_Q2_3d_BE<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BE<6>  {typedef trunk_Q2_3d_BE<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BE<7>  {typedef trunk_Q2_3d_BE<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BE<8>  {typedef trunk_Q2_3d_BE<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BE<9>  {typedef trunk_Q2_3d_BE<8>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BE<10> {typedef trunk_Q2_3d_BE<9>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BE<11> {typedef trunk_Q2_3d_BE<10> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BE<12> {typedef trunk_Q2_3d_BE<11> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BE    {static const UInt S = 12; typedef trunk_Q2_3d_BE<12> ROOT;};

template<Int N> struct trunk_Q2_3d_BF;
template<>      struct trunk_Q2_3d_BF<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BF<2> {typedef trunk_Q2_3d_BF<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BF<3> {typedef trunk_Q2_3d_BF<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BF<4> {typedef trunk_Q2_3d_BF<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BF<5> {typedef trunk_Q2_3d_BF<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BF<6> {typedef trunk_Q2_3d_BF<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BF<7> {typedef trunk_Q2_3d_BF<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BF<8> {typedef trunk_Q2_3d_BF<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BF    {static const UInt S = 8; typedef trunk_Q2_3d_BF<8> ROOT;};

template<Int N> struct trunk_Q2_3d_BG;
template<>      struct trunk_Q2_3d_BG<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BG<2> {typedef trunk_Q2_3d_BG<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BG<3> {typedef trunk_Q2_3d_BG<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BG<4> {typedef trunk_Q2_3d_BG<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BG<5> {typedef trunk_Q2_3d_BG<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BG<6> {typedef trunk_Q2_3d_BG<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BG<7> {typedef trunk_Q2_3d_BG<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BG<8> {typedef trunk_Q2_3d_BG<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BG    {static const UInt S = 8; typedef trunk_Q2_3d_BG<8> ROOT;};

template<Int N> struct trunk_Q2_3d_BH;
template<>      struct trunk_Q2_3d_BH<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BH<2>  {typedef trunk_Q2_3d_BH<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BH<3>  {typedef trunk_Q2_3d_BH<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BH<4>  {typedef trunk_Q2_3d_BH<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BH<5>  {typedef trunk_Q2_3d_BH<4>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BH<6>  {typedef trunk_Q2_3d_BH<5>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BH<7>  {typedef trunk_Q2_3d_BH<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BH<8>  {typedef trunk_Q2_3d_BH<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BH<9>  {typedef trunk_Q2_3d_BH<8>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BH<10> {typedef trunk_Q2_3d_BH<9>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BH<11> {typedef trunk_Q2_3d_BH<10> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BH<12> {typedef trunk_Q2_3d_BH<11> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BH     {static const UInt S = 12; typedef trunk_Q2_3d_BH<12> ROOT;};

template<Int N> struct trunk_Q2_3d_BI;
template<>      struct trunk_Q2_3d_BI<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BI<2> {typedef trunk_Q2_3d_BI<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BI<3> {typedef trunk_Q2_3d_BI<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BI<4> {typedef trunk_Q2_3d_BI<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_BI<5> {typedef trunk_Q2_3d_BI<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BI<6> {typedef trunk_Q2_3d_BI<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BI<7> {typedef trunk_Q2_3d_BI<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_BI<8> {typedef trunk_Q2_3d_BI<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
struct                       Q2_3d_BI    {static const UInt S = 8; typedef trunk_Q2_3d_BI<8> ROOT;};

template<Int N> struct trunk_Q2_3d_CA;
template<>      struct trunk_Q2_3d_CA<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<2>  {typedef trunk_Q2_3d_CA<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<3>  {typedef trunk_Q2_3d_CA<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<4>  {typedef trunk_Q2_3d_CA<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<5>  {typedef trunk_Q2_3d_CA<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<6>  {typedef trunk_Q2_3d_CA<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<7>  {typedef trunk_Q2_3d_CA<6>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<8>  {typedef trunk_Q2_3d_CA<7>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<9>  {typedef trunk_Q2_3d_CA<8>  SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CA<10> {typedef trunk_Q2_3d_CA<9>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<11> {typedef trunk_Q2_3d_CA<10> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<12> {typedef trunk_Q2_3d_CA<11> SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<13> {typedef trunk_Q2_3d_CA<12> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<14> {typedef trunk_Q2_3d_CA<13> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<15> {typedef trunk_Q2_3d_CA<14> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<16> {typedef trunk_Q2_3d_CA<15> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<17> {typedef trunk_Q2_3d_CA<16> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CA<18> {typedef trunk_Q2_3d_CA<17> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CA     {static const UInt S = 18; typedef trunk_Q2_3d_CA<18> ROOT;};

template<Int N> struct trunk_Q2_3d_CB;
template<>      struct trunk_Q2_3d_CB<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CB<2>  {typedef trunk_Q2_3d_CB<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CB<3>  {typedef trunk_Q2_3d_CB<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CB<4>  {typedef trunk_Q2_3d_CB<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CB<5>  {typedef trunk_Q2_3d_CB<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CB<6>  {typedef trunk_Q2_3d_CB<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CB<7>  {typedef trunk_Q2_3d_CB<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CB<8>  {typedef trunk_Q2_3d_CB<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CB<9>  {typedef trunk_Q2_3d_CB<8>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CB<10> {typedef trunk_Q2_3d_CB<9>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CB<11> {typedef trunk_Q2_3d_CB<10> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CB<12> {typedef trunk_Q2_3d_CB<11> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CB     {static const UInt S = 12; typedef trunk_Q2_3d_CB<12> ROOT;};

template<Int N> struct trunk_Q2_3d_CC;
template<>      struct trunk_Q2_3d_CC<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CC<2> {typedef trunk_Q2_3d_CC<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CC<3> {typedef trunk_Q2_3d_CC<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CC<4> {typedef trunk_Q2_3d_CC<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CC<5> {typedef trunk_Q2_3d_CC<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CC<6> {typedef trunk_Q2_3d_CC<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CC<7> {typedef trunk_Q2_3d_CC<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CC<8> {typedef trunk_Q2_3d_CC<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CC    {static const UInt S = 8; typedef trunk_Q2_3d_CC<8> ROOT;};

template<Int N> struct trunk_Q2_3d_CD;
template<>      struct trunk_Q2_3d_CD<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CD<2>  {typedef trunk_Q2_3d_CD<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CD<3>  {typedef trunk_Q2_3d_CD<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CD<4>  {typedef trunk_Q2_3d_CD<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CD<5>  {typedef trunk_Q2_3d_CD<4>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CD<6>  {typedef trunk_Q2_3d_CD<5>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CD<7>  {typedef trunk_Q2_3d_CD<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CD<8>  {typedef trunk_Q2_3d_CD<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CD<9>  {typedef trunk_Q2_3d_CD<8>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CD<10> {typedef trunk_Q2_3d_CD<9>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CD<11> {typedef trunk_Q2_3d_CD<10> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CD<12> {typedef trunk_Q2_3d_CD<11> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CD     {static const UInt S = 12; typedef trunk_Q2_3d_CD<12> ROOT;};

template<Int N> struct trunk_Q2_3d_CE;
template<>      struct trunk_Q2_3d_CE<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CE<2>  {typedef trunk_Q2_3d_CE<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CE<3>  {typedef trunk_Q2_3d_CE<2>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CE<4>  {typedef trunk_Q2_3d_CE<3>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CE<5>  {typedef trunk_Q2_3d_CE<4>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CE<6>  {typedef trunk_Q2_3d_CE<5>  SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CE<7>  {typedef trunk_Q2_3d_CE<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CE<8>  {typedef trunk_Q2_3d_CE<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CE<9>  {typedef trunk_Q2_3d_CE<8>  SON; static SInt Px = 2; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CE<10> {typedef trunk_Q2_3d_CE<9>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CE<11> {typedef trunk_Q2_3d_CE<10> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CE<12> {typedef trunk_Q2_3d_CE<11> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CE     {static const UInt S = 12; typedef trunk_Q2_3d_CE<12> ROOT;};

template<Int N> struct trunk_Q2_3d_CF;
template<>      struct trunk_Q2_3d_CF<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CF<2> {typedef trunk_Q2_3d_CF<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CF<3> {typedef trunk_Q2_3d_CF<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CF<4> {typedef trunk_Q2_3d_CF<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CF<5> {typedef trunk_Q2_3d_CF<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CF<6> {typedef trunk_Q2_3d_CF<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CF<7> {typedef trunk_Q2_3d_CF<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CF<8> {typedef trunk_Q2_3d_CF<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CF    {static const UInt S = 8; typedef trunk_Q2_3d_CF<8> ROOT;};

template<Int N> struct trunk_Q2_3d_CG;
template<>      struct trunk_Q2_3d_CG<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CG<2> {typedef trunk_Q2_3d_CG<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CG<3> {typedef trunk_Q2_3d_CG<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CG<4> {typedef trunk_Q2_3d_CG<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CG<5> {typedef trunk_Q2_3d_CG<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CG<6> {typedef trunk_Q2_3d_CG<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CG<7> {typedef trunk_Q2_3d_CG<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CG<8> {typedef trunk_Q2_3d_CG<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CG    {static const UInt S = 8; typedef trunk_Q2_3d_CG<8> ROOT;};

template<Int N> struct trunk_Q2_3d_CH;
template<>      struct trunk_Q2_3d_CH<1>  {                                static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CH<2>  {typedef trunk_Q2_3d_CH<1>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CH<3>  {typedef trunk_Q2_3d_CH<2>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CH<4>  {typedef trunk_Q2_3d_CH<3>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CH<5>  {typedef trunk_Q2_3d_CH<4>  SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CH<6>  {typedef trunk_Q2_3d_CH<5>  SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CH<7>  {typedef trunk_Q2_3d_CH<6>  SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CH<8>  {typedef trunk_Q2_3d_CH<7>  SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CH<9>  {typedef trunk_Q2_3d_CH<8>  SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CH<10> {typedef trunk_Q2_3d_CH<9>  SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CH<11> {typedef trunk_Q2_3d_CH<10> SON; static SInt Px = 0; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CH<12> {typedef trunk_Q2_3d_CH<11> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CH     {static const UInt S = 12; typedef trunk_Q2_3d_CH<12> ROOT;};

template<Int N> struct trunk_Q2_3d_CI;
template<>      struct trunk_Q2_3d_CI<1> {                               static SInt Px = 2; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CI<2> {typedef trunk_Q2_3d_CI<1> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CI<3> {typedef trunk_Q2_3d_CI<2> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CI<4> {typedef trunk_Q2_3d_CI<3> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 2; static SReal C; };
template<>      struct trunk_Q2_3d_CI<5> {typedef trunk_Q2_3d_CI<4> SON; static SInt Px = 2; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CI<6> {typedef trunk_Q2_3d_CI<5> SON; static SInt Px = 2; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CI<7> {typedef trunk_Q2_3d_CI<6> SON; static SInt Px = 1; static SInt Py = 2; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_Q2_3d_CI<8> {typedef trunk_Q2_3d_CI<7> SON; static SInt Px = 1; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       Q2_3d_CI    {static const UInt S = 8; typedef trunk_Q2_3d_CI<8> ROOT;};



/*! Polynomial D0 - 3D */
template<Int N> struct trunk_D0_3d_A;
template<>      struct trunk_D0_3d_A<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       D0_3d_A    {static const UInt S = 1; typedef trunk_D0_3d_A<1> ROOT;};

template<Int N> struct trunk_D0_3d_B;
template<>      struct trunk_D0_3d_B<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       D0_3d_B    {static const UInt S = 1; typedef trunk_D0_3d_B<1> ROOT;};

template<Int N> struct trunk_D0_3d_C;
template<>      struct trunk_D0_3d_C<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       D0_3d_C    {static const UInt S = 1; typedef trunk_D0_3d_C<1> ROOT;};
      
template<Int N> struct trunk_D0_3d_D;
template<>      struct trunk_D0_3d_D<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       D0_3d_D    {static const UInt S = 1; typedef trunk_D0_3d_D<1> ROOT;};
      

/*! Polynomial P1BP1 - 3D */
template<Int N> struct trunk_P1BP1_3d_A;
template<>      struct trunk_P1BP1_3d_A<1> {                                 static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1BP1_3d_A<2> {typedef trunk_P1BP1_3d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1BP1_3d_A<3> {typedef trunk_P1BP1_3d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1BP1_3d_A<4> {typedef trunk_P1BP1_3d_A<3> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P1BP1_3d_A    {static const UInt S = 4; typedef trunk_P1BP1_3d_A<4> ROOT;};

template<Int N> struct trunk_P1BP1_3d_B;
template<>      struct trunk_P1BP1_3d_B<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P1BP1_3d_B    {static const UInt S = 1; typedef trunk_P1BP1_3d_B<1> ROOT;};

template<Int N> struct trunk_P1BP1_3d_C;
template<>      struct trunk_P1BP1_3d_C<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P1BP1_3d_C    {static const UInt S = 1; typedef trunk_P1BP1_3d_C<1> ROOT;};

template<Int N> struct trunk_P1BP1_3d_D;
template<>      struct trunk_P1BP1_3d_D<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P1BP1_3d_D    {static const UInt S = 1; typedef trunk_P1BP1_3d_D<1> ROOT;}; 

template<Int N> struct trunk_P1BP1_3d_E;
template<>      struct trunk_P1BP1_3d_E<1> {                                 static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1BP1_3d_E<2> {typedef trunk_P1BP1_3d_E<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1BP1_3d_E<3> {typedef trunk_P1BP1_3d_E<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_P1BP1_3d_E<4> {typedef trunk_P1BP1_3d_E<3> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P1BP1_3d_E    {static const UInt S = 4; typedef trunk_P1BP1_3d_E<4> ROOT;};

template<Int N> struct trunk_P1BP1_3d_F;
template<>      struct trunk_P1BP1_3d_F<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       P1BP1_3d_F    {static const UInt S = 1; typedef trunk_P1BP1_3d_F<1> ROOT;};
      
template<Int N> struct trunk_P1BP1_3d_G;
template<>      struct trunk_P1BP1_3d_G<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       P1BP1_3d_G    {static const UInt S = 1; typedef trunk_P1BP1_3d_G<1> ROOT;};
      
template<Int N> struct trunk_P1BP1_3d_H;
template<>      struct trunk_P1BP1_3d_H<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       P1BP1_3d_H    {static const UInt S = 1; typedef trunk_P1BP1_3d_H<1> ROOT;}; 


/*! Polynomial RT0 - 3D */
template<Int N> struct trunk_RT0_3d_Ax;
template<>      struct trunk_RT0_3d_Ax<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_Ax    {static const UInt S = 1; typedef trunk_RT0_3d_Ax<1> ROOT;}; 

template<Int N> struct trunk_RT0_3d_Ay;
template<>      struct trunk_RT0_3d_Ay<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_Ay    {static const UInt S = 1; typedef trunk_RT0_3d_Ay<1> ROOT;}; 

template<Int N> struct trunk_RT0_3d_Az;
template<>      struct trunk_RT0_3d_Az<1> {                                static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_RT0_3d_Az<2> {typedef trunk_RT0_3d_Az<1> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       RT0_3d_Az    {static const UInt S = 2; typedef trunk_RT0_3d_Az<2> ROOT;}; 

template<Int N> struct trunk_RT0_3d_Bx;
template<>      struct trunk_RT0_3d_Bx<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_Bx    {static const UInt S = 1; typedef trunk_RT0_3d_Bx<1> ROOT;};

template<Int N> struct trunk_RT0_3d_By;
template<>      struct trunk_RT0_3d_By<1> {                                static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_RT0_3d_By<2> {typedef trunk_RT0_3d_By<1> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_By    {static const UInt S = 2; typedef trunk_RT0_3d_By<2> ROOT;}; 
      
template<Int N> struct trunk_RT0_3d_Bz;
template<>      struct trunk_RT0_3d_Bz<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       RT0_3d_Bz    {static const UInt S = 1; typedef trunk_RT0_3d_Bz<1> ROOT;};      
      
template<Int N> struct trunk_RT0_3d_Cx;
template<>      struct trunk_RT0_3d_Cx<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_Cx    {static const UInt S = 1; typedef trunk_RT0_3d_Cx<1> ROOT;};  
      
template<Int N> struct trunk_RT0_3d_Cy;
template<>      struct trunk_RT0_3d_Cy<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_Cy    {static const UInt S = 1; typedef trunk_RT0_3d_Cy<1> ROOT;}; 
      
template<Int N> struct trunk_RT0_3d_Cz;
template<>      struct trunk_RT0_3d_Cz<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       RT0_3d_Cz    {static const UInt S = 1; typedef trunk_RT0_3d_Cz<1> ROOT;}; 
      
template<Int N> struct trunk_RT0_3d_Dx;
template<>      struct trunk_RT0_3d_Dx<1> {                                static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_RT0_3d_Dx<2> {typedef trunk_RT0_3d_Dx<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_Dx    {static const UInt S = 2; typedef trunk_RT0_3d_Dx<2> ROOT;};       

template<Int N> struct trunk_RT0_3d_Dy;
template<>      struct trunk_RT0_3d_Dy<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       RT0_3d_Dy    {static const UInt S = 1; typedef trunk_RT0_3d_Dy<1> ROOT;};       
 
template<Int N> struct trunk_RT0_3d_Dz;
template<>      struct trunk_RT0_3d_Dz<1> {static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       RT0_3d_Dz    {static const UInt S = 1; typedef trunk_RT0_3d_Dz<1> ROOT;};


/*! Polynomial W1 - 3D */
template<Int N> struct trunk_W1_3d_A;
template<>      struct trunk_W1_3d_A<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_W1_3d_A<2> {typedef trunk_W1_3d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_W1_3d_A<3> {typedef trunk_W1_3d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_W1_3d_A<4> {typedef trunk_W1_3d_A<3> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_W1_3d_A<5> {typedef trunk_W1_3d_A<4> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_W1_3d_A<6> {typedef trunk_W1_3d_A<5> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       W1_3d_A    {static const UInt S = 6; typedef trunk_W1_3d_A<6> ROOT;};

template<Int N> struct trunk_W1_3d_B;
template<>      struct trunk_W1_3d_B<1> {                              static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_W1_3d_B<2> {typedef trunk_W1_3d_B<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       W1_3d_B    {static const UInt S = 2; typedef trunk_W1_3d_B<2> ROOT;};

template<Int N> struct trunk_W1_3d_C;
template<>      struct trunk_W1_3d_C<1> {                              static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_W1_3d_C<2> {typedef trunk_W1_3d_C<1> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       W1_3d_C    {static const UInt S = 2; typedef trunk_W1_3d_C<2> ROOT;};

template<Int N> struct trunk_W1_3d_D;
template<>      struct trunk_W1_3d_D<1> {                              static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_W1_3d_D<2> {typedef trunk_W1_3d_D<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
template<>      struct trunk_W1_3d_D<3> {typedef trunk_W1_3d_D<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       W1_3d_D    {static const UInt S = 3; typedef trunk_W1_3d_D<3> ROOT;};

template<Int N> struct trunk_W1_3d_E;
template<>      struct trunk_W1_3d_E<1> {static SInt Px = 1; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       W1_3d_E    {static const UInt S = 1; typedef trunk_W1_3d_E<1> ROOT;};

template<Int N> struct trunk_W1_3d_F;
template<>      struct trunk_W1_3d_F<1> {static SInt Px = 0; static SInt Py = 1; static SInt Pz = 1; static SReal C; };
struct                       W1_3d_F    {static const UInt S = 1; typedef trunk_W1_3d_F<1> ROOT;};


/*! Polynomial CR1 - 3D */
template<Int N> struct trunk_CR1_3d_A;
template<>      struct trunk_CR1_3d_A<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_CR1_3d_A<2> {typedef trunk_CR1_3d_A<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_CR1_3d_A<3> {typedef trunk_CR1_3d_A<2> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_CR1_3d_A<4> {typedef trunk_CR1_3d_A<3> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       CR1_3d_A    {static const UInt S = 4; typedef trunk_CR1_3d_A<4> ROOT;};

template<Int N> struct trunk_CR1_3d_B;
template<>      struct trunk_CR1_3d_B<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_CR1_3d_B<2> {typedef trunk_CR1_3d_B<1> SON; static SInt Px = 1; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
struct                       CR1_3d_B    {static const UInt S = 2; typedef trunk_CR1_3d_B<2> ROOT;};

template<Int N> struct trunk_CR1_3d_C;
template<>      struct trunk_CR1_3d_C<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_CR1_3d_C<2> {typedef trunk_CR1_3d_C<1> SON; static SInt Px = 0; static SInt Py = 1; static SInt Pz = 0; static SReal C; };
struct                       CR1_3d_C    {static const UInt S = 2; typedef trunk_CR1_3d_C<2> ROOT;};

template<Int N> struct trunk_CR1_3d_D;
template<>      struct trunk_CR1_3d_D<1> {                               static SInt Px = 0; static SInt Py = 0; static SInt Pz = 0; static SReal C; };
template<>      struct trunk_CR1_3d_D<2> {typedef trunk_CR1_3d_D<1> SON; static SInt Px = 0; static SInt Py = 0; static SInt Pz = 1; static SReal C; };
struct                       CR1_3d_D    {static const UInt S = 2; typedef trunk_CR1_3d_D<2> ROOT;};

#endif
