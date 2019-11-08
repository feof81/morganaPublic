/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "polyCards.h"


/*! Polynomial P0 - 1D */
SReal trunk_P0_1d_A<1>::C = 1.0;

/*! Polynomial P1 - 1D */
SReal trunk_P1_1d_A<1>::C =  1.0;
SReal trunk_P1_1d_A<2>::C = -1.0;

SReal trunk_P1_1d_B<1>::C = 1.0;

/*! Polynomial P2 - 1D */
SReal trunk_P2_1d_A<1>::C =  1.0;
SReal trunk_P2_1d_A<2>::C = -3.0;
SReal trunk_P2_1d_A<3>::C =  2.0;

SReal trunk_P2_1d_B<1>::C = -1.0;
SReal trunk_P2_1d_B<2>::C =  2.0;

SReal trunk_P2_1d_C<1>::C =  4.0;
SReal trunk_P2_1d_C<2>::C = -4.0;


/*! Polynomial D0 - 1D */
SReal trunk_D0_1d_A<1>::C = 1.0;


/*! Polynomial DG0 - 1D */
SReal trunk_DG0_1d_A<1>::C = 1.0;


/*! Polynomial DG1 - 1D */
SReal trunk_DG1_1d_A<1>::C =  1.0;

SReal trunk_DG1_1d_B<1>::C = -1.0;
SReal trunk_DG1_1d_B<2>::C =  2.0;




/*! -------------------------------------------------------------------------------------------- */
/*!                                         POLYCARDS 2D                                         */
/*! -------------------------------------------------------------------------------------------- */

/*! Polynomial P0 - 2D */
SReal trunk_P0_2d_A<1>::C = 1.0;   


/*! Polynomial P1 - 2D */
SReal trunk_P1_2d_A<1>::C =  1.0;
SReal trunk_P1_2d_A<2>::C = -1.0;
SReal trunk_P1_2d_A<3>::C = -1.0;

SReal trunk_P1_2d_B<1>::C =  1.0;

SReal trunk_P1_2d_C<1>::C =  1.0;


/*! Polynomial P2 - 2D */
SReal trunk_P2_2d_A<1>::C =  1.0;
SReal trunk_P2_2d_A<2>::C = -3.0;
SReal trunk_P2_2d_A<3>::C = -3.0;
SReal trunk_P2_2d_A<4>::C =  2.0;
SReal trunk_P2_2d_A<5>::C =  2.0;
SReal trunk_P2_2d_A<6>::C =  4.0;

SReal trunk_P2_2d_B<1>::C =  2.0;
SReal trunk_P2_2d_B<2>::C = -1.0;

SReal trunk_P2_2d_C<1>::C =  2.0;
SReal trunk_P2_2d_C<2>::C = -1.0;    

SReal trunk_P2_2d_D<1>::C =  4.0;
SReal trunk_P2_2d_D<2>::C = -4.0;
SReal trunk_P2_2d_D<3>::C = -4.0;

SReal trunk_P2_2d_E<1>::C =  4.0;     

SReal trunk_P2_2d_F<1>::C =  4.0;
SReal trunk_P2_2d_F<2>::C = -4.0;
SReal trunk_P2_2d_F<3>::C = -4.0;       


/*! Polynomial Q0 - 2D */
SReal trunk_Q0_2d_A<1>::C = 1.0;


/*! Polynomial Q1 - 2D */
SReal trunk_Q1_2d_A<1>::C =  1.0;
SReal trunk_Q1_2d_A<2>::C = -1.0;
SReal trunk_Q1_2d_A<3>::C = -1.0;
SReal trunk_Q1_2d_A<4>::C =  1.0;

SReal trunk_Q1_2d_B<1>::C =  1.0;
SReal trunk_Q1_2d_B<2>::C = -1.0;

SReal trunk_Q1_2d_C<1>::C =  1.0;

SReal trunk_Q1_2d_D<1>::C =  1.0;
SReal trunk_Q1_2d_D<2>::C = -1.0;


/*! Polynomial Q2 - 2D */
SReal trunk_Q2_2d_A<1>::C =  4.0;
SReal trunk_Q2_2d_A<2>::C = -6.0;
SReal trunk_Q2_2d_A<3>::C =  2.0;
SReal trunk_Q2_2d_A<4>::C = -6.0;
SReal trunk_Q2_2d_A<5>::C =  9.0;
SReal trunk_Q2_2d_A<6>::C = -3.0;
SReal trunk_Q2_2d_A<7>::C =  2.0;
SReal trunk_Q2_2d_A<8>::C = -3.0;
SReal trunk_Q2_2d_A<9>::C =  1.0;

SReal trunk_Q2_2d_B<1>::C =  4.0;
SReal trunk_Q2_2d_B<2>::C = -6.0;
SReal trunk_Q2_2d_B<3>::C =  2.0;
SReal trunk_Q2_2d_B<4>::C = -2.0;
SReal trunk_Q2_2d_B<5>::C =  3.0;
SReal trunk_Q2_2d_B<6>::C = -1.0;

SReal trunk_Q2_2d_C<1>::C =  4.0;
SReal trunk_Q2_2d_C<2>::C = -2.0;
SReal trunk_Q2_2d_C<3>::C = -2.0;
SReal trunk_Q2_2d_C<4>::C =  1.0;

SReal trunk_Q2_2d_D<1>::C =  4.0;
SReal trunk_Q2_2d_D<2>::C = -2.0;
SReal trunk_Q2_2d_D<3>::C = -6.0;
SReal trunk_Q2_2d_D<4>::C =  3.0;
SReal trunk_Q2_2d_D<5>::C =  2.0;
SReal trunk_Q2_2d_D<6>::C = -1.0;

SReal trunk_Q2_2d_E<1>::C =  -8.0;
SReal trunk_Q2_2d_E<2>::C =  12.0;
SReal trunk_Q2_2d_E<3>::C =  -4.0;
SReal trunk_Q2_2d_E<4>::C =   8.0;
SReal trunk_Q2_2d_E<5>::C = -12.0;
SReal trunk_Q2_2d_E<6>::C =   4.0;

SReal trunk_Q2_2d_F<1>::C = -8.0;
SReal trunk_Q2_2d_F<2>::C =  8.0;
SReal trunk_Q2_2d_F<3>::C =  4.0;
SReal trunk_Q2_2d_F<4>::C = -4.0;

SReal trunk_Q2_2d_G<1>::C = -8.0;
SReal trunk_Q2_2d_G<2>::C =  4.0;
SReal trunk_Q2_2d_G<3>::C =  8.0;
SReal trunk_Q2_2d_G<4>::C = -4.0;

SReal trunk_Q2_2d_H<1>::C =  -8.0;
SReal trunk_Q2_2d_H<2>::C =   8.0;
SReal trunk_Q2_2d_H<3>::C =  12.0;
SReal trunk_Q2_2d_H<4>::C = -12.0;
SReal trunk_Q2_2d_H<5>::C =  -4.0;
SReal trunk_Q2_2d_H<6>::C =   4.0;

SReal trunk_Q2_2d_I<1>::C =  16.0;
SReal trunk_Q2_2d_I<2>::C = -16.0;
SReal trunk_Q2_2d_I<3>::C = -16.0;
SReal trunk_Q2_2d_I<4>::C =  16.0;


/*! Polynomial D0 - 2D */
SReal trunk_D0_2d_A<1>::C =  1.0;    

SReal trunk_D0_2d_B<1>::C =  1.0;
      
SReal trunk_D0_2d_C<1>::C =  1.0;   
      

/*! Polynomial RT0 - 2D */
SReal trunk_RT0_2d_Ax<1>::C =  1.0;

SReal trunk_RT0_2d_Ay<1>::C = -1.0;
SReal trunk_RT0_2d_Ay<2>::C =  1.0;

SReal trunk_RT0_2d_Az<1>::C =  0.0;

SReal trunk_RT0_2d_Bx<1>::C =  1.0;

SReal trunk_RT0_2d_By<1>::C =  1.0;

SReal trunk_RT0_2d_Bz<1>::C =  0.0;

SReal trunk_RT0_2d_Cx<1>::C = -1.0;
SReal trunk_RT0_2d_Cx<2>::C =  1.0;
      
SReal trunk_RT0_2d_Cy<1>::C =  1.0;  

SReal trunk_RT0_2d_Cz<1>::C =  0.0;     


/*! Polynomial DG0 - 2D */
SReal trunk_DG0_2d_A<1>::C =  1.0;


/*! Polynomial DG1 - 2D */
SReal trunk_DG1_2d_A<1>::C =  1.0;
SReal trunk_DG1_2d_A<2>::C = -2.0;

SReal trunk_DG1_2d_B<1>::C = -1.0;
SReal trunk_DG1_2d_B<2>::C =  2.0;
SReal trunk_DG1_2d_B<3>::C =  2.0;

SReal trunk_DG1_2d_C<1>::C =  1.0;
SReal trunk_DG1_2d_C<2>::C = -2.0;




/*! -------------------------------------------------------------------------------------------- */
/*!                                         POLYCARDS 3D                                         */
/*! -------------------------------------------------------------------------------------------- */

/*! Polynomial P0 - 3D */
SReal trunk_P0_3d_A<1>::C =  1.0;

/*! Polynomial P1 - 3D */
SReal trunk_P1_3d_A<1>::C =  1.0;
SReal trunk_P1_3d_A<2>::C = -1.0;
SReal trunk_P1_3d_A<3>::C = -1.0;
SReal trunk_P1_3d_A<4>::C = -1.0;

SReal trunk_P1_3d_B<1>::C =  1.0;

SReal trunk_P1_3d_C<1>::C =  1.0;

SReal trunk_P1_3d_D<1>::C =  1.0;

      
/*! Polynomial P2 - 3D */
SReal trunk_P2_3d_A<1>::C =  1.0;
SReal trunk_P2_3d_A<2>::C = -3.0;
SReal trunk_P2_3d_A<3>::C = -3.0;
SReal trunk_P2_3d_A<4>::C = -3.0;
SReal trunk_P2_3d_A<5>::C =  2.0;
SReal trunk_P2_3d_A<6>::C =  2.0;
SReal trunk_P2_3d_A<7>::C =  2.0;
SReal trunk_P2_3d_A<8>::C =  4.0;
SReal trunk_P2_3d_A<9>::C =  4.0;
SReal trunk_P2_3d_A<10>::C =  4.0;

SReal trunk_P2_3d_B<1>::C =  2.0;
SReal trunk_P2_3d_B<2>::C = -1.0;

SReal trunk_P2_3d_C<1>::C =  2.0;
SReal trunk_P2_3d_C<2>::C = -1.0;    

SReal trunk_P2_3d_D<1>::C =  2.0;
SReal trunk_P2_3d_D<2>::C = -1.0;     

SReal trunk_P2_3d_E<1>::C =  4.0;
SReal trunk_P2_3d_E<2>::C = -4.0;
SReal trunk_P2_3d_E<3>::C = -4.0;
SReal trunk_P2_3d_E<4>::C = -4.0;   

SReal trunk_P2_3d_F<1>::C =  4.0;
     
SReal trunk_P2_3d_G<1>::C =  4.0;
SReal trunk_P2_3d_G<2>::C = -4.0;
SReal trunk_P2_3d_G<3>::C = -4.0;
SReal trunk_P2_3d_G<4>::C = -4.0;
      
SReal trunk_P2_3d_H<1>::C =  4.0;
SReal trunk_P2_3d_H<2>::C = -4.0;
SReal trunk_P2_3d_H<3>::C = -4.0;
SReal trunk_P2_3d_H<4>::C = -4.0;
      
SReal trunk_P2_3d_I<1>::C =  4.0;

SReal trunk_P2_3d_L<1>::C =  4.0;


/*! Polynomial Q0 - 3D */
SReal trunk_Q0_3d_A<1>::C =  1.0;


/*! Polynomial Q1 - 3D */
SReal trunk_Q1_3d_A<1>::C =  1.0;
SReal trunk_Q1_3d_A<2>::C = -1.0;
SReal trunk_Q1_3d_A<3>::C = -1.0;
SReal trunk_Q1_3d_A<4>::C =  1.0;
SReal trunk_Q1_3d_A<5>::C = -1.0;
SReal trunk_Q1_3d_A<6>::C =  1.0;
SReal trunk_Q1_3d_A<7>::C =  1.0;
SReal trunk_Q1_3d_A<8>::C = -1.0;

SReal trunk_Q1_3d_B<1>::C =  1.0;
SReal trunk_Q1_3d_B<2>::C = -1.0;
SReal trunk_Q1_3d_B<3>::C = -1.0;
SReal trunk_Q1_3d_B<4>::C =  1.0;

SReal trunk_Q1_3d_C<1>::C =  1.0;
SReal trunk_Q1_3d_C<2>::C = -1.0;

SReal trunk_Q1_3d_D<1>::C =  1.0;
SReal trunk_Q1_3d_D<2>::C = -1.0;
SReal trunk_Q1_3d_D<3>::C = -1.0;
SReal trunk_Q1_3d_D<4>::C =  1.0;

SReal trunk_Q1_3d_E<1>::C =  1.0;
SReal trunk_Q1_3d_E<2>::C = -1.0;
SReal trunk_Q1_3d_E<3>::C = -1.0;
SReal trunk_Q1_3d_E<4>::C =  1.0;

SReal trunk_Q1_3d_F<1>::C =  1.0;
SReal trunk_Q1_3d_F<2>::C = -1.0;

SReal trunk_Q1_3d_G<1>::C =  1.0;

SReal trunk_Q1_3d_H<1>::C =  1.0;
SReal trunk_Q1_3d_H<2>::C = -1.0;


/*! Polynomial Q2 - 3D */
SReal trunk_Q2_3d_AA<1>::C =  4.0;
SReal trunk_Q2_3d_AA<2>::C = -6.0;
SReal trunk_Q2_3d_AA<3>::C =  2.0;
SReal trunk_Q2_3d_AA<4>::C = -6.0;
SReal trunk_Q2_3d_AA<5>::C =  9.0;
SReal trunk_Q2_3d_AA<6>::C = -3.0;
SReal trunk_Q2_3d_AA<7>::C =  2.0;
SReal trunk_Q2_3d_AA<8>::C = -3.0;
SReal trunk_Q2_3d_AA<9>::C =  1.0;
SReal trunk_Q2_3d_AA<10>::C = -4.0 * 3.0;
SReal trunk_Q2_3d_AA<11>::C =  6.0 * 3.0;
SReal trunk_Q2_3d_AA<12>::C = -2.0 * 3.0;
SReal trunk_Q2_3d_AA<13>::C =  6.0 * 3.0;
SReal trunk_Q2_3d_AA<14>::C = -9.0 * 3.0;
SReal trunk_Q2_3d_AA<15>::C =  3.0 * 3.0;
SReal trunk_Q2_3d_AA<16>::C = -2.0 * 3.0;
SReal trunk_Q2_3d_AA<17>::C =  3.0 * 3.0;
SReal trunk_Q2_3d_AA<18>::C = -1.0 * 3.0;
SReal trunk_Q2_3d_AA<19>::C =  4.0 * 2;
SReal trunk_Q2_3d_AA<20>::C = -6.0 * 2;
SReal trunk_Q2_3d_AA<21>::C =  2.0 * 2;
SReal trunk_Q2_3d_AA<22>::C = -6.0 * 2;
SReal trunk_Q2_3d_AA<23>::C =  9.0 * 2;
SReal trunk_Q2_3d_AA<24>::C = -3.0 * 2;
SReal trunk_Q2_3d_AA<25>::C =  2.0 * 2;
SReal trunk_Q2_3d_AA<26>::C = -3.0 * 2;
SReal trunk_Q2_3d_AA<27>::C =  1.0 * 2;

SReal trunk_Q2_3d_AB<1>::C =  4.0;
SReal trunk_Q2_3d_AB<2>::C = -6.0;
SReal trunk_Q2_3d_AB<3>::C =  2.0;
SReal trunk_Q2_3d_AB<4>::C = -2.0;
SReal trunk_Q2_3d_AB<5>::C =  3.0;
SReal trunk_Q2_3d_AB<6>::C = -1.0;
SReal trunk_Q2_3d_AB<7>::C =  4.0 * (-3.0);
SReal trunk_Q2_3d_AB<8>::C = -6.0 * (-3.0);
SReal trunk_Q2_3d_AB<9>::C =  2.0 * (-3.0);
SReal trunk_Q2_3d_AB<10>::C = -2.0 * (-3.0);
SReal trunk_Q2_3d_AB<11>::C =  3.0 * (-3.0);
SReal trunk_Q2_3d_AB<12>::C = -1.0 * (-3.0);
SReal trunk_Q2_3d_AB<13>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_AB<14>::C = -6.0 * 2.0;
SReal trunk_Q2_3d_AB<15>::C =  2.0 * 2.0;
SReal trunk_Q2_3d_AB<16>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_AB<17>::C =  3.0 * 2.0;
SReal trunk_Q2_3d_AB<18>::C = -1.0 * 2.0;

SReal trunk_Q2_3d_AC<1>::C =  4.0;
SReal trunk_Q2_3d_AC<2>::C = -2.0;
SReal trunk_Q2_3d_AC<3>::C = -2.0;
SReal trunk_Q2_3d_AC<4>::C =  1.0;
SReal trunk_Q2_3d_AC<5>::C =  4.0 * (-3.0);
SReal trunk_Q2_3d_AC<6>::C = -2.0 * (-3.0);
SReal trunk_Q2_3d_AC<7>::C = -2.0 * (-3.0);
SReal trunk_Q2_3d_AC<8>::C =  1.0 * (-3.0);
SReal trunk_Q2_3d_AC<9>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_AC<10>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_AC<11>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_AC<12>::C =  1.0 * 2.0;

SReal trunk_Q2_3d_AD<1>::C =  4.0;
SReal trunk_Q2_3d_AD<2>::C = -2.0;
SReal trunk_Q2_3d_AD<3>::C = -6.0;
SReal trunk_Q2_3d_AD<4>::C =  3.0;
SReal trunk_Q2_3d_AD<5>::C =  2.0;
SReal trunk_Q2_3d_AD<6>::C = -1.0;
SReal trunk_Q2_3d_AD<7>::C =  4.0 * (-3.0);
SReal trunk_Q2_3d_AD<8>::C = -2.0 * (-3.0);
SReal trunk_Q2_3d_AD<9>::C = -6.0 * (-3.0);
SReal trunk_Q2_3d_AD<10>::C =  3.0 * (-3.0);
SReal trunk_Q2_3d_AD<11>::C =  2.0 * (-3.0);
SReal trunk_Q2_3d_AD<12>::C = -1.0 * (-3.0);
SReal trunk_Q2_3d_AD<13>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_AD<14>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_AD<15>::C = -6.0 * 2.0;
SReal trunk_Q2_3d_AD<16>::C =  3.0 * 2.0;
SReal trunk_Q2_3d_AD<17>::C =  2.0 * 2.0;
SReal trunk_Q2_3d_AD<18>::C = -1.0 * 2.0;


SReal trunk_Q2_3d_AE<1>::C =  -8.0;
SReal trunk_Q2_3d_AE<2>::C =  12.0;
SReal trunk_Q2_3d_AE<3>::C =  -4.0;
SReal trunk_Q2_3d_AE<4>::C =   8.0;
SReal trunk_Q2_3d_AE<5>::C = -12.0;
SReal trunk_Q2_3d_AE<6>::C =   4.0;
SReal trunk_Q2_3d_AE<7>::C =  -8.0 * (-3.0);
SReal trunk_Q2_3d_AE<8>::C =  12.0 * (-3.0);
SReal trunk_Q2_3d_AE<9>::C =  -4.0 * (-3.0);
SReal trunk_Q2_3d_AE<10>::C =   8.0 * (-3.0);
SReal trunk_Q2_3d_AE<11>::C = -12.0 * (-3.0);
SReal trunk_Q2_3d_AE<12>::C =   4.0 * (-3.0);
SReal trunk_Q2_3d_AE<13>::C =  -8.0 * 2.0;
SReal trunk_Q2_3d_AE<14>::C =  12.0 * 2.0;
SReal trunk_Q2_3d_AE<15>::C =  -4.0 * 2.0;
SReal trunk_Q2_3d_AE<16>::C =   8.0 * 2.0;
SReal trunk_Q2_3d_AE<17>::C = -12.0 * 2.0;
SReal trunk_Q2_3d_AE<18>::C =   4.0 * 2.0;

SReal trunk_Q2_3d_AF<1>::C = -8.0;
SReal trunk_Q2_3d_AF<2>::C =  8.0;
SReal trunk_Q2_3d_AF<3>::C =  4.0;
SReal trunk_Q2_3d_AF<4>::C = -4.0;
SReal trunk_Q2_3d_AF<5>::C = -8.0 * (-3.0);
SReal trunk_Q2_3d_AF<6>::C =  8.0 * (-3.0);
SReal trunk_Q2_3d_AF<7>::C =  4.0 * (-3.0);
SReal trunk_Q2_3d_AF<8>::C = -4.0 * (-3.0);
SReal trunk_Q2_3d_AF<9>::C = -8.0 * 2.0;
SReal trunk_Q2_3d_AF<10>::C =  8.0 * 2.0;
SReal trunk_Q2_3d_AF<11>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_AF<12>::C = -4.0 * 2.0;

SReal trunk_Q2_3d_AG<1>::C = -8.0;
SReal trunk_Q2_3d_AG<2>::C =  4.0;
SReal trunk_Q2_3d_AG<3>::C =  8.0;
SReal trunk_Q2_3d_AG<4>::C = -4.0;
SReal trunk_Q2_3d_AG<5>::C = -8.0 * (-3.0);
SReal trunk_Q2_3d_AG<6>::C =  4.0 * (-3.0);
SReal trunk_Q2_3d_AG<7>::C =  8.0 * (-3.0);
SReal trunk_Q2_3d_AG<8>::C = -4.0 * (-3.0);
SReal trunk_Q2_3d_AG<9>::C = -8.0 * 2.0;
SReal trunk_Q2_3d_AG<10>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_AG<11>::C =  8.0 * 2.0;
SReal trunk_Q2_3d_AG<12>::C = -4.0 * 2.0;
 
SReal trunk_Q2_3d_AH<1>::C =  -8.0;
SReal trunk_Q2_3d_AH<2>::C =   8.0;
SReal trunk_Q2_3d_AH<3>::C =  12.0;
SReal trunk_Q2_3d_AH<4>::C = -12.0;
SReal trunk_Q2_3d_AH<5>::C =  -4.0;
SReal trunk_Q2_3d_AH<6>::C =   4.0;
SReal trunk_Q2_3d_AH<7>::C =  -8.0 * (-3.0);
SReal trunk_Q2_3d_AH<8>::C =   8.0 * (-3.0);
SReal trunk_Q2_3d_AH<9>::C =  12.0 * (-3.0);
SReal trunk_Q2_3d_AH<10>::C = -12.0 * (-3.0);
SReal trunk_Q2_3d_AH<11>::C =  -4.0 * (-3.0);
SReal trunk_Q2_3d_AH<12>::C =   4.0 * (-3.0);
SReal trunk_Q2_3d_AH<13>::C =  -8.0 * 2.0;
SReal trunk_Q2_3d_AH<14>::C =   8.0 * 2.0;
SReal trunk_Q2_3d_AH<15>::C =  12.0 * 2.0;
SReal trunk_Q2_3d_AH<16>::C = -12.0 * 2.0;
SReal trunk_Q2_3d_AH<17>::C =  -4.0 * 2.0;
SReal trunk_Q2_3d_AH<18>::C =   4.0 * 2.0;

SReal trunk_Q2_3d_AI<1>::C =  16.0;
SReal trunk_Q2_3d_AI<2>::C = -16.0;
SReal trunk_Q2_3d_AI<3>::C = -16.0;
SReal trunk_Q2_3d_AI<4>::C =  16.0;
SReal trunk_Q2_3d_AI<5>::C =  16.0 * (-3.0);
SReal trunk_Q2_3d_AI<6>::C = -16.0 * (-3.0);
SReal trunk_Q2_3d_AI<7>::C = -16.0 * (-3.0);
SReal trunk_Q2_3d_AI<8>::C =  16.0 * (-3.0);
SReal trunk_Q2_3d_AI<9>::C =  16.0 * 2.0;
SReal trunk_Q2_3d_AI<10>::C = -16.0 * 2.0;
SReal trunk_Q2_3d_AI<11>::C = -16.0 * 2.0;
SReal trunk_Q2_3d_AI<12>::C =  16.0 * 2.0;

SReal trunk_Q2_3d_BA<1>::C =  4.0 * 4.0;
SReal trunk_Q2_3d_BA<2>::C = -6.0 * 4.0;
SReal trunk_Q2_3d_BA<3>::C =  2.0 * 4.0;
SReal trunk_Q2_3d_BA<4>::C = -6.0 * 4.0;
SReal trunk_Q2_3d_BA<5>::C =  9.0 * 4.0;
SReal trunk_Q2_3d_BA<6>::C = -3.0 * 4.0;
SReal trunk_Q2_3d_BA<7>::C =  2.0 * 4.0;
SReal trunk_Q2_3d_BA<8>::C = -3.0 * 4.0;
SReal trunk_Q2_3d_BA<9>::C =  1.0 * 4.0;
SReal trunk_Q2_3d_BA<10>::C =  4.0 * (-4.0);
SReal trunk_Q2_3d_BA<11>::C = -6.0 * (-4.0);
SReal trunk_Q2_3d_BA<12>::C =  2.0 * (-4.0);
SReal trunk_Q2_3d_BA<13>::C = -6.0 * (-4.0);
SReal trunk_Q2_3d_BA<14>::C =  9.0 * (-4.0);
SReal trunk_Q2_3d_BA<15>::C = -3.0 * (-4.0);
SReal trunk_Q2_3d_BA<16>::C =  2.0 * (-4.0);
SReal trunk_Q2_3d_BA<17>::C = -3.0 * (-4.0);
SReal trunk_Q2_3d_BA<18>::C =  1.0 * (-4.0);

SReal trunk_Q2_3d_BB<1>::C =  4.0 * 4.0;
SReal trunk_Q2_3d_BB<2>::C = -6.0 * 4.0;
SReal trunk_Q2_3d_BB<3>::C =  2.0 * 4.0;
SReal trunk_Q2_3d_BB<4>::C = -2.0 * 4.0;
SReal trunk_Q2_3d_BB<5>::C =  3.0 * 4.0;
SReal trunk_Q2_3d_BB<6>::C = -1.0 * 4.0;
SReal trunk_Q2_3d_BB<7>::C =  4.0 * (-4.0);
SReal trunk_Q2_3d_BB<8>::C = -6.0 * (-4.0);
SReal trunk_Q2_3d_BB<9>::C =  2.0 * (-4.0);
SReal trunk_Q2_3d_BB<10>::C = -2.0 * (-4.0);
SReal trunk_Q2_3d_BB<11>::C =  3.0 * (-4.0);
SReal trunk_Q2_3d_BB<12>::C = -1.0 * (-4.0);

SReal trunk_Q2_3d_BC<1>::C =  4.0 * 4.0;
SReal trunk_Q2_3d_BC<2>::C = -2.0 * 4.0;
SReal trunk_Q2_3d_BC<3>::C = -2.0 * 4.0;
SReal trunk_Q2_3d_BC<4>::C =  1.0 * 4.0;
SReal trunk_Q2_3d_BC<5>::C =  4.0 * (-4.0);
SReal trunk_Q2_3d_BC<6>::C = -2.0 * (-4.0);
SReal trunk_Q2_3d_BC<7>::C = -2.0 * (-4.0);
SReal trunk_Q2_3d_BC<8>::C =  1.0 * (-4.0);

SReal trunk_Q2_3d_BD<1>::C =  4.0 * 4.0;
SReal trunk_Q2_3d_BD<2>::C = -2.0 * 4.0;
SReal trunk_Q2_3d_BD<3>::C = -6.0 * 4.0;
SReal trunk_Q2_3d_BD<4>::C =  3.0 * 4.0;
SReal trunk_Q2_3d_BD<5>::C =  2.0 * 4.0;
SReal trunk_Q2_3d_BD<6>::C = -1.0 * 4.0;
SReal trunk_Q2_3d_BD<7>::C =  4.0 * (-4.0);
SReal trunk_Q2_3d_BD<8>::C = -2.0 * (-4.0);
SReal trunk_Q2_3d_BD<9>::C = -6.0 * (-4.0);
SReal trunk_Q2_3d_BD<10>::C =  3.0 * (-4.0);
SReal trunk_Q2_3d_BD<11>::C =  2.0 * (-4.0);
SReal trunk_Q2_3d_BD<12>::C = -1.0 * (-4.0);

SReal trunk_Q2_3d_BE<1>::C =  -8.0 * 4.0;
SReal trunk_Q2_3d_BE<2>::C =  12.0 * 4.0;
SReal trunk_Q2_3d_BE<3>::C =  -4.0 * 4.0;
SReal trunk_Q2_3d_BE<4>::C =   8.0 * 4.0;
SReal trunk_Q2_3d_BE<5>::C = -12.0 * 4.0;
SReal trunk_Q2_3d_BE<6>::C =   4.0 * 4.0;
SReal trunk_Q2_3d_BE<7>::C =  -8.0 * (-4.0);
SReal trunk_Q2_3d_BE<8>::C =  12.0 * (-4.0);
SReal trunk_Q2_3d_BE<9>::C =  -4.0 * (-4.0);
SReal trunk_Q2_3d_BE<10>::C =   8.0 * (-4.0);
SReal trunk_Q2_3d_BE<11>::C = -12.0 * (-4.0);
SReal trunk_Q2_3d_BE<12>::C =   4.0 * (-4.0);

SReal trunk_Q2_3d_BF<1>::C = -8.0 * 4.0;
SReal trunk_Q2_3d_BF<2>::C =  8.0 * 4.0;
SReal trunk_Q2_3d_BF<3>::C =  4.0 * 4.0;
SReal trunk_Q2_3d_BF<4>::C = -4.0 * 4.0;
SReal trunk_Q2_3d_BF<5>::C = -8.0 * (-4.0);
SReal trunk_Q2_3d_BF<6>::C =  8.0 * (-4.0);
SReal trunk_Q2_3d_BF<7>::C =  4.0 * (-4.0);
SReal trunk_Q2_3d_BF<8>::C = -4.0 * (-4.0);

SReal trunk_Q2_3d_BG<1>::C = -8.0 * 4.0;
SReal trunk_Q2_3d_BG<2>::C =  4.0 * 4.0;
SReal trunk_Q2_3d_BG<3>::C =  8.0 * 4.0;
SReal trunk_Q2_3d_BG<4>::C = -4.0 * 4.0;
SReal trunk_Q2_3d_BG<5>::C = -8.0 * (-4.0);
SReal trunk_Q2_3d_BG<6>::C =  4.0 * (-4.0);
SReal trunk_Q2_3d_BG<7>::C =  8.0 * (-4.0);
SReal trunk_Q2_3d_BG<8>::C = -4.0 * (-4.0);

SReal trunk_Q2_3d_BH<1>::C =  -8.0 * 4.0;
SReal trunk_Q2_3d_BH<2>::C =   8.0 * 4.0;
SReal trunk_Q2_3d_BH<3>::C =  12.0 * 4.0;
SReal trunk_Q2_3d_BH<4>::C = -12.0 * 4.0;
SReal trunk_Q2_3d_BH<5>::C =  -4.0 * 4.0;
SReal trunk_Q2_3d_BH<6>::C =   4.0 * 4.0;
SReal trunk_Q2_3d_BH<7>::C =  -8.0 * (-4.0);
SReal trunk_Q2_3d_BH<8>::C =   8.0 * (-4.0);
SReal trunk_Q2_3d_BH<9>::C =  12.0 * (-4.0);
SReal trunk_Q2_3d_BH<10>::C = -12.0 * (-4.0);
SReal trunk_Q2_3d_BH<11>::C =  -4.0 * (-4.0);
SReal trunk_Q2_3d_BH<12>::C =   4.0 * (-4.0);

SReal trunk_Q2_3d_BI<1>::C =  16.0 * 4.0;
SReal trunk_Q2_3d_BI<2>::C = -16.0 * 4.0;
SReal trunk_Q2_3d_BI<3>::C = -16.0 * 4.0;
SReal trunk_Q2_3d_BI<4>::C =  16.0 * 4.0;
SReal trunk_Q2_3d_BI<5>::C =  16.0 * (-4.0);
SReal trunk_Q2_3d_BI<6>::C = -16.0 * (-4.0);
SReal trunk_Q2_3d_BI<7>::C = -16.0 * (-4.0);
SReal trunk_Q2_3d_BI<8>::C =  16.0 * (-4.0);

SReal trunk_Q2_3d_CA<1>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_CA<2>::C = -6.0 * 2.0;
SReal trunk_Q2_3d_CA<3>::C =  2.0 * 2.0;
SReal trunk_Q2_3d_CA<4>::C = -6.0 * 2.0;
SReal trunk_Q2_3d_CA<5>::C =  9.0 * 2.0;
SReal trunk_Q2_3d_CA<6>::C = -3.0 * 2.0;
SReal trunk_Q2_3d_CA<7>::C =  2.0 * 2.0;
SReal trunk_Q2_3d_CA<8>::C = -3.0 * 2.0;
SReal trunk_Q2_3d_CA<9>::C =  1.0 * 2.0;
SReal trunk_Q2_3d_CA<10>::C =  4.0 * (-1.0);
SReal trunk_Q2_3d_CA<11>::C = -6.0 * (-1.0);
SReal trunk_Q2_3d_CA<12>::C =  2.0 * (-1.0);
SReal trunk_Q2_3d_CA<13>::C = -6.0 * (-1.0);
SReal trunk_Q2_3d_CA<14>::C =  9.0 * (-1.0);
SReal trunk_Q2_3d_CA<15>::C = -3.0 * (-1.0);
SReal trunk_Q2_3d_CA<16>::C =  2.0 * (-1.0);
SReal trunk_Q2_3d_CA<17>::C = -3.0 * (-1.0);
SReal trunk_Q2_3d_CA<18>::C =  1.0 * (-1.0);

SReal trunk_Q2_3d_CB<1>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_CB<2>::C = -6.0 * 2.0;
SReal trunk_Q2_3d_CB<3>::C =  2.0 * 2.0;
SReal trunk_Q2_3d_CB<4>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_CB<5>::C =  3.0 * 2.0;
SReal trunk_Q2_3d_CB<6>::C = -1.0 * 2.0;
SReal trunk_Q2_3d_CB<7>::C =  4.0 * (-1.0);
SReal trunk_Q2_3d_CB<8>::C = -6.0 * (-1.0);
SReal trunk_Q2_3d_CB<9>::C =  2.0 * (-1.0);
SReal trunk_Q2_3d_CB<10>::C = -2.0 * (-1.0);
SReal trunk_Q2_3d_CB<11>::C =  3.0 * (-1.0);
SReal trunk_Q2_3d_CB<12>::C = -1.0 * (-1.0);

SReal trunk_Q2_3d_CC<1>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_CC<2>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_CC<3>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_CC<4>::C =  1.0 * 2.0;
SReal trunk_Q2_3d_CC<5>::C =  4.0 * (-1.0);
SReal trunk_Q2_3d_CC<6>::C = -2.0 * (-1.0);
SReal trunk_Q2_3d_CC<7>::C = -2.0 * (-1.0);
SReal trunk_Q2_3d_CC<8>::C =  1.0 * (-1.0);

SReal trunk_Q2_3d_CD<1>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_CD<2>::C = -2.0 * 2.0;
SReal trunk_Q2_3d_CD<3>::C = -6.0 * 2.0;
SReal trunk_Q2_3d_CD<4>::C =  3.0 * 2.0;
SReal trunk_Q2_3d_CD<5>::C =  2.0 * 2.0;
SReal trunk_Q2_3d_CD<6>::C = -1.0 * 2.0;
SReal trunk_Q2_3d_CD<7>::C =  4.0 * (-1.0);
SReal trunk_Q2_3d_CD<8>::C = -2.0 * (-1.0);
SReal trunk_Q2_3d_CD<9>::C = -6.0 * (-1.0);
SReal trunk_Q2_3d_CD<10>::C =  3.0 * (-1.0);
SReal trunk_Q2_3d_CD<11>::C =  2.0 * (-1.0);
SReal trunk_Q2_3d_CD<12>::C = -1.0 * (-1.0);

SReal trunk_Q2_3d_CE<1>::C =  -8.0 * 2.0;
SReal trunk_Q2_3d_CE<2>::C =  12.0 * 2.0;
SReal trunk_Q2_3d_CE<3>::C =  -4.0 * 2.0;
SReal trunk_Q2_3d_CE<4>::C =   8.0 * 2.0;
SReal trunk_Q2_3d_CE<5>::C = -12.0 * 2.0;
SReal trunk_Q2_3d_CE<6>::C =   4.0 * 2.0;
SReal trunk_Q2_3d_CE<7>::C =  -8.0 * (-1.0);
SReal trunk_Q2_3d_CE<8>::C =  12.0 * (-1.0);
SReal trunk_Q2_3d_CE<9>::C =  -4.0 * (-1.0);
SReal trunk_Q2_3d_CE<10>::C =   8.0 * (-1.0);
SReal trunk_Q2_3d_CE<11>::C = -12.0 * (-1.0);
SReal trunk_Q2_3d_CE<12>::C =   4.0 * (-1.0);

SReal trunk_Q2_3d_CF<1>::C = -8.0 * 2.0;
SReal trunk_Q2_3d_CF<2>::C =  8.0 * 2.0;
SReal trunk_Q2_3d_CF<3>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_CF<4>::C = -4.0 * 2.0;
SReal trunk_Q2_3d_CF<5>::C = -8.0 * (-1.0);
SReal trunk_Q2_3d_CF<6>::C =  8.0 * (-1.0);
SReal trunk_Q2_3d_CF<7>::C =  4.0 * (-1.0);
SReal trunk_Q2_3d_CF<8>::C = -4.0 * (-1.0);

SReal trunk_Q2_3d_CG<1>::C = -8.0 * 2.0;
SReal trunk_Q2_3d_CG<2>::C =  4.0 * 2.0;
SReal trunk_Q2_3d_CG<3>::C =  8.0 * 2.0;
SReal trunk_Q2_3d_CG<4>::C = -4.0 * 2.0;
SReal trunk_Q2_3d_CG<5>::C = -8.0 * (-1.0);
SReal trunk_Q2_3d_CG<6>::C =  4.0 * (-1.0);
SReal trunk_Q2_3d_CG<7>::C =  8.0 * (-1.0);
SReal trunk_Q2_3d_CG<8>::C = -4.0 * (-1.0);

SReal trunk_Q2_3d_CH<1>::C =  -8.0 * 2.0;
SReal trunk_Q2_3d_CH<2>::C =   8.0 * 2.0;
SReal trunk_Q2_3d_CH<3>::C =  12.0 * 2.0;
SReal trunk_Q2_3d_CH<4>::C = -12.0 * 2.0;
SReal trunk_Q2_3d_CH<5>::C =  -4.0 * 2.0;
SReal trunk_Q2_3d_CH<6>::C =   4.0 * 2.0;
SReal trunk_Q2_3d_CH<7>::C =  -8.0 * (-1.0);
SReal trunk_Q2_3d_CH<8>::C =   8.0 * (-1.0);
SReal trunk_Q2_3d_CH<9>::C =  12.0 * (-1.0);
SReal trunk_Q2_3d_CH<10>::C = -12.0 * (-1.0);
SReal trunk_Q2_3d_CH<11>::C =  -4.0 * (-1.0);
SReal trunk_Q2_3d_CH<12>::C =   4.0 * (-1.0);

SReal trunk_Q2_3d_CI<1>::C =  16.0 * 2.0;
SReal trunk_Q2_3d_CI<2>::C = -16.0 * 2.0;
SReal trunk_Q2_3d_CI<3>::C = -16.0 * 2.0;
SReal trunk_Q2_3d_CI<4>::C =  16.0 * 2.0;
SReal trunk_Q2_3d_CI<5>::C =  16.0 * (-1.0);
SReal trunk_Q2_3d_CI<6>::C = -16.0 * (-1.0);
SReal trunk_Q2_3d_CI<7>::C = -16.0 * (-1.0);
SReal trunk_Q2_3d_CI<8>::C =  16.0 * (-1.0);



/*! Polynomial D0 - 3D */
SReal trunk_D0_3d_A<1>::C =  1.0;

SReal trunk_D0_3d_B<1>::C =  1.0;

SReal trunk_D0_3d_C<1>::C =  1.0;
      
SReal trunk_D0_3d_D<1>::C =  1.0;
      

/*! Polynomial P1BP1 - 3D */
SReal trunk_P1BP1_3d_A<1>::C =  1.0;
SReal trunk_P1BP1_3d_A<2>::C = -1.0;
SReal trunk_P1BP1_3d_A<3>::C = -1.0;
SReal trunk_P1BP1_3d_A<4>::C = -1.0;

SReal trunk_P1BP1_3d_B<1>::C =  1.0;

SReal trunk_P1BP1_3d_C<1>::C =  1.0;

SReal trunk_P1BP1_3d_D<1>::C =  1.0;

SReal trunk_P1BP1_3d_E<1>::C =  4.0;
SReal trunk_P1BP1_3d_E<2>::C = -4.0;
SReal trunk_P1BP1_3d_E<3>::C = -4.0;
SReal trunk_P1BP1_3d_E<4>::C = -4.0;

SReal trunk_P1BP1_3d_F<1>::C =  4.0;
      
SReal trunk_P1BP1_3d_G<1>::C =  4.0;
      
SReal trunk_P1BP1_3d_H<1>::C =  4.0;


/*! Polynomial RT0 - 3D */
SReal trunk_RT0_3d_Ax<1>::C =  1.0;

SReal trunk_RT0_3d_Ay<1>::C =  1.0;

SReal trunk_RT0_3d_Az<1>::C = -1.0;
SReal trunk_RT0_3d_Az<2>::C =  1.0; 

SReal trunk_RT0_3d_Bx<1>::C =  1.0;

SReal trunk_RT0_3d_By<1>::C = -1.0;
SReal trunk_RT0_3d_By<2>::C =  1.0;
      
SReal trunk_RT0_3d_Bz<1>::C =  1.0;  
      
SReal trunk_RT0_3d_Cx<1>::C =  1.0; 
      
SReal trunk_RT0_3d_Cy<1>::C =  1.0;
      
SReal trunk_RT0_3d_Cz<1>::C =  1.0;
      
SReal trunk_RT0_3d_Dx<1>::C = -1.0;
SReal trunk_RT0_3d_Dx<2>::C =  1.0;     

SReal trunk_RT0_3d_Dy<1>::C =  1.0;     
 
SReal trunk_RT0_3d_Dz<1>::C =  1.0;


/*! Polynomial W1 - 2D */
SReal trunk_W1_3d_A<1>::C =  1.0;
SReal trunk_W1_3d_A<2>::C = -1.0;
SReal trunk_W1_3d_A<3>::C = -1.0;
SReal trunk_W1_3d_A<4>::C = -1.0;
SReal trunk_W1_3d_A<5>::C =  1.0;
SReal trunk_W1_3d_A<6>::C =  1.0;

SReal trunk_W1_3d_B<1>::C =  1.0;
SReal trunk_W1_3d_B<2>::C = -1.0;

SReal trunk_W1_3d_C<1>::C =  1.0;
SReal trunk_W1_3d_C<2>::C = -1.0;

SReal trunk_W1_3d_D<1>::C =  1.0;
SReal trunk_W1_3d_D<2>::C = -1.0;
SReal trunk_W1_3d_D<3>::C = -1.0;

SReal trunk_W1_3d_E<1>::C =  1.0;

SReal trunk_W1_3d_F<1>::C =  1.0;

/*! Polynomial CR1 - 3D */
SReal trunk_CR1_3d_A<1>::C = -2.0;
SReal trunk_CR1_3d_A<2>::C =  3.0;
SReal trunk_CR1_3d_A<3>::C =  3.0;
SReal trunk_CR1_3d_A<4>::C =  3.0;

SReal trunk_CR1_3d_B<1>::C =  1.0;
SReal trunk_CR1_3d_B<2>::C = -3.0;

SReal trunk_CR1_3d_C<1>::C =  1.0;
SReal trunk_CR1_3d_C<2>::C = -3.0;

SReal trunk_CR1_3d_D<1>::C =  1.0;
SReal trunk_CR1_3d_D<2>::C = -3.0;
