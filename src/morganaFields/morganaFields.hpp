/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef MORGANAFIELDS_HPP
#define MORGANAFIELDS_HPP

enum FEClass { feStatic, feDynamic, feDoubleLayer };

enum FEGridType {primal=1, dual=2};

enum FEBaseType{scalar=1, vectorial=2, komplexx=3};

enum FELabel {FE_P0_3d, FE_P1_3d, FE_P2_3d, FE_Q0_3d, FE_Q1_3d, FE_Q2_3d, FE_D0LT_3d, FE_P1B1_3d, FE_RT0LT_3d, FE_RT0LOC_3d, FE_FVFLUX_3d,
              FE_P0_2d, FE_P1_2d, FE_P2_2d, FE_Q0_2d, FE_Q1_2d, FE_Q2_2d, FE_D0LT_2d,             FE_RT0LT_2d, FE_RT0LOC_2d, FE_FVFLUX_2d,
              FE_P0_1d, FE_P1_1d, FE_P2_1d,
              FE_XF_3d, FE_SK_3d, FE_HY_3d, FE_HYB_3d, FE_OSCHYB_3d,
                        FE_SK_2d, FE_HY_2d,
                                  FE_SK_1d,
              FE_CR1_3d,
              FE_OS1P0_3d, FE_OS1P1_3d,
              FE_OS0P0_2d, FE_OS0P1_2d,
              FE_OSCR1_3d,
              FE_DGN0_3d,  FE_DGN1_3d};

enum FEBaseLabel {BL_Pr_3d, BL_Qr_3d, BL_D0LT_3d, BL_P1B1_3d, BL_RT0LT_3d, BL_RT0LOC_3d, BL_FVFLUX_3d,
                  BL_Pr_2d, BL_Qr_2d, BL_D0LT_2d,             BL_RT0LT_2d, BL_RT0LOC_2d, BL_FVFLUX_2d,
                  BL_Pr_1d,
                  BL_XF_3d, BL_SK_3d, BL_HY_3d, BL_HYB_3d, BL_OSCHYB_3d,
                            BL_SK_2d, BL_HY_2d,
                                      BL_SK_1d,
                  BL_CRr_3d,
                  BL_OSwPr_3d, BL_OSwPr_2d, BL_OSCRr_3d,
                  BL_DGNr_3d};

enum FECardLabel {CR_HY_2d, CR_HY_3d, CR_HYB_3d, CR_OS_3d, CR_OSCR_3d, CR_OSHYB_3d, CR_RT0LT_2d, CR_RT0LT_3d, CR_SK_1d, CR_SK_2d, CR_SK_3d, CR_XF_3d};

#endif
