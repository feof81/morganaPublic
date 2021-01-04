/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "time.h"
#include <cmath>
#include <iostream>

#include "typesInterface.hpp"
#include "polyCards.h"
#include "polyStatic.hpp"

/*! The Q2 basis forms a partition of the unity. Therefore the first and second derivatives should be null. */
int main(int argc, char *argv[])
{
  point3d P(1.0,0.3,0.5);
  static const UInt dx = 0, dy = 1, dz = 1;
  
  
  polyStatic<Q2_3d_AA> polyEval_AA;
  polyStatic<Q2_3d_AB> polyEval_AB;
  polyStatic<Q2_3d_AC> polyEval_AC;
  polyStatic<Q2_3d_AD> polyEval_AD;
  polyStatic<Q2_3d_AE> polyEval_AE;
  polyStatic<Q2_3d_AF> polyEval_AF;
  polyStatic<Q2_3d_AG> polyEval_AG;
  polyStatic<Q2_3d_AH> polyEval_AH;
  polyStatic<Q2_3d_AI> polyEval_AI;
  
  polyStatic<Q2_3d_BA> polyEval_BA;
  polyStatic<Q2_3d_BB> polyEval_BB;
  polyStatic<Q2_3d_BC> polyEval_BC;
  polyStatic<Q2_3d_BD> polyEval_BD;
  polyStatic<Q2_3d_BE> polyEval_BE;
  polyStatic<Q2_3d_BF> polyEval_BF;
  polyStatic<Q2_3d_BG> polyEval_BG;
  polyStatic<Q2_3d_BH> polyEval_BH;
  polyStatic<Q2_3d_BI> polyEval_BI;
  
  polyStatic<Q2_3d_CA> polyEval_CA;
  polyStatic<Q2_3d_CB> polyEval_CB;
  polyStatic<Q2_3d_CC> polyEval_CC;
  polyStatic<Q2_3d_CD> polyEval_CD;
  polyStatic<Q2_3d_CE> polyEval_CE;
  polyStatic<Q2_3d_CF> polyEval_CF;
  polyStatic<Q2_3d_CG> polyEval_CG;
  polyStatic<Q2_3d_CH> polyEval_CH;
  polyStatic<Q2_3d_CI> polyEval_CI;
  
  
  Real u =
  polyEval_AA.evaluateStatic(P) + polyEval_AB.evaluateStatic(P) + polyEval_AC.evaluateStatic(P) +
  polyEval_AD.evaluateStatic(P) + polyEval_AE.evaluateStatic(P) + polyEval_AF.evaluateStatic(P) +
  polyEval_AG.evaluateStatic(P) + polyEval_AH.evaluateStatic(P) + polyEval_AI.evaluateStatic(P) +
  
  polyEval_BA.evaluateStatic(P) + polyEval_BB.evaluateStatic(P) + polyEval_BC.evaluateStatic(P) +
  polyEval_BD.evaluateStatic(P) + polyEval_BE.evaluateStatic(P) + polyEval_BF.evaluateStatic(P) +
  polyEval_BG.evaluateStatic(P) + polyEval_BH.evaluateStatic(P) + polyEval_BI.evaluateStatic(P) +
  
  polyEval_CA.evaluateStatic(P) + polyEval_CB.evaluateStatic(P) + polyEval_CC.evaluateStatic(P) +
  polyEval_CD.evaluateStatic(P) + polyEval_CE.evaluateStatic(P) + polyEval_CF.evaluateStatic(P) +
  polyEval_CG.evaluateStatic(P) + polyEval_CH.evaluateStatic(P) + polyEval_CI.evaluateStatic(P);
  
  cout << "eval: " << u << endl;
  
  
  point3d G_AA = polyEval_AA.evaluateGradient(P);
  point3d G_AB = polyEval_AB.evaluateGradient(P);
  point3d G_AC = polyEval_AC.evaluateGradient(P);
  point3d G_AD = polyEval_AD.evaluateGradient(P);
  point3d G_AE = polyEval_AE.evaluateGradient(P);
  point3d G_AF = polyEval_AF.evaluateGradient(P);
  point3d G_AG = polyEval_AG.evaluateGradient(P);
  point3d G_AH = polyEval_AH.evaluateGradient(P);
  point3d G_AI = polyEval_AI.evaluateGradient(P);
  
  point3d G_BA = polyEval_BA.evaluateGradient(P);
  point3d G_BB = polyEval_BB.evaluateGradient(P);
  point3d G_BC = polyEval_BC.evaluateGradient(P);
  point3d G_BD = polyEval_BD.evaluateGradient(P);
  point3d G_BE = polyEval_BE.evaluateGradient(P);
  point3d G_BF = polyEval_BF.evaluateGradient(P);
  point3d G_BG = polyEval_BG.evaluateGradient(P);
  point3d G_BH = polyEval_BH.evaluateGradient(P);
  point3d G_BI = polyEval_BI.evaluateGradient(P);
  
  point3d G_CA = polyEval_CA.evaluateGradient(P);
  point3d G_CB = polyEval_CB.evaluateGradient(P);
  point3d G_CC = polyEval_CC.evaluateGradient(P);
  point3d G_CD = polyEval_CD.evaluateGradient(P);
  point3d G_CE = polyEval_CE.evaluateGradient(P);
  point3d G_CF = polyEval_CF.evaluateGradient(P);
  point3d G_CG = polyEval_CG.evaluateGradient(P);
  point3d G_CH = polyEval_CH.evaluateGradient(P);
  point3d G_CI = polyEval_CI.evaluateGradient(P);
  
  
  cout << "G_AA " << G_AA;
  cout << "G_AB " << G_AB;
  cout << "G_AC " << G_AC;
  cout << "G_AD " << G_AD;
  cout << "G_AE " << G_AE;
  cout << "G_AF " << G_AF;
  cout << "G_AG " << G_AG;
  cout << "G_AH " << G_AH;
  cout << "G_AI " << G_AI;
  
  cout << "G_BA " << G_BA;
  cout << "G_BB " << G_BB;
  cout << "G_BC " << G_BC;
  cout << "G_BD " << G_BD;
  cout << "G_BE " << G_BE;
  cout << "G_BF " << G_BF;
  cout << "G_BG " << G_BG;
  cout << "G_BH " << G_BH;
  cout << "G_BI " << G_BI;
  
  cout << "G_CA " << G_CA;
  cout << "G_CB " << G_CB;
  cout << "G_CC " << G_CC;
  cout << "G_CD " << G_CD;
  cout << "G_CE " << G_CE;
  cout << "G_CF " << G_CF;
  cout << "G_CG " << G_CG;
  cout << "G_CH " << G_CH;
  cout << "G_CI " << G_CI;
  
  point3d GA = G_AA + G_AB + G_AC + G_AD + G_AE + G_AF + G_AG + G_AH + G_AI;
  point3d GB = G_BA + G_BB + G_BC + G_BD + G_BE + G_BF + G_BG + G_BH + G_BI;
  point3d GC = G_CA + G_CB + G_CC + G_CD + G_CE + G_CF + G_CG + G_CH + G_CI;
  
  cout << "Final gradient: " << GA + GB + GC << endl;
  
  
  Real H_AA = polyEval_AA.evaluateDerivative<dx,dy,dz>(P);
  Real H_AB = polyEval_AB.evaluateDerivative<dx,dy,dz>(P);
  Real H_AC = polyEval_AC.evaluateDerivative<dx,dy,dz>(P);
  Real H_AD = polyEval_AD.evaluateDerivative<dx,dy,dz>(P);
  Real H_AE = polyEval_AE.evaluateDerivative<dx,dy,dz>(P);
  Real H_AF = polyEval_AF.evaluateDerivative<dx,dy,dz>(P);
  Real H_AG = polyEval_AG.evaluateDerivative<dx,dy,dz>(P);
  Real H_AH = polyEval_AH.evaluateDerivative<dx,dy,dz>(P);
  Real H_AI = polyEval_AI.evaluateDerivative<dx,dy,dz>(P);
  
  Real H_BA = polyEval_BA.evaluateDerivative<dx,dy,dz>(P);
  Real H_BB = polyEval_BB.evaluateDerivative<dx,dy,dz>(P);
  Real H_BC = polyEval_BC.evaluateDerivative<dx,dy,dz>(P);
  Real H_BD = polyEval_BD.evaluateDerivative<dx,dy,dz>(P);
  Real H_BE = polyEval_BE.evaluateDerivative<dx,dy,dz>(P);
  Real H_BF = polyEval_BF.evaluateDerivative<dx,dy,dz>(P);
  Real H_BG = polyEval_BG.evaluateDerivative<dx,dy,dz>(P);
  Real H_BH = polyEval_BH.evaluateDerivative<dx,dy,dz>(P);
  Real H_BI = polyEval_BI.evaluateDerivative<dx,dy,dz>(P);
  
  Real H_CA = polyEval_CA.evaluateDerivative<dx,dy,dz>(P);
  Real H_CB = polyEval_CB.evaluateDerivative<dx,dy,dz>(P);
  Real H_CC = polyEval_CC.evaluateDerivative<dx,dy,dz>(P);
  Real H_CD = polyEval_CD.evaluateDerivative<dx,dy,dz>(P);
  Real H_CE = polyEval_CE.evaluateDerivative<dx,dy,dz>(P);
  Real H_CF = polyEval_CF.evaluateDerivative<dx,dy,dz>(P);
  Real H_CG = polyEval_CG.evaluateDerivative<dx,dy,dz>(P);
  Real H_CH = polyEval_CH.evaluateDerivative<dx,dy,dz>(P);
  Real H_CI = polyEval_CI.evaluateDerivative<dx,dy,dz>(P);
  
  cout << "H_AA " << H_AA << endl;
  cout << "H_AB " << H_AB << endl;
  cout << "H_AC " << H_AC << endl;
  cout << "H_AD " << H_AD << endl;
  cout << "H_AE " << H_AE << endl;
  cout << "H_AF " << H_AF << endl;
  cout << "H_AG " << H_AG << endl;
  cout << "H_AH " << H_AH << endl;
  cout << "H_AI " << H_AI << endl;
  
  cout << "H_BA " << H_BA << endl;
  cout << "H_BB " << H_BB << endl;
  cout << "H_BC " << H_BC << endl;
  cout << "H_BD " << H_BD << endl;
  cout << "H_BE " << H_BE << endl;
  cout << "H_BF " << H_BF << endl;
  cout << "H_BG " << H_BG << endl;
  cout << "H_BH " << H_BH << endl;
  cout << "H_BI " << H_BI << endl;
  
  cout << "H_CA " << H_CA << endl;
  cout << "H_CB " << H_CB << endl;
  cout << "H_CC " << H_CC << endl;
  cout << "H_CD " << H_CD << endl;
  cout << "H_CE " << H_CE << endl;
  cout << "H_CF " << H_CF << endl;
  cout << "H_CG " << H_CG << endl;
  cout << "H_CH " << H_CH << endl;
  cout << "H_CI " << H_CI << endl;
  
  cout << "Hessian final: " <<
  H_AA + H_AB + H_AC + H_AD + H_AE + H_AF + H_AG + H_AH + H_AI +
  H_BA + H_BB + H_BC + H_BD + H_BE + H_BF + H_BG + H_BH + H_BI + 
  H_CA + H_CB + H_CC + H_CD + H_CE + H_CF + H_CG + H_CH + H_CI << endl;
}
