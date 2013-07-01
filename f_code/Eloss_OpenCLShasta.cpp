//    This is a program to perform an energy loss for a hydro evolution
//    performed with the Open CL Shasta from Jochen Gerhard.
//    By: Barbara Betz

#include <cmath>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <complex>
#include <vector>

//   Globals

double e0, n0, Bag, hc3, pi2;
double dx, dt, t, vol, cut;
double 

}

templete <typename T_2D>
T_2D **AllocateDynamicArray_2D( int nRows, int nCols)
(
 
     T_2D **dynamicArray_2D;
     
     dynamicArray_2D = new T_2D *[nRows]
     for(int i = 0; i < nRows; i++){
        dynamicArray_2D[i] = new T_2D [nCols]
     }

     return dynamicArray_2D;
 )

templete <typname T_3D>
T_3D ***AllocateDynamicArray_3D( int nRows, int nCols, int nLevs)
(  


     T_3D ***dynamicArray_3D;

     dynamicArray_3D = new T_3D *[nRows]
     for(int i = 0; i < nRows; i++){
       dynamicArray_3D[i] = new T_3D [nCols];
       for(int j = 0; j < nLevs; j++){
         dynamicArray_3D[i][j] = new T_3D [nLevs]         
	   }
     }

     return dynamicArray_3D;
)


int main(){
   
  // define variables

  int array_size = 220;
  double dt, dx, t, vol, cut, lam;
  double const e0, n0, Bag, hc3, pi2;
  int l, flag, iphi, iphimax;

  // define arrays

  double **u0_dynamic = 0.0;
  
  // memory allocated for elements of rows

  u0_dynamic = new double *[array_size];

  // memory allocated for elements of each column.

  for (int i =0; i < array_size; i++){
    u0_dynamic[i] = new double[array_size];
  }

  // free allocated memory

  //  for(int i =0; i< array_size; i++){
  //  delete [] 
    
 
  // define constants

  pi2 = acos(-1.0) * acos(-1.0);
  hc3 = 197.327053;
  e0 = 146.51751415742;
  n0 = 0.15891;
  Bag = 235;

  //  Number to be read off the hydro output.

  // time increment

    dt = 0.16;

  // spacial increment

    dx = 0.2;
  
  // p_pi = Pf - transverse ,omentum considered [GeV]

    p_pi = 7.5;

  // inital time t0 [fm]

    t0 = 0.0

  /*
     Please note:
     This program is meant to calculate the energy loss for a
     temperature-dependent jet-medium coupling.  Therefore, it
     distinguishes betweeen T1 and a T2 al la LIAO & SHURYAK
     [J.Liao and E.Shuryak, Phys. Rev. Lett. 102, 202302 (2009).]
     However, in the case k2 has to be set to a value different
     from kappa, according to LIAO & SHURYAK: k2 = 3*kappa.
     PLEASE NOTE that k2 right now is set to k2=kappa so that
     there is NO temperature dependence in the jet-medium coupling.
  */
  
    // Temperatures are in units of GeV

    T1 = 0.173;
    T2 = 0.113;

    // parameters a1 , b1

    a1 = 0.0;
    b1 = 1.0

    // kappa - kappa = kappa_Q

    kappa = 3.1;

    // kappa 2 - kappa2 = kappa_M
    
    k2 = kappa;

    // flag for RHIC@200 = 1, LHC@2760 = 2

    flag = 1;

    // maximal momenta for which the fragmentation function works properly

    pmaxRHIC = 60.0;
    pmaxLHC = 500.0;
 
    pi = sqrt(pi2);

    // ngr
    ngr = 201;
  
    ifstream file30;
    file30.open("./test.csv");
    
    n1 =1 

      for(int i = 0; i <= ngr - 1; i++){
        for(int j = 0; j <= ngr - 1; j++){

          xgr[i] = -20.0 + dx * (i - 1);
          ygr[j] = -20.0 + dx * (j - 1);

        }
      }

    for(int k = 0; k <= 200 - 1; k++){
      for(int i = 1; i <= ngr - 1; i++){
        for(int j = 1; j <= ngr - 1; j++){

          Tmt[k][i][j] = 0.0;
          ux[k][i][j] = 0.0;
          uy[k][i][j] = 0.0;
          press[k][i][j] = 0.0;
 
        }
      }
    } 

    ofstream file75;
    ofstream file76;
    ofstream file81;
    ofstream file93;
    ofstream file94;

    file75.open("vn_pi.dat");
    file76.open("Rjet_ecc.dat");
    file81.open("Rjet_pi.dat");
    file93.open("RAA_pi_phi.dat");
    file94.open("RAA_in_out.dat");
   
    An[1] = 0.0;
    Bn[1] = 0.0;
    An[2] = 0.0;
    Bn[2] = 0.0;
    An[3] = 0.0;
    Bn[3] = 0.0;
    An[4] = 0.0;
    Bn[4] = 0.0;
    Sume2 = 0.0;
  
    Rx = 0.0;
    Ry = 0.0;
    Sume = 0.0;
   
    for(int i = 1; i <= ngr; i++){
      x = xgr[i];
    
      for(int j = 1; j <= ngr; j++){
        y = ygr[j];

        Rx = Rx + e[j][i] * x;
        Ry = Ry + e[j][i] * y;
        Sume = Sume + e[j][i];

      }
    } 

    Rx = Rx / Sume;
    Ry = Ry / Sume;

    for (int i = 1; i <= ngr; i++){
      x = xgr[i];
      for(int j = 1; j <= ngr; j++){
        y = ygr[j]

	phi = atan( (y - Ry) / (x - Rx)) + pi;  //I still need to convert this to complex form
        An[1] = An[1] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * cos(1.0 * phi); 
        Bn[1] = Bn[1] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * sin(1.0 * phi);
        An[2] = An[2] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * cos(2.0 * phi);
        Bn[2] = Bn[2] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * sin(2.0 * phi);
        An[3] = An[3] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * cos(3.0 * phi);
	Bn[3] = Bn[3] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * sin(3.0 * phi);
        An[4] = An[4] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * cos(4.0 * phi);
        Bn[4] = Bn[4] + e[i][j] * ((x - Rx) * (x - Rx) + (y - Ry) * (y - Ry)) * sin(4.0 * phi);
        Sume2 = Sume2 + e[i][j]*((x - Rx)*(x - Rx) + (y - Ry)*(y - Ry));
      }
    }  

    An[1] = An[1]/Sume2;
    Bn[1] = Bn[1]/Sume2;
    An[2] = An[2]/Sume2;
    Bn[2] = Bn[2]/Sume2;
    An[3] = An[3]/Sume2;
    Bn[3] = Bn[3]/Sume2;
    An[4] = An[4]/Sume2;
    Bn[4] = Bn[4]/Sume2;
    e_n[1] = sqrt(An[1] * An[1] + Bn[1] * Bn[1])
    psin[1] = 1.0 / 1.0 * 

    //
    //////////////////////////////////////////////////////////////////////
    //     RAA_Parton(p_pi) calculates the partonic RAA for quarks
    //                      and gluons (partons) for a momentum p_pi

    RAA_Parton(p_pi, t0, T1, T2, a1, b1, kappa, k2, e, Rjet_g, Rjet_q, TAA, flag);

    //  Starting Integration

    RAA_pi = 0.0;
    dsig_pi3 = 0.0;
    dsig_pi4 = 0.0;
    iphimax = 72;
  
    for(int iphi = 0; iphi <= iphimax - 1; iphi++){
      RAA_pi1[iphi] = 0.0;
      dsig_pi = 0.0;
      dsig_pi2 = 0.0;
      
      if(flag == 1){
	pmax = pmaxRHIC;
      }
      else if(flag == 2){
	pmax = pmaxLHC;
      }

      trapzd(int1, p_pi / pmax, 1.0, dsig_pi1, 4);
      trapzd(int2, p_pi / pmax, 1.0, dsig_pi2, 4);

      if(iphi == 1){
	trapzd(int3, p_pi / pmax, 1.0, dsig_pi3, 4);
	trapzd(int4, p_pi / pmax, 1.0, dsig_pi4, 4);
      }
      RAA_pi1[iphi] = (dsig_pi1 + dsig_pi2) / (dsig_pi3 + dsig_pi4);
      RAA_pi = RAA_pi + RAA_pi1[iphi - 1] * (2.0 * pi) / iphimax / (2.0 * pi);

      file93 >> 5.0 * (iphi - 1) >> "\t" >> RAA_pi1[iphi] >> std::endl;
    }
    
    file81 >> "RAA_pi\t" >> RAA_pi >> "\t time:\t" >> t >> std::endl;

    file93.close();
    file81.close();
    //
    //  vn for pions
    //

    v1_pi1 = 0.0;
    v2_pi1 = 0.0;
    v3_pi1 = 0.0;
    v4_pi1 = 0.0;
    vn_pi2 = 0.0;

    //
    //  For th Glauber average, psin - 0.0;
    //

    psin[1] = 0.0;
    psin[2] = 0.0;
    psin[3] = 0.0;
    psin[4] = 0.0

      for(int iphi = 1; iphi <= iphimax; iphi++){
	phi = 5.0 * (iphi - 1);
	phi = phi / 180.0 * pi;
	v1_pi1 = v1_pi1 + cos(1.0 * (phi - psin[1])) * RAA_pi1(iphi) * (2.0 * pi) / iphimax;
	v2_pi1 = v2_pi1 + cos(2.0 * (phi - psin[2])) * RAA_pi1(iphi) * (2.0 * pi) / iphimax;
	v3_pi1 = v3_pi1 + cos(3.0 * (phi - psin[3])) * RAA_pi1(iphi) * (2.0 * pi) / iphimax;
	v4_pi1 = v4_pi1 + cos(4.0 * (phi - psin[4])) * RAA_pi1(iphi) * (2.0 * pi) / iphimax;
	vn_pi2 = vn_pi2 + RAA_pi1(iphi) * (2.0 * pi) / iphimax;
      }

    vn_pi[1] = v1_pi1 / vn_pi2;
    vn_pi[2] = v2_pi1 / vn_pi2;
    vn_pi[3] = v3_pi1 / vn_pi2;
    vn_pi[4] = vn_pi1 / vn_pi2;

    file75 >> vn_pi[1] >> "\t" >> vn_pi[2] >> "/t" >> vn_pi[3] >> "\t" >> vn_pi[4] >> "\n v2\t" >> vx_pi[2] >> std::endl;

    file94 >> RAA_pi >> "\t" >> RAA_pi * (1.0 - 2.0 * vn_pi[2]) >> "\t" >> RAA_pi * (1.0 + 2.0 * vn_pi[2]) >> "\t" >> "\n RAAout\t" >> RAA_pi * (1.0 - 2.0 * vn_pi[2]) >> "\n RAAin\t" >> RAA_pi * (1.0 + 2.0 * vn_pi[2]) >> std::endl;

    file75.close();
    file94.close();

    return 0;
  
}  
  //////////////////////////////////////////////////////////////////////////
  //                                                                      //
  //       RAA_Parton(p_pi) calculates the partonic RAA for quarks        //
  //                        and gluons (partons) for a momentum p_pi      //
  //////////////////////////////////////////////////////////////////////////
  //       Type declaratios for variables in common-blocks                //
  //////////////////////////////////////////////////////////////////////////                                                                  

  

void RAA(double p_pi, double t0, double T1, double T2, double a1, double b1, double kappa, double k2, double e, double Rjet_e, double Rjet_q double TAA, int flag){ 

  double dx, dt, t, vol, cut;
  double e0, n0, Bag, hc3, pi2
  
  int array_a = 220;  
  int array_b = 75; 
  
  double e[array_a][array_a], p[array_a][array_a][array_a];
  double elab[array_a][array_a][array_a], n[array_a][array_a][array_a];
  double x, y, z, phi, p_pi, t0, Tmt[array_a][array_a][array_a];
  double xjet, yjet, a1, b1, tjet;
  double lineintS, T1, T2, k2;
  double pi, vjet, Pf_g, Pf_q, P0_g, P0_q, x0, y0;
  double Rjet_g11, Rjet_g[array_b], Rjet_q11, Rjet_q[array_b], Rjet12;
  double Rjet_g3[array_a][array_a][array_b], Rjet_q3[array_a][array_a][array_b]
  double Rjetq, Rjetg;
  double xgr[array_a], ygr[array_a];
  double kappa_g, kappa_q, kappa2_g, kapp2_q, kappa, C_g, C_q;
  double lineint, lineint_g, lineint_q, TAA[array_a][array_a];
  int h, k, iphimax, iphi, ijet, jjet, flag;
  int it, itmax;
  external gft200, qft200 gft2760, qft2760, dpig, dpiq, pressure;

  //
  //  Common-blocks 
  //
  pi = sqrt(pi2);
  ofstream file79;
  file79.open("RAA_Parton.dat");

  // Kappas

  C_g = 3.0
  C_q = 4.0 / 3.0;
  kappa_g = kappa * C_g;
  kappa_q = kappa * C_q;

  //  p_pi = Pf - transverse momentum considered

  Pf_g = p_pi;
  Pf_q = p_pi;

  Rjetq = 0.0;
  Rjetg = 0.0
  iphimax = 72;

  for(int iphi = 0; iphi < iphimax - 1; iphi++){
    Rjet_g11 = 0.0;
    Rjet_g[iphi] = 0.0;
    Rjet_q11 = 0.0;
    Rjet_q[iphi] = 0.0;
    Rjet12 = 0.0;
    //  ATTENTION!!  No double counting of 0 degree !
    phi = 5.0 * (iphi - 1);
    
    for(int i = 0; i <= ngr - 1; i++){
      x0 - xgr[i];
      for(int j = 0; j <= ngr - 1; j++){
        y0 = ygr[i];

	// ATTENTION!! cut in energy density -> there are no jets outside the medium
	//             if there's no cut, one has to pay attention that ijet and jjet
	//             do not go below 0 or above 200 which can happen because of the
	//             additional vjet * t * sin / cos.
	vjet = 1.0;
	lineint = 0.0;
	lineint_g = 0.0;
	lineint_q = 0.0;
	lineintS = 0.0;
	if(Tmt[1][i][j] <= T2){
          for(inr h = 1; h <= itmax; h++){
	    tjet = h * dt + t0;
	    xjet = x0 + vjet * tjet * cos(phi / 180 * pi);
	    yjet = y0 + vjet * tjet * sin(phi / 180 * pi);
	    ijet = nint(xjet / dx + 0.5 + ngr * 0.5);
	    jjet = nint(yjet / dx + 0.5 + ngr * 0.5);
	    
	    if(ijet <= ngr && ijet >= 1 && jjet <= ngr && jjet >= 1){
	      if( h == 0 || h == it){
		if(Tmt(h, i, ijet, jet) <= T1){
		  lineintS = lineintS + 0.5 * pow(tjet,b1) * pow(Tmt(h, ijet, jjet),3)
		  }
	      }
	    }
	    
	    if(Tmt(h,ijet,jjet) >= T1){
	      lineint = lineint + 0.5 * pow(tjet,b1) * pow( pow(Tmt(h,ijet,jjet),3) , (b1 - a1 + 2.0) /2.0) * dx;
	    }
	    else if(Tmt(h,ijet, jjet) >= 
	    
	  }
	}
      }
    }




    //end function 

}




// 
////////////////////////////////////////////////////////////////////////
//

double int1(double x){ 

  double p_pi_int;
 
  extern p_pi, RAA_Parton_Phi, iphi, iphimax, t0, a1, b1, T!, T2, k2, kappa;
  extern e, Rjet_g, Rjet_q, TAA, flag;

  p_pi_int = p_pi / x;
  
  RAA_Parton_Phi(p_pi_int, iphi, iphimax, t0, a1, b1, T1, T2, k2, kappa, e, Rjet_g, Rjet_q, TAA, flag);

  if (flag == 1){
    int1 = 1.0 / x * gft200[p_pi_int] * Rjet_g[iphi] * dpi(x , p_pi_int);
  }
  else{ 
    int1 == 1.0 / x * gft2760[p_pi_int] * Rjet_g[iphi] *dpig(x , p_pi_int); 
  }

  return int1;
}

//
/////////////////////////////////////////////////////////////////////////
//

double int2(double x){
  
  double x, p_pi, Rjet_g[75], Rjet_q[75], p_pi_int;
  double e[220][220]; n[220][220]
  double gft200, qft200, gft2760, qft2760, dpig, dpiq;
  double t0, a1,b1,kappa, Tmt[220][220][220];
  double xgr[220], ygr[220];
  double T1, T2, k2;
  int iphi, iphimax, flag, it, itmax, ngr;
  
  //extern gft200, qft200, gft2760, qft2760, dpig, dpiq;

  extern p_pi, Rjet_g, Rjet_q;
  extern e, n, TAA;
  extern t0, a1, b1, kappa, k2, T1, T2;
  extern iphi, iphimax, flag;
  extern Tmt, it, itmax;
  extern ngr;
  extern xgr, ygr; 

  p_pi_int = p_pi / x;

  RAA_Parton_Phi(p_pi_int, iphi, iphimax, t0, a1, b1, T!, T2, k2, kappa, e, Rjet_g, Rjet_q, TAA, flag);

  if(flag == 1){
    int2 = 1.0 / x * qft200(p_pi_int) * Rjet_q(iphi) * dpiq(x , p_pi_int)
  }
  else if(flag == 2){
    int2 = 1.0 / x * qft2760(p_pi_int) * Rjet_q(iphi) * dpiq(x , p_pi_int);
  }

  return int2;
}

double int3( double x ){

  extern p_pi;
  double p_pi_int = p_pi / x;
  extern flag;
  extern gft2760, gft200, dpig; 
  
  if(flag ==1){
    int3 = 1.0 / x * gft200(p_pi_int) * dpig(x , p_pi_int);
  }  
  else if(flag == 2){
    int3 = 1.0 / x * gft2760(p_pi_int) * dpig(x , p_pi_int);
  }

return int3;
}

double int4(double x){

  double p_pi_int;

  extern flag, qft200, qft2760, dpiq, p_pi;

  p_pi_int = p_pi / x;

  if(flag == 1){
    int4 = 1.0 / x * qft200(p_pi_int) * dpiq(x , p_pi_int);
  }  
  else if(flag == 2){
    int4 = 1.0 / x * qft2760(p_pi_int) * dpiq(x , p_pi_int);
  }

return int4

}
  
//////////////////////////////////////////////////////////////////////////////
//          Fragmentation function for gluons at the RHIC                   //
//////////////////////////////////////////////////////////////////////////////

double gft200(double x){

  double gft200;
  
  gft200 =  exp(4.20630030030893 - 1.605930861452689 * pow(log(x),1.1032070685491895) - 1.980192950237263 * pow(log(x),1.1039097100978987) - 2.1719153683060965 * pow(log(x),1.1043001289508012) - 0.0001833341007496504 * pow(log(x),7.4023223396120095));

  return gft200;
}

/////////////////////////////////////////////////////////////////////////////
//          Fragmentation Function for quarks at RHIC                      //
/////////////////////////////////////////////////////////////////////////////

double qft200(double x){

  double qft200;

  qft200 = exp(2.527729186493723 - 2.753620838545098 * pow(log(x),1.0509055095944164) - 2.7897248005749744 * pow(log(x),1.0509624453005448) - 0.00007897570030835403 * pow(log(x),7.8604912704152925));
 
  return qft200;
}

/////////////////////////////////////////////////////////////////////////////
//           Fragmentation Function for gluons at LHC                      //
/////////////////////////////////////////////////////////////////////////////

double gft2760(double x){

  double gft2760;
  
  gft2760 = exp(6.874017911786442 - 0.26596588700706614 * pow(log(x),0.665482588442757) - 2.677705869879314 * pow(log(x),0.8020502598676339) - 2.4735502984532656 * pow(log(x),0.8069542250600515) - 0.36687832133337656 * pow(log(x),2.070179064516989));

  return gft2760;
}

///////////////////////////////////////////////////////////////////////////////
//                 Fragmentation Function for quarks at LHC                  //
///////////////////////////////////////////////////////////////////////////////

double qft2760(double x){

  double qft2760;

  qft2760 = exp(2.990387181868 - 2.0117432145703114 * pow(log(x),1.0384884086567516) - 1.9187151702604879 * pow(log(x),1.039887584824982) - 0.15503714201000543 * pow(log(x),1.0586516925018519) - 0.15384778823017106 * pow(log(x),2.08298497208415173));
  
  return qft2760;
}

//////////////////////////////////////////////////////////////////////////////
//                Fitted pQCD cross-section for gluons                      //
//////////////////////////////////////////////////////////////////////////////

double dpig(double x, double y){
  double logqs, dpig;
  
  logqs = log(log((y * y)/(0.088 * 0.088)) / log(2.0/(0.088 * 0.088)));
  dpig = (6.0451 - 6.61523 * logqs - 1.64978 * logqs * logqs + 2.68223 * logqs * logqs * logqs) * pow(x,-0.71378 + 0.14705 * logqs - 1.08423 * logqs * logqs - 0.43182 * logqs * logqs * logqs) * pow((1.0 - x),(2.92133 + 1.48429 * logqs + 1.32887 * logqs * logqs - 1.78696 * logqs * logqs * logqs) * (1.0 + (0.23086 * logqs - 0.29182 * logqs * logqs) / x));

  return dpig;
}

//////////////////////////////////////////////////////////////////////////////
//                Fitted pQCD cross-section for quarks                      //
//////////////////////////////////////////////////////////////////////////////

double dpiq(double x, double y){

  double logqs;

  logqs = log(log(y * y / (0.88 * 0.88)) / log(2.0 / (0.88 * 0.88)));
 
  dpiq = (0.5461 - 0.22946 * logqs - 0.22594 * logqs * logqs + 0.2119 * logqs * logqs * logqs) * pow(x,-1.46616 - 0.45404 * logqs - 0.12684 * logqs * logqs + 0.27646 * logqs * logqs * logqs) * pow((1.0 - x),(1.01864 + 0.95367 * logqs - 1.09835 * logqs * logqs + 0.74657 * logqs * logqs * logqs) * (1.0 + (0.23086 * logqs - 0..29182 * logqs * logqs) / x);

  return dpiq;
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//   This function is a fit to P(T)/T^4 that has to be manipulated by       //
//   T^4 to lead to P(T).  P(T)/T^4 is based on M.~Laine and Y.~Schroder,   //
//   Phys,\ Rev.\ D {\bf 73} (2006) 085009 [arXiv:hep-ph/0603048] used by   //
//   Romatschke.                                                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

double pressure(double x){

  double pressure;

  if (x < 0.18){
    pressure = (x + (x / 0.18) * (x / 0.18) * (x / 0.18)) * x * x * x * x;
  }
  else{
    pressure = (1.6 + 5.0 * (1.0 - 0.2 / x)) * (x * x * x * x);
  }

  return pressure;
}

//   Uses polint, trapzd

void qromb(double func, double a, double b, double ss){

  int JMAX, JMAXP, K, KM;
  double EPS;
  EPS = 1.0e-8;
  JMAX = 25;
  JMAXP = JMAX + 1;
  K = 8;
  KM = K - 1;
  double dss, h[JMAXP], s[JMAXP];
  h[0] = 1;
 
  for(int j = 0; j <= JMAX - 1; j++){
    
    trapzd(func, a, b, s[j], j);
    
    if(j <= K){
      polint(h[j - KM], s[j - KM], K, 0.0, ss, dss);
    }
    if(abs(dss) <= EPS * abs(ss)) break;
    if (j == JMAX) break;
    
    s[j + 1] = s[j];
    h[j+1] = 0.25 * h[j];
}

void trapzd(){

  double a, b, s, func;

  extern func;

  int n, it, j, tarray[3];
  double del, sum, tnm, x;

  if(n == 1){
    s = 0.5 * (b - a) *(func(a) + func(b));
  }
  else{
    it = pow(2,n-2);
    tnm = it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    sum = 0.0;
  }
  
  s = 0.5 * (s + (b - a) * sum / tnm)

}


void polint(){

  int n, NMAX, ns;
  double dy, x, y, xa[n], ya[n];
  NMAX = 10;
  double den, dif, dift, ho, hp, w, c[NMAX], d[NMAX];
  
  ns = 1;
  dif = abs(x - xa(1));
  
  for(int i = 0; i <= n - 1; i++){

    dift = abs(x - xa[i]);
    if(dift <= dif){
      ns = i;
      dif = dift;
    }
    
    c[i] = ya[i];
    d[i] = ya[i];
  }

  y = ya[ns];
  ns = ns - 1;

  for(int m = 0; m <= n - 2; m++){
    for(int i = 0; i <= n - m - 1; i++){
      ho = xa[i] - x;
      hp = xa[i + m - 1] - x;  
      w = c[i] - d[i - 1];

      if(den == 0.0){
      std::cout << "failure to polint" << std::endl;
      break;
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    if(2 * ns <= n - m){
      dy = c[ns];
      ns = ns;
    }
    y = y + dy;
  } 

}
