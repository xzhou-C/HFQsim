#include <vector>//Vector hearder, defines the vector container class
#include <iostream>//Header that defines the standard input/output stream objects, e.g. cin, cout, cerr
#include <fstream> //Dynamic input
#include <string.h>//functions to manipulate C strings and arrays
#include <stdio.h>//file related input/output
#include <algorithm>//vector copy function, binary search
#include <stdlib.h>// srand, rand
#include <random> // C++11 library for higher-quality randomness
#include <numeric> // std::iota
#include <iterator> // std::distance
#include <omp.h> // openmp
#include <math.h> 
#include "physical_constants.h"
#include "input_parameters.h" //Input value 
#include "fileio.h"
#include "file_array.h"
#include "temperature_profile.h"
#include "strength_profile.h"
#include "creep_law.h"
#include "side_creep.h"
#include "stiffness_matrix.h"
#include "fault.h"
#include "earthquake.h"
#include "romberg_integration.h"
#include "fault_cooling.h"
#include "activ_energy.h"
#include "units.h"

int Romberg::JMAX = 20;
int Romberg::K = 5;
double Romberg::EPS = 1.0e-10;
/** Make operator for multiplying a vector by constant scalar */
std::vector<double> operator*( const std::vector<double>& v, double alpha)
{
  std::vector<double> v2 = v;
  for( unsigned int i = 0; i<v2.size(); i++ ) v2[i] *= alpha; 
  return v2;
}

//-------------------------------------------------------------------------
int main(int argc, char * argv[])
{
  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;
  /* Geometry parameters */
  int nL, nD; // grid size along strike, depth
  double faultLength, faultDepth;
  /* Fritional parameters */
  double fs; // static coefficient of friction
  double tau0; // cohesion
  double StrengthChangeCoef; // strength change coefficient 
  /* Creep parameters */
  double activationEnergy;
  double activationEnergyDamaged; // ratio of activation energy of damaged rock
  double tau_min; // minimum healing time for a full recovery of activation energy
  double tau_max; // maximum healing time 
  int RandomSeed; // random seed for generating activation energy fluctuation
  /* Thermal parameters */
  double SlipZoneWidth; // controls co-seismic temperature change
  double Cooling_distance; // controls the decay of temperature increase
  /* Initial distribution of fields at cell centers */
  std::string ifile_strengthdrops; // stress drop distribution 
  std::string ifile_initSlipDef; // initial slip deficit distribution
  /* Running parameters */
  double maxTime;
  double tSave_min;
  double tSave_max;
  std::string Output_file_prefix; // prefix of the name of all output files

  /* input file, each data separated by \n, order MATTERS! */
  std::ifstream inputFile("input_parameters.txt"); 
  /* check if the file opened correctly */
  if( !inputFile ){
    std::cerr << "Unable to open file input_parameters.txt" << std::endl;
    return 1; // exit with an error code
  }
  /* read file */
  if( inputFile >> nL >> nD >> faultLength >> faultDepth >> fs >> tau0 >> StrengthChangeCoef >>
      activationEnergy >> activationEnergyDamaged >> tau_min >> tau_max >> RandomSeed >> 
      SlipZoneWidth >> Cooling_distance >> ifile_strengthdrops >> ifile_initSlipDef >> 
      maxTime >> tSave_min >> tSave_max >> Output_file_prefix){
      /* Adjust units */
    faultLength = faultLength * Units::km;
    faultDepth = faultDepth * Units::km;
    tau0 = tau0 * Units::MPa;
    activationEnergy = activationEnergy * Units::kJ/Units::mole;
    tau_min = tau_min * Units::day;
    tau_max = tau_max * Units::day;
    SlipZoneWidth = SlipZoneWidth * Units::m;
    Cooling_distance = Cooling_distance * Units::m;
    maxTime = maxTime * Units::year;
    tSave_min = tSave_min * Units::year;
    tSave_max = tSave_max * Units::year;
    std::cerr << "1. nL = " << nL << std::endl;
    std::cerr << "2. nD = " << nD << std::endl;
    std::cerr << "3. Fault Length = " << faultLength/Units::km << " (km)" << std::endl;
    std::cerr << "4. Fault Depth = " << faultDepth/Units::km << " (km)" << std::endl;
    std::cerr << "5. Static Coefficiet of friction = " << fs << std::endl;
    std::cerr << "6. Cohesion = " << tau0/Units::MPa << " (MPa)" << std::endl;
    std::cerr << "7. Strength Change Parameter = " << StrengthChangeCoef << std::endl;
    std::cerr << "8. Activation Energy = " << activationEnergy/(Units::kJ/Units::mole)<< " (kJ/mol)" << std::endl;
    std::cerr << "9. Activation Energy ratio of damaged rock = " << activationEnergyDamaged << std::endl;
    std::cerr << "10. Heal Time Min = " << tau_min/Units::day << " (day)" << std::endl;
    std::cerr << "11. Heal Time Max = " << tau_max/Units::day << " (day)" << std::endl;
    std::cerr << "12. Random Seed value = " << RandomSeed << std::endl;
    std::cerr << "13. Slip Zone Width = " << SlipZoneWidth/Units::m<< " (m)" << std::endl;
    std::cerr << "14. Cooling distance = " << Cooling_distance/Units::m << " (m)" << std::endl;
    std::cerr << "15. Stress drop distribution: "<< ifile_strengthdrops << std::endl;
    std::cerr << "16. Initial slip deficit distribution: "<< ifile_initSlipDef << std::endl;
    std::cerr << "17. Maximum simulation time: " << maxTime/Units::year << " (yr)" << std::endl;
    std::cerr << "18. Initial time to save time evolutions of model parameters: " << tSave_min/Units::year<< " (yr)" << std::endl;
    std::cerr << "19. End time to save time evolutions of model parameters: " << tSave_max/Units::year << " (yr)" << std::endl;
    std::cerr << "20. Output files Prefix: " << Output_file_prefix <<std::endl;
  }else{
    std::cerr << "Error reading input file." << std::endl;
  }
  inputFile.close();

  const double cellLength = faultLength/(double)nL;
  const double cellHeight = faultDepth/(double)nD;

  /*Background Temperatre*/
  BackgroundTemperatureProfile T( In::Tsurface, In::Tgradient);

  /*Static friction coefficient profile*/
  int nCells = nL*nD; 
  std::vector<double> vector_fs(nCells,fs); // profile with a single uniform value
  /* load pre-generated fs distribution */
  // if( !txtfile2faultvec( vector_fs, In::ifile_frictionCoef, nL, nD)){
  //   return -1;
  // }
  
  /* Static strength profile */
  StaticStrengthProfile StatStrength(tau0, vector_fs, In::dsigmadz, nL, nD, cellHeight);
  std::vector<double> taus(nCells,0.0);
  StatStrength.getStaticStrength(taus); 

  /*Dynamic stress drop*/
  std::vector<double> strengthDrops (nCells,0.0); // τs-τa
  /*  load pre-generated stress drop distribution */
  if( !txtfile2faultvec( strengthDrops, ifile_strengthdrops, nL, nD ) ){   
    return -1;
  }

  /* Initial slip deficit */
  std::vector<double> slipDeficit(nCells, 0.0);  // start with zone slip deficit
  /* loading pre-generated slip deficit distribution */
  if( !txtfile2faultvec( slipDeficit, ifile_initSlipDef, nL, nD) ){
    return -1;
  }

  /*Stiffness matrix*/
  StiffnessMatrix K(cellLength, nL, cellHeight, nD, In::rigidity); //cellLength = 70 km/nL, cellHeight = 17.5 km/nD, rigidity of the half space 90 GPa
  FaultActivationEnergy FAE;

  /* set up creep parameters */
  std::vector<double> faultE(nCells, activationEnergy);
  std::vector<double> faultA(nCells, In::arrhAmplitude);
  std::vector<double> faultn(nCells, In::stressExponent );
  std::vector<double> faultdE(nCells); // activation energy after EQ

  /* setup a linear activation energy gradient with depth until zBD*/
  FAE.setAE_linear( faultE, activationEnergy, activationEnergy, In::zBD,nL, nD, cellHeight);

  /* Initialize fault strength profile */
  std::vector<double> strengthProfile(nD);
  std::vector<double> StrengthChangeMatrix(nCells);
  for( int j=0; j<nD; j++){
    double depth = ((double)j+0.5)*cellHeight;
    for( int i=0; i<nL; i++){
      int iCell = i*nD+j;
      StrengthChangeMatrix[iCell] = StrengthChangeCoef;
      CreepLaw C0( faultA[iCell], faultE[iCell], faultn[iCell] );
      strengthProfile[j] = std::min( taus[i*nD+j],
                                     C0.stress( In::loadingStrainRate, T(depth) ) );
    }
  }

  /* Get the  background temperature */
  std::vector<double> faultBkgdT = T.faultgrid_co( nL, nD, cellHeight );
  /* Side creep to smooth out creep motion near fault boundary*/
  sidecreepmask( faultA, faultE, faultn, faultBkgdT,
		                      nL, nD, cellLength, In::xBD);

  /* add random fluctuation to activation energy. do this after side creep */
  std::random_device rd; // non-deterministic random number generator
  if( RandomSeed < 0){
    RandomSeed = rd(); // Mersenne Twister engine
  }
  std::mt19937 gen(RandomSeed);
  std::uniform_int_distribution<> dis(1,5); //Uniform distribution of integer in range [1,5]
  for( int j = 0; j<nD; j++){
    for( int i = 0; i<nL; i++){
      int fluc = dis(gen);
      int iCell = i*nD+j;
      if( fluc == 1 ){
	    faultE[iCell] = (1.0-In::AEasp_ratio) * faultE[iCell];
      }
      else if( fluc == 2 ){
	    faultE[iCell] = (1.0+In::AEasp_ratio) * faultE[iCell];
      }
      else{
        faultE[iCell] = faultE[iCell];
      }
      faultdE[iCell] = activationEnergyDamaged * faultE[iCell];
    }
  }  

  /* fault background temperature, initial stress, arrest stress, initial activation energy */
  std::vector<double> bkgdT(nCells), initStress(nCells), taua(nCells), faultE0(nCells);

  /* update initial stress for plotting only */
  for( int j=0; j<nD; j++ ){
    double depth = ((double)j + 0.5)*cellHeight;
    for( int i=0; i<nL; i++ ){
      int iCell = i*nD+j;
	    taua[iCell] = taus[iCell] - strengthDrops[iCell]; // arrest stress
      faultE0[iCell] = faultE[iCell];
      CreepLaw C( faultA[iCell], faultE[iCell], faultn[iCell] );
	    bkgdT[iCell] = T(depth);
	    initStress[iCell] = std::min( taus[iCell],C.stress( In::loadingStrainRate, bkgdT[iCell] ));
	  }
  }

  /* make healing mechanism stronger with depth, linear dependence on normal stress(depth) */
  std::vector<double> tauHeal(nCells, 0.0);
  for( unsigned int j=0; j<nD; j++){
    //double depth = ((double)j + 0.5)*In::cellHeight;
    for(unsigned int i=0; i<nL; i++){
      tauHeal[i*nD+j] = tau_max - (tau_max-tau_min)*(double)j/(double)nD; // depth has unit of km
    }
  }

  /* Save the initial temperature */
  const std::string ofilename_faultBkgdT(Output_file_prefix + In::ofilesuffix_faultBkgdT );
  faultvec2txtfile( ofilename_faultBkgdT, faultBkgdT*(1.0/Units::K), nL,nD );

  /* Save the Arrhenius amplitude and activation energy */
  const std::string ofilename_Arrheamp( Output_file_prefix + In::ofilesuffix_Arrheamp );
  faultvec2txtfile( ofilename_Arrheamp, faultA*(1.0/( pow(Units::Pa, 3) * Units::second)), nL, nD );

  const std::string ofilename_activE( Output_file_prefix + In::ofilesuffix_activE );
  faultvec2txtfile( ofilename_activE, faultE*(1.0/Units::kJ/Units::mole), nL, nD );

  /* Save the Strength Change Parameter distribution */
  const std::string ofilename_StrengthChange( Output_file_prefix + In::ofilesuffix_StrengthChange );
  faultvec2txtfile( ofilename_StrengthChange, StrengthChangeMatrix, nL, nD);

  /* save tau_Heal distribution */
  const std::string ofilename_tau_heal( Output_file_prefix + In::ofilesuffix_tau_heal );
  faultvec2txtfile( ofilename_tau_heal, tauHeal*(1.0/Units::year), nL, nD );

  /* Save the strength profile */
  const std::string ofilename_initStress( Output_file_prefix + In::ofilesuffix_initStress );
  faultvec2txtfile( ofilename_initStress, initStress*(1.0/Units::MPa), nL, nD );

  /* dynamic strength */
  std::vector<double> taud(nCells,0.0);
  for( int k=0; k<nCells; k++){
    taud[k] = taus[k] - ( taus[k] - taua[k] ) * StrengthChangeMatrix[k];
  }

  /* set initial stress to zero, calculate stress only from slip deficit */
  for( int j=0; j<nD; j++ ){
    double depth = ((double)j + 0.5)*cellHeight;
    for( int i=0; i<nL; i++ ){
      int iCell = i*nD+j;
      initStress[iCell] = 0.0;
    }
  }

  /** Initialize the fault */
  Fault F( nL, nD, cellLength, cellHeight, In::faultWidth,
           taus, taua, taud, StrengthChangeMatrix, faultA, faultn, K, initStress );
  
  F.setSlipDeficit(slipDeficit);

  /* Set containers for storing the fault stress and creep velocity */
  std::vector<double> stress(nCells), creepVel(nCells);
  F.getStress( stress ); 

  /* container for the time of last EQ at each cell */
  std::vector<double> fault_EQtime(nCells, 0.0);

  /* Containers for recording the cumulative slip */
  std::vector<double> totalCreepSlip( nCells, 0.0), totalEqkSlip( nCells, 0.0);

  /* Containers for recording the temperature on the fault. Each loop T_start = bkgdT + fault_dT */
  std::vector<double> T_start(nCells);
  
  /* starting time of each event */
  std::vector<double> tHist;  

  /* slip deficit on the fault at each step  */
  std::vector<double> slipdefHist;

  /* duration of each event */
  std::vector<double> tStep;

  /*container for recording dT on the fault*/
  std::vector<double> fault_dT(nCells, 0.0);

  /*dT on the fault from each time step*/
  std::vector<double> hist_dT;

  /*stress on the fault at the end of each step*/
  std::vector<double> hist_tau;

  /*creep velocity on the fault at each step*/
  std::vector<double> hist_vc;

  /*heat rate at each cell*/
  std::vector<double> dqdt(nCells, 0.0);

  /*record each step as either EQ or creep slip*/
  std::vector<int> slip_type; //1 for creep, 2 for EQ

  /*record slip on the fault at each time step*/
  std::vector<double> faultSlip;

  /* activation energy evolution on the fault at each step */
  std::vector<double> histActE;
  
  /* Container of stress on the fault before loading */
  std::vector<double> stressPre(nCells,0.0);

  /* Container of creep slip */
  std::vector<double> c_slip(nCells, 0.0);
  std::vector<double> eqkSlip( nCells, 0.0 );

  /* Container of slip deficit before and after an EQ */
  std::vector<double> slipDeficitPre;
  std::vector<double> slipDeficitPost;
  
  /* Calculate the stress */
  F.getStress( stress );

  /* Initialize the cooling model */
  FaultCoolingModel FCM( In::density, In::diffusivity, In::specHeat, 
                         SlipZoneWidth, Cooling_distance);

  /* Print header for the earthquake catalog */
  std::cout << "Time_yr, x_km, z_km, Mag_L, Mag_P, Mag_W, Area_km2, StressDrop_MPa\n";

  /* Run the algorithm */
  int i = 0;
  double time = 0.0;
  double timePre; // time before EQ
  /* index for saving files */
  int iSave = 0; // saving iteration
  int iFile = 1; // output file label number
  bool isSaving = false;
  while( i<In::nTimeMax && time < maxTime )
  {
    // check the condition to save parameter history
    if( time >= tSave_min && time <= tSave_max){
      isSaving = true;
    }else
    {
      isSaving = false;
    }
    /* update temperature */
    for(unsigned int k=0; k<nCells; k++){
      T_start[k] = bkgdT[k] + fault_dT[k];
    }
    /* update activation energy */
    FAE.logHeal(faultE,faultdE,faultE0,fault_EQtime,tauHeal,time,nL,nD);
    /* Get the creep velocity on the fault */
    F.getCreepVelocity( creepVel, stress, T_start, faultE);
 
    /* Compute the time to failure */
    double timeStep = F.estimateTimeToFailure( stress, creepVel, In::plateVelocity );
    /* Adjust the time step */
    if( timeStep < 0 ){
      /* Negative implies something went wrong */
        std::cerr << "ERROR: negative time step\n";
        return -1;
      } else if( timeStep > In::maxTimeStep ){
        /* Don't let it get too big or creep rates will be inaccurate */
        timeStep = In::maxTimeStep;
	      //std::cerr << "Routine Report: Using Max time step. "<< std::endl;
      } else if( timeStep < In::minTimeStep ){
        /* Don't let it get too small or we will be waiting forever */
        timeStep = In::minTimeStep;
	      //std::cerr << "Routine Report: Using Min time step. " << std::endl;
      } else{
        //std::cerr << "Routine Report: Time step = " << timeStep << std::endl;
    }
    std::cerr << "Routine Report: Model Time = " << time/Units::year << " (yr)" << std::endl;

    /* Load the fault */
    F.loadFault( In::plateVelocity, creepVel, timeStep );
    /* update slip deficit */
    F.getSlipDeficit( slipDeficit );

    /* update total creep slip */
    for( unsigned int k=0; k<nCells; k++ ){
      c_slip[k] = creepVel[k]*timeStep;
      totalCreepSlip[k] += c_slip[k];
    }

    /* Get the new stress */
    stressPre = stress;
    F.getStress( stress );

    /* set heat generated by creep to zero */
    for( unsigned int k=0; k<nCells; k++) {
	    // dqdt[k] = 0.0*(stressPre[k] + stress[k]) * creepVel[k] * timeStep;
      dqdt[k] = 0.0;
    }
    /* update temperature on the fault */
    FCM.finite_cooling(fault_dT, timeStep, dqdt, nCells);

    /* update model time*/
    time += timeStep;
    
    if( isSaving ){
      slip_type.push_back(1); // step count +1, creep slip type = 1
      histActE.insert(histActE.end(),faultE.begin(),faultE.end());
      slipdefHist.insert(slipdefHist.end(),slipDeficit.begin(),slipDeficit.end());
      hist_vc.insert(hist_vc.end(),creepVel.begin(),creepVel.end());
      faultSlip.insert(faultSlip.end(),c_slip.begin(),c_slip.end());
      hist_tau.insert(hist_tau.end(),stress.begin(),stress.end());
      hist_dT.insert(hist_dT.end(),fault_dT.begin(),fault_dT.end());
      tHist.push_back(time); // record time, time step
      tStep.push_back(timeStep);
      iSave++;
    }

    /* Calculate whether there are any hypocenters */
    unsigned int iHypo, jHypo;
    unsigned int nCrit = F.nCriticalCells(iHypo, jHypo, stress );
    if( nCrit > 1 ){
      std::cerr << "\nWARNING: Multiple hypocenters (" << nCrit << ") ";
    }
    if( nCrit > 0 ){
      /* Output a dot for progress to terminal */
      std::cerr <<"Earthquake Report: Time = "<< time/Units::year<< "  **EQ** "<<std::endl;

      /* Compute the slip deficit before the earthquake */
      F.getSlipDeficit( slipDeficitPre );

      /* Store the prior stress and time */
      stressPre = stress;
      timePre = time;

      /* Compute the earthquake */
      F.computeEarthquake( stress, time, iHypo, jHypo, In::slipVelocity, eqkSlip, fault_EQtime, In::recordTime );

      /* Compute the new slip deficit */
      F.getSlipDeficit( slipDeficitPost );

      /* Compute the earthquake properties */
      Earthquake thisEQ( iHypo, jHypo, timePre, cellLength*cellHeight,
                    slipDeficitPre, slipDeficitPost, stressPre, stress );
        
      /* record EQ slip */
      for( unsigned int k=0; k<nCells; k++ ){
	      totalEqkSlip[k] += eqkSlip[k];
        dqdt[k] = stress[k]*eqkSlip[k];
        
      }
      /* update temperature change on the fault */
      FCM.finite_cooling(fault_dT, time-timePre, dqdt, nCells);

      /* update temperature */
      for(unsigned int k=0; k<nCells; k++){
        T_start[k] = bkgdT[k] + fault_dT[k];
      }
      /* check max on-fault temperature increase, useful for debugging */
      /*std::cerr <<"Routine Report: Current max dT on the fault, index= "
        <<std::max_element(fault_dT.begin(),fault_dT.end())-fault_dT.begin()
        <<" , dT = "<<*std::max_element(fault_dT.begin(),fault_dT.end())<<std::endl;*/

      /* update activation energy */
      FAE.logHeal(faultE,faultdE,faultE0,fault_EQtime,tauHeal,time,nL,nD);
      /* Get the creep velocity on the fault */
      F.getCreepVelocity( creepVel, stress, T_start, faultE);

      if( isSaving ){
        /* record the duration*/
        tStep.push_back(time-timePre);
        /* record the post seismic time */
        tHist.push_back(time);
        /* EQ slip */
        faultSlip.insert(faultSlip.end(),eqkSlip.begin(),eqkSlip.end());
        /* record activation energy */
        histActE.insert(histActE.end(),faultE.begin(),faultE.end());
        slipdefHist.insert(slipdefHist.end(),slipDeficitPost.begin(),slipDeficitPost.end());
        /* record type of slip, EQ slip type = 2 */
        slip_type.push_back(2);
        /* record temperature and stress*/
        hist_dT.insert(hist_dT.end(),fault_dT.begin(),fault_dT.end());
        hist_tau.insert(hist_tau.end(),stress.begin(),stress.end());
        /* record post-seismic creep velocity */
        hist_vc.insert(hist_vc.end(),creepVel.begin(),creepVel.end());
        iSave++;
      } /* end if isSaving */
        
      if( i > 0 ){
        /* Output the earthquake catalog: time, x, z, magL, magP, magW, area, stressdrop*/
        printf( "%9.6E, %6.2f, %5.2f, %4.2f, %4.2f, %4.2f, %6.2f, %5.2f\n",
                  timePre/Units::year,
                  (iHypo+0.5)*cellLength/Units::km ,
                  (jHypo+0.5)*cellHeight/Units::km ,
                  thisEQ.localMagnitude() ,
                  thisEQ.potMagnitude() ,
                  thisEQ.momMagnitude( In::rigidity ) ,
                  thisEQ.ruptureArea() / (Units::km*Units::km) ,
                  thisEQ.staticStressDrop()/Units::MPa );
        } /* end if output */
        
    } /*END earthquake */
      
    /* save output files and clear memory */
    if(iSave > 0 && ( (iSave+1) % 30000 == 0 || time >= tSave_max) ){     
      std::cerr << " ############### Saving Outputs ############### " << std::endl;
      /* model time */
      const std::string ofilename_timeStep( Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_timeStep );
      x_txtfile( ofilename_timeStep, tHist*(1.0/Units::year));

	    /* slip deficit */        
	    const std::string ofilename_slipdefHist( Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_slipdefHist );
      x_txtfile( ofilename_slipdefHist, slipdefHist*(1.0/Units::m));

	    /* slip */
      const std::string ofilename_faultSlip( Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_faultSlip );
      x_txtfile( ofilename_faultSlip, faultSlip*(1.0/Units::m));

	    /* stress */
      const std::string ofilename_hist_tau( Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_hist_tau );
      x_txtfile( ofilename_hist_tau, hist_tau*(1.0/Units::MPa));

	    /* temperature change */
      const std::string ofilename_hist_dT( Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_hist_dT );
      x_txtfile( ofilename_hist_dT, hist_dT*(1.0/Units::K));

	    /* activation energy */
      const std::string ofilename_hist_actE( Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_hist_actE );
      x_txtfile( ofilename_hist_actE, histActE*(1.0/Units::kJ/Units::mole));

	    /* creep velocity */
      const std::string ofilename_hist_vc( Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_hist_vc );
      x_txtfile( ofilename_hist_vc, hist_vc*(1.0*Units::second/Units::m));

      /* slip type */
	    const std::string ofilename_slip_type(Output_file_prefix + std::to_string(iFile) + In::ofilesuffix_slip_type);
	    x_txtfile( ofilename_slip_type, slip_type);

	    tHist.clear();
      slipdefHist.clear();
      faultSlip.clear();
      hist_tau.clear();
      hist_dT.clear();
      histActE.clear();
      hist_vc.clear();
      tHist.clear();
      slip_type.clear();
      tStep.clear();
      /* Reset iSave */
	    iSave = 0;
      iFile++;
    }/* END saving */

    i++;

  } /*END fault algorithm */

  /* Break the line */
  std::cerr << std::endl;

  /* Report why we finished the loop */
  if( i == In::nTimeMax ) std::cerr << "Finishing Report: Max number of iterations reached\n";
  else std::cerr << "Finishing Report: Max time reached\n";

  /* Output files */
  std::cerr << "Saving output files. " << std::endl; 
  /* total creep slip */
  const std::string ofilename_creepSlip( Output_file_prefix + In::ofilesuffix_creepSlip );
  faultvec2txtfile( ofilename_creepSlip, totalCreepSlip*(1.0/Units::m), nL, nD );

  /* total EQ slip */
  const std::string ofilename_eqSlip( Output_file_prefix + In::ofilesuffix_eqSlip );
  faultvec2txtfile( ofilename_eqSlip, totalEqkSlip*(1.0/Units::m), nL, nD );

  /* final slip deficit */
  F.getSlipDeficit( slipDeficit );
  const std::string ofilename_slipdef( Output_file_prefix + In::ofilesuffix_slipdef );
  faultvec2txtfile( ofilename_slipdef, slipDeficit*(1.0/Units::m), nL, nD );

  /* final stress */
  const std::string ofilename_finalStress( Output_file_prefix + In::ofilesuffix_finalStress );
  faultvec2txtfile( ofilename_finalStress, stress*(1.0/Units::MPa), nL, nD );
}

/* END */         
                 











  
