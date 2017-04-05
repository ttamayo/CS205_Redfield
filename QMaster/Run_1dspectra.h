/*
 * Run_1dspectra.h
 *
 *  Created on: Dec 16, 2013
 *      Author: christoph
 */

#ifndef RUN_1DSPECTRA_H_
#define RUN_1DSPECTRA_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "../headers.h"

class Run_1dSpectra
{
public:
  System0ex1ex* system;
  OpenCL_Init* OpenCLinfo;
  Param_HEOM_Simulation* param_heom_sim;
  Param_Absorption_Spectra* param_1d_spec;
  Dipole_Operator_0ex1ex* dipole_operator;
  const char *fn_output;
  double diff;
private:
  SigmaHelpMatrices *sigma_Host;
  int Ntuples;
  struct timeval t0, t1;
  int NTP;
  bool use_specified_laserDirection;
public:
  Run_1dSpectra(System0ex1ex &INsystem,
                          OpenCL_Init& INOpenCLinfo,
                          Param_HEOM_Simulation &INparam_heom_sim,
                          Param_Absorption_Spectra &INparam_1d_spec,
                          Dipole_Operator_0ex1ex &INdipole_operator,
                          const char *INfn_output
                         );
  void execute_1dspectrum();
};
Run_1dSpectra::Run_1dSpectra(System0ex1ex &INsystem,
        OpenCL_Init& INOpenCLinfo,
        Param_HEOM_Simulation &INparam_heom_sim,
        Param_Absorption_Spectra &INparam_1d_spec,
        Dipole_Operator_0ex1ex &INdipole_operator,
        const char *INfn_output) :
	system(&INsystem),
	OpenCLinfo(&INOpenCLinfo),
	param_heom_sim (&INparam_heom_sim),
	param_1d_spec(&INparam_1d_spec),
	dipole_operator(&INdipole_operator),
	fn_output(INfn_output)
{

  fprintf(stdout, "#\n# initialize parameter for absorption spectrum\n");
  use_specified_laserDirection=false;
  diff = 0;
  Ntuples = 0;
  sigma_Host = NULL;
  NTP = param_heom_sim->Nsteps * param_1d_spec->padding; // steps along time including padding
  //max energy range
  double Eminall = -NTP * 0.5 * const_hbar / const_invcmtomeV * 2.*M_PI / ((NTP - 1) * param_heom_sim->dt);
  double Emaxall = +NTP * 0.5 * const_hbar / const_invcmtomeV * 2.*M_PI / ((NTP - 1) * param_heom_sim->dt);
  if(param_1d_spec->specify_Enrange == true)
  {
    // check if specified Erange fits into [Eminall, Emaxall]
    if(param_1d_spec->Emin < Eminall)
    {
      fprintf(stderr, "error: specified En_range Emin=%e must be larger than %e (constraint by FFT grid)\n", param_1d_spec->Emin, Eminall);
      fprintf(stderr, "change either Emin or reduce the time step dt\n");
      exit(1);
    }
    if(param_1d_spec->Emax > Emaxall)
    {
      fprintf(stderr, "error: specified En_range Emax=%e must be smaller than %e (constraint by FFT grid)\n", param_1d_spec->Emax, Emaxall);
      fprintf(stderr, "change either Emax or reduce the time step dt\n");
      exit(1);
    }
  }
  else
  {
	  param_1d_spec->Emax=system->Eigenvals[system->Nsites1ex-1]/const_invcmtomeV+system->Ediffmax_invcm;
	  param_1d_spec->Emin=system->Eigenvals[0]/const_invcmtomeV-system->Ediffmax_invcm;
  }
}

void Run_1dSpectra::execute_1dspectrum()
{
  fftw_plan fftplan1d;
  int sizeN = NTP * sizeof(fftw_complex);
  char fn[500];
  int num_shots;
  FILE *fo;

  fftw_complex *data_t, *data_f, *data_fsum;

  // arrays to hold time line and its FFT to energy space
  data_t = (fftw_complex*)fftw_malloc(sizeN);
  data_f = (fftw_complex*)fftw_malloc(sizeN);
  data_fsum = (fftw_complex*)fftw_malloc(sizeN);
  for(int i = 0; i < NTP; i++)
  {
    data_t[i][0] = 0.0;
    data_t[i][1] = 0.0;
    data_fsum[i][0] = 0.0;
    data_fsum[i][1] = 0.0;
  }

  fftplan1d = fftw_plan_dft_1d(NTP, (fftw_complex*)data_t, (fftw_complex*)data_f, FFTW_BACKWARD, FFTW_ESTIMATE);

  if(param_1d_spec->rotational_average == true)
  {
    num_shots = 3;
  }
  else
  {
    num_shots = 1;
    dipole_operator->set_new_laser_direction(param_1d_spec->Pol[0], param_1d_spec->Pol[1], param_1d_spec->Pol[2]); // default Pol=(1,0,0) if not default Pol is specified by keyword laser_polarization
  }
  fprintf(stdout, "# specified dipole moments: %s\n", param_1d_spec->dipole_moments);

//   int Nmod = param_heom_sim->Nsteps / 10; // every 10% progress report to stdout
//   fprintf(stdout, "# initialize propagation ...\n");
//   sigma_Host = new SigmaHelpMatrices(system->Nsites, param_heom_sim->Nmax, system->Ncoupled);
//   Ntuples = sigma_Host->AllsigmaTuples.size();
//   Propagation evolution(*system, *OpenCLinfo, *sigma_Host, param_heom_sim->Nsteps, param_heom_sim->dt, param_heom_sim->Nmax);
//   fprintf(stdout, "# ... done\n");
//   LiouvilleHExciton liouHExc(*system, *OpenCLinfo, Ntuples);
//   LiouvillePhononLowTemp liouPhon(*system, *OpenCLinfo,  *sigma_Host, param_heom_sim->Nmax);
//   //
//   // AddSites to LiouvillPhononLowTemp
//   for(int i = 0; i < system->all_spectral_dens.size(); i++)
//     {
//       system->all_spectral_dens[i]->plot(fn_output, system->Ediffmax_invcm);
//       liouPhon.AddSite(system->all_spectral_dens[i]->CoupledSite, // the shift +1 is taken care of in system0exRC1ex.h
//                        system->all_spectral_dens[i]->Npeaks,
//                        system->all_spectral_dens[i]->lambda,
//                        system->all_spectral_dens[i]->gamma);
//     }
//   fprintf(stdout, "# initialize coupling to the phonon-bath ...\n");
//   liouPhon.initialize();
//   fprintf(stdout, "# ... done\n");
//   evolution.Add_Liouville(liouHExc);
//   evolution.Add_Liouville(liouPhon);
  for(int shot = 0; shot < num_shots; shot++)
  {
    //if(shot==0) continue;
  int Nmod = param_heom_sim->Nsteps / 10; // every 10% progress report to stdout
  fprintf(stdout, "# initialize propagation ...\n");
  sigma_Host = new SigmaHelpMatrices(system->Nsites, param_heom_sim->Nmax, system->Ncoupled);
  Ntuples = sigma_Host->AllsigmaTuples.size();
  Propagation evolution(*system, *OpenCLinfo, *sigma_Host, param_heom_sim->Nsteps, param_heom_sim->dt, param_heom_sim->Nmax);
  fprintf(stdout, "# ... done\n");
  LiouvilleHExciton liouHExc(*system, *OpenCLinfo, Ntuples);
  LiouvillePhononLowTemp liouPhon(*system, *OpenCLinfo,  *sigma_Host, param_heom_sim->Nmax);
  //
  // AddSites to LiouvillPhononLowTemp
  for(int i = 0; i < system->all_spectral_dens.size(); i++)
    {
      system->all_spectral_dens[i]->plot(fn_output, system->Ediffmax_invcm);
      liouPhon.AddSite(system->all_spectral_dens[i]->CoupledSite, // the shift +1 is taken care of in system0exRC1ex.h
                       system->all_spectral_dens[i]->Npeaks,
                       system->all_spectral_dens[i]->lambda,
                       system->all_spectral_dens[i]->gamma);
    }
  fprintf(stdout, "# initialize coupling to the phonon-bath ...\n");
  liouPhon.initialize();
  fprintf(stdout, "# ... done\n");
  evolution.Add_Liouville(liouHExc);
  evolution.Add_Liouville(liouPhon);
    
    
    
    
    if (param_1d_spec->rotational_average == false and (strcmp(param_1d_spec->dipole_moments, "effective dipole moments eigen") != 0) and (strcmp(param_1d_spec->dipole_moments, "effective dipole moments site") != 0))
    {
      fprintf(stdout, "# laser polarization along l=(%e,%e,%e)\n", param_1d_spec->Pol[0], param_1d_spec->Pol[1], param_1d_spec->Pol[2]);
      use_specified_laserDirection = true;
    }
    if ((param_1d_spec->rotational_average  == true) && (shot == 0))
    {
    	dipole_operator->set_new_laser_direction(1.0, 0.0, 0.0);
    }
    if ((param_1d_spec->rotational_average  == true) && (shot == 1))
    {
    	dipole_operator->set_new_laser_direction(0.0, 1.0, 0.0);
    }
    if ((param_1d_spec->rotational_average  == true) && (shot == 2))
    {
    	dipole_operator->set_new_laser_direction(0.0, 0.0, 1.0);
    }
    if (param_1d_spec->rotational_average  == true)
    {
      fprintf(stdout, "#\n# perform rotational average shot=%i\n", shot + 1);
    }
    dipole_operator->create_mum_mup();
    Matrix rho0(system->Nsites);
    Matrix rhoInit(system->Nsites);
    Matrix Help3(system->Nsites);
    rho0.assignEntry(0, 0, 1.0);
    MatrixMuliplication(&rhoInit, dipole_operator->mup, &rho0, system->Nsites);
    evolution.initialize(rhoInit);
    // returns the trace of mu_minus * rho, set DNT=1
    ReturnTrace_mumRho SpectrumOut(*system, * OpenCLinfo, 1, dipole_operator->mum, &Help3, 0);
    evolution.Add_ManipulatingObservable(SpectrumOut);
    fprintf(stdout, "#\n# starting propagation ...\n");
    gettimeofday(&t0, 0);
    evolution.propagate(Nmod);
    cout << "# ... done" << endl;
    gettimeofday(&t1, 0);
    double diff_local = (1000000.0 * (t1.tv_sec - t0.tv_sec) + t1.tv_usec - t0.tv_usec) / 1000000.0;
    diff += diff_local;
    for(int i = 0; i < param_heom_sim->Nsteps; i++)
    {
      data_t[i][0] += SpectrumOut.tracer[i];
      data_t[i][1] += SpectrumOut.tracei[i];
    }
    evolution.clear_Sigma();
    evolution.clear_Observable();
  cout << "# free Device-Memory ..." << endl;
  liouPhon.freeDeviceMemory();
  evolution.freeDeviceMemory();
  liouHExc.freeDeviceMemory();
    cout << "# ... done" << endl;
   }
  cout << "# free OpenCl ..." << endl;
  clFlush(OpenCLinfo->queue);
  clFinish(OpenCLinfo->queue);
  cout << "# ... done" << endl;
  /*
    *********************************************************
     * Save output values to file
     *********************************************************
    */
  sprintf(fn, "%s_spectrum1d.dat", fn_output);
  fo = fopen(fn, "w");
  if (fo == NULL)
  {
    fprintf(stderr, "Cannot open output parameter file \"%s\". ABORT.", fn);
    exit(1);
  }
  time_t t;
  time(&t);
  fprintf(fo, "# output generated by QMaster-%s on %s", INPUT_FORMAT_VERSION, ctime(&t));
  if(strncmp(param_1d_spec->dipole_moments, "effective dipole moments eigen", strlen(param_1d_spec->dipole_moments)) == 0)
  {
    fprintf(fo, "# %s\n", param_1d_spec->dipole_moments);
  }
  if(strncmp(param_1d_spec->dipole_moments, "effective dipole moments site", strlen(param_1d_spec->dipole_moments)) == 0)
  {
    fprintf(fo, "# %s\n", param_1d_spec->dipole_moments);
  }
  if(param_1d_spec->rotational_average == true)
  {
    fprintf(fo, "# xyz-rotational average\n");
  }
  if(use_specified_laserDirection == true)
  {
    fprintf(fo, "# laser polarization along l=(%e,%e,%e)\n", param_1d_spec->Pol[0], param_1d_spec->Pol[1], param_1d_spec->Pol[2]);
  }
  fprintf(fo, "# output: absorption spectrum, time signal\n#\n");
  fprintf(fo, "#  t in s     Re[signal]   Im[signal] \n");
  for(int i = 0; i < NTP; i++)
  {
    fprintf(fo, "%e %e %e\n", i * param_heom_sim->dt, data_t[i][0], data_t[i][1]);
  }
  fclose(fo);
  // perform FFT
  data_t[0][0] = 0.5 * (data_t[0][0] + data_t[NTP - 1][0]);
  data_t[0][1] = 0.5 * (data_t[0][1] + data_t[NTP - 1][1]);
  fftw_execute(fftplan1d);
//   for(int i=0;i<NTP;i++) { data_fsum[i][0]+=data_f[i][0]; data_fsum[i][1]+=data_f[i][1]; }
  sprintf(fn, "%s_spectrum1d_freq.dat", fn_output);
  fo = fopen(fn, "w");
  if (fo == NULL)
  {
    fprintf(stderr, "Cannot open output parameter file \"%s\". ABORT.", fn);
    exit(1);
  }
  time(&t);
  fprintf(fo, "# output generated by QMaster-%s on %s", INPUT_FORMAT_VERSION, ctime(&t));
  if(strncmp(param_1d_spec->dipole_moments, "effective dipole moments eigen", strlen(param_1d_spec->dipole_moments)) == 0)
  {
    fprintf(fo, "# %s\n", param_1d_spec->dipole_moments);
  }
  if(strncmp(param_1d_spec->dipole_moments, "effective dipole moments site", strlen(param_1d_spec->dipole_moments)) == 0)
  {
    fprintf(fo, "# %s\n", param_1d_spec->dipole_moments);
  }
  if(param_1d_spec->rotational_average == true)
  {
    fprintf(fo, "# xyz-rotational average\n");
  }
  if(use_specified_laserDirection == true)
  {
    fprintf(fo, "# laser polarization along l=(%e,%e,%e)\n", param_1d_spec->Pol[0], param_1d_spec->Pol[1], param_1d_spec->Pol[2]);
  }
  fprintf(fo, "# output: absorption spectrum\n#\n");
  fprintf(fo, "#  w in cm^-1   Re[signal]  Im[signal] \n");
  for(int w = 0; w < NTP; w++)
  {
    double w_invcm;
    int i = w;
    if (w < NTP / 2)  i += NTP / 2; // w=0..NTP/2-1 gives the negative freq
    if (w >= NTP / 2) i -= NTP / 2; // w=NTP/2..NTP gives the positive freq
    w_invcm = double(w - NTP / 2.0) * const_hbar / const_invcmtomeV * 2.0 * M_PI / ((NTP - 1) * param_heom_sim->dt);
    if(w_invcm > param_1d_spec->Emin and w_invcm < param_1d_spec->Emax)
    {
      fprintf(fo, "%e %e %e\n", w_invcm, data_f[i][0], data_f[i][1]);
    }
  }
  fclose(fo);
}

#endif /* GPU_RUN1DSPECTRA_H_ */
