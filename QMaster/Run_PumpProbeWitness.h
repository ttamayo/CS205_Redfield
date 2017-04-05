/*
 * Run_PumpProbeWitness.h
 *
 *  Created on: Feb 11, 2014
 *      Author: christoph kreisbeck
 */

#ifndef Run_PUMPPROBEWITNESS_H_
#define RUN_PUMPPROBEWITNESS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "headers.h"
#include <iostream>

class calculate_PumpProbeWitness
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Nmax;
   double dt;
   int Nsteps;
   int DNT;
   double** sd_lambda_invcm;
   double** sd_gamma_fs;
   double** sd_gammapos_invcm;
   vector<int> sd_CoupledSite;
   vector<int> sd_Npeaks;
   const char *fn_output;
   int GBRP;
   int GBNR;
   int SERP;
   int SENR;
   int ESARP;
   int ESANR;
   int Ndata;
   double_complex *data_SERP_t;
   double_complex *data_SENR_t;
   double_complex *data_GBRP_t;
   double_complex *data_GBNR_t;
   double_complex *data_ESARP_t;
   double_complex *data_ESANR_t;
   double_complex *data_sum_t;
   char Fileheader_calculatedPath[255];
private:
   SigmaHelpMatrices *sigma_Host;
   int Ntuples;
   double* lambda;
   double_complex* gamma;
   double Ediffmax_invcm;
public:
   calculate_PumpProbeWitness(System& INsystem,
                              OpenCL_Init& INOpenCLinfo,
                              int INNmax,
                              double INdt,
                              int INNsteps,
                              int INDNT,
                              double** INsd_lambda_invcm,
                              double** INsd_gamma_fs,
                              double** INsd_gammapos_invcm,
                              vector<int> INsd_CoupledSite,
                              vector<int> INsd_Npeaks,
                              const char* INfn_output,
                              int INGBRP,
                              int INGBNR,
                              int INSERP,
                              int INSENR,
                              int INESARP,
                              int INESANR
                             );
   void execute_calculate_PumpProbeWitness();
   void execute_calculate_Redfield_PumpProbeWitness(const char* method);
   void execute_calculate_PumpProbeWitness(const char* dipole_moments, const char* method);
   void execute_calculate_PumpProbeWitness_rotav_4shots(const char* method);
   void execute_calculate_PumpProbeWitness_rotav_10shots(const char* method);
   void write_data(const char* info_shotdipole);
};


calculate_PumpProbeWitness::calculate_PumpProbeWitness(System& INsystem,
      OpenCL_Init& INOpenCLinfo,
      int INNmax,
      double INdt,
      int INNsteps,
      int INDNT,
      double** INsd_lambda_invcm,
      double** INsd_gamma_fs,
      double** INsd_gammapos_invcm,
      vector<int> INsd_CoupledSite,
      vector<int> INsd_Npeaks,
      const char* INfn_output,
      int INGBRP,
      int INGBNR,
      int INSERP,
      int INSENR,
      int INESARP,
      int INESANR
                                                      )  :
   system(&INsystem),
   OpenCLinfo(&INOpenCLinfo),
   Nmax(INNmax),
   dt(INdt),
   Nsteps(INNsteps),
   DNT(INDNT),
   sd_lambda_invcm(INsd_lambda_invcm),
   sd_gamma_fs(INsd_gamma_fs),
   sd_gammapos_invcm(INsd_gammapos_invcm),
   sd_CoupledSite(INsd_CoupledSite),
   sd_Npeaks(INsd_Npeaks),
   fn_output(INfn_output),
   GBRP(INGBRP),
   GBNR(INGBNR),
   SERP(INSERP),
   SENR(INSENR),
   ESARP(INESARP),
   ESANR(INESANR)
{
   fprintf(stdout, "#\n# Initialize parameter for pump-probe witness\n");
   sprintf(Fileheader_calculatedPath,"# output: sum of all specified pathways:");
   // initialize data for time lines
   Ndata=Nsteps/DNT+1;
   if(SERP==1)
   {
      sprintf(Fileheader_calculatedPath,"%s SE-RP,",Fileheader_calculatedPath);
      data_SERP_t=new double_complex[Ndata];
      for(int i=0; i<Ndata; i++)
      {
         data_SERP_t[i]=double_complex(0.,0.);
      }
   }
   if(SENR==1)
   {
      sprintf(Fileheader_calculatedPath,"%s SE-NR,",Fileheader_calculatedPath);
      data_SENR_t=new double_complex[Ndata];
      for(int i=0; i<Ndata; i++)
      {
         data_SENR_t[i]=double_complex(0.,0.);
      }
   }
   if(GBRP==1)
   {
      sprintf(Fileheader_calculatedPath,"%s GB-RP,",Fileheader_calculatedPath);
      data_GBRP_t=new double_complex[Ndata];
      for(int i=0; i<Ndata; i++)
      {
         data_GBRP_t[i]=double_complex(0.,0.);
      }
   }
   if(GBNR==1)
   {
      sprintf(Fileheader_calculatedPath,"%s GB-NR,",Fileheader_calculatedPath);
      data_GBNR_t=new double_complex[Ndata];
      for(int i=0; i<Ndata; i++)
      {
         data_GBNR_t[i]=double_complex(0.,0.);
      }
   }
   if(ESARP==1)
   {
      sprintf(Fileheader_calculatedPath,"%s ESA-RP,",Fileheader_calculatedPath);
      data_ESARP_t=new double_complex[Ndata];
      for(int i=0; i<Ndata; i++)
      {
         data_ESARP_t[i]=double_complex(0.,0.);
      }
   }
   if(ESANR==1)
   {
      sprintf(Fileheader_calculatedPath,"%s ESA-NR,",Fileheader_calculatedPath);
      data_ESANR_t=new double_complex[Ndata];
      for(int i=0; i<Ndata; i++)
      {
         data_ESANR_t[i]=double_complex(0.,0.);
      }
   }
   data_sum_t=new double_complex[Ndata];;
   for(int i=0; i<Ndata; i++)
   {
      data_sum_t[i]=double_complex(0.,0.);;
   }

   // specify maximal difference between to exction eigen energies to specify plot range of the spectral density (plot energy range where differences in exciton
   // energies live (extend by factor of 1.5 to be sure that important part is plotted))
   system->computeEigenvects_Eigenvals();
   Ediffmax_invcm=2.*(system->Eigenvals[system->Nsites1-1]-system->Eigenvals[0])/const_invcmtomeV;
   // for monomer difference of eigenenergies=0
   // need additional automatic tool specify the range of the spectral density
   // posistion of maximal shifted peak
   double peakposmax=0.;
   for(int i=0; i<sd_CoupledSite.size(); i++)
   {
      for(int k=0; k<sd_Npeaks[i]; k++)
      {
         if(sd_gammapos_invcm[i][k]>peakposmax)
         {
            peakposmax=sd_gammapos_invcm[i][k];
         }
      }
   }
   Ediffmax_invcm+=peakposmax;
   //sometimes, monomer and unshifted DL spectral density, last criteria use reorganization energy
   double lambdamax=0.;
   for(int i=0; i<system->Nsites; i++)
   {
      if(system->lambda[i]>lambdamax)
      {
         lambdamax=system->lambda[i];
      }
   }
   Ediffmax_invcm+=10.*lambdamax/const_invcmtomeV;
}


void calculate_PumpProbeWitness::execute_calculate_PumpProbeWitness(const char* dipole_moments, const char* method)
{
   fprintf(stdout,"# specified dipole moments: %s\n",dipole_moments);
   if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
   {
      if(strncmp(method,"HEOM",strlen(method))==0)
      {
         execute_calculate_PumpProbeWitness();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_PumpProbeWitness(method);
      }
   }
   else if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
   {
      Vector *mu;
      mu=new Vector(system->Nsites1);
      Vector *help_vec;
      help_vec=new Vector(system->Nsites1);
      for(int i=0; i < system->Nsites1; i++)
      {
         // if effective dipole moments in eigenbasis. they are stored in x-component of the DipoleMoments
         mu->entry[i]=system->DipoleMoments[3*i+0];
      }
      MatrixVectorMuliplication(help_vec,system->UT, mu, system->Nsites1);
      for(int i=0; i < system->Nsites1; i++)
      {
         system->DipoleOperator[i+1]=help_vec->entry[i];
         system->DipoleMoments[3*i+0]=help_vec->entry[i];
      }
      if(strncmp(method,"HEOM",strlen(method))==0)
      {
         execute_calculate_PumpProbeWitness();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_PumpProbeWitness(method);
      }
   }
   else
   {
      fprintf(stderr,"effective dipole moments not specified correctly. ABORT.\n");
      exit(1);
   }
   write_data(dipole_moments);
}



void calculate_PumpProbeWitness::write_data(const char* info_shotdipole)
{
   char fn[500];
   FILE *fo;
   time_t t;
   // FIXME: in the end w1 and w3 seem to be swapped ??? ugly hack to write out transpose: FIXME: same for (t1,t3) timeline!
   if(SERP==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_delta-pulse_witness_SERP.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: delta-pulse pump-probe witness pathway SE-RP\n#\n");
      fprintf(fo,"#  T_delay in s    Re[signal]    Im[signal]\n");
      for (int i=0; i<Ndata; i++)
      {
         fprintf(fo, "%e %e %e\n",i*dt*DNT, real(data_SERP_t[i]), imag(data_SERP_t[i]));
      }
      fclose(fo);
      delete [] data_SERP_t;
   }
   //
   if(GBRP==1)
   {
      sprintf(fn,"%s_delta-pulse_witness_GBRP.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: delta-pulse pump-probe witness pathway GB-RP\n#\n");
      fprintf(fo,"#  T_delay in s    Re[signal]    Im[signal]\n");
      for (int i=0; i<Ndata; i++)
      {
         fprintf(fo, "%e %e %e\n",i*dt*DNT, real(data_GBRP_t[i]), imag(data_GBRP_t[i]));
      }
      fclose(fo);
      delete [] data_GBRP_t;
   }
   //
   if(ESARP==1)
   {
      sprintf(fn,"%s_delta-pulse_witness_ESARP.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: delta-pulse pump-probe witness pathway ESA-RP\n#\n");
      fprintf(fo,"#  T_delay in s    Re[signal]    Im[signal]\n");
      for (int i=0; i<Ndata; i++)
      {
         fprintf(fo, "%e %e %e\n",i*dt*DNT, real(data_ESARP_t[i]), imag(data_ESARP_t[i]));
      }
      fclose(fo);
      delete [] data_ESARP_t;
   }

   if(SENR==1)
   {
      sprintf(fn,"%s_delta-pulse_witness_SENR.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: delta-pulse pump-probe witness pathway SE-NR\n#\n");
      fprintf(fo,"#  T_delay in s    Re[signal]    Im[signal]\n");
      for (int i=0; i<Ndata; i++)
      {
         fprintf(fo, "%e %e %e\n",i*dt*DNT, real(data_SENR_t[i]), imag(data_SENR_t[i]));
      }
      fclose(fo);
      delete [] data_SENR_t;
   }
   //
   if(GBNR==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_delta-pulse_witness_GBNR.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: delta-pulse pump-probe witness pathway GB-NR\n#\n");
      fprintf(fo,"#  T_delay in s    Re[signal]    Im[signal]\n");
      for (int i=0; i<Ndata; i++)
      {
         fprintf(fo, "%e %e %e\n",i*dt*DNT, real(data_GBNR_t[i]), imag(data_GBNR_t[i]));
      }
      fclose(fo);
      delete [] data_GBNR_t;
   }
   //
   if(ESANR==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_delta-pulse_witness_ESANR.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: delta-pulse pump-probe witness pathway ESA-NR\n#\n");
      fprintf(fo,"#  T_delay in s    Re[signal]    Im[signal]\n");
      for (int i=0; i<Ndata; i++)
      {
         fprintf(fo, "%e %e %e\n",i*dt*DNT, real(data_ESANR_t[i]), imag(data_ESANR_t[i]));
      }
      fclose(fo);
      delete [] data_ESANR_t;
   }
// TOTAL Witness data_sum_t
   /// TIME line output
   sprintf(fn,"%s_delta-pulse_witness_sum.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time(&t);
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# %s\n",info_shotdipole);
   fprintf(fo,"%s time signal\n#\n",Fileheader_calculatedPath);
   fprintf(fo,"#  T_delay in s    Re[signal]    Im[signal]\n");
   for (int i=0; i<Ndata; i++)
   {
      fprintf(fo, "%e %e %e\n",i*dt*DNT, real(data_sum_t[i]), imag(data_sum_t[i]));
   }
   fclose(fo);
   delete [] data_sum_t;
}


void calculate_PumpProbeWitness::execute_calculate_PumpProbeWitness()
{
   system->create_mum_mup();
   fprintf(stdout, "# # initialize propagation method=HEOM ...\n");
   sigma_Host=new SigmaHelpMatrices(system->Nsites, Nmax, system->Ncoupled);
   Ntuples=sigma_Host->AllsigmaTuples.size();
   fprintf(stdout, "# ... done\n");
   LiouvilleHExciton liouHexci(*system, *OpenCLinfo, Ntuples);
   LiouvillePhononLowTemp liouPhon(*system, *OpenCLinfo, *sigma_Host, Nmax);
   for(int i=0; i<sd_CoupledSite.size(); i++)
   {
      //handle original DL spectral density
      if(sd_Npeaks[i]==1 and sd_gammapos_invcm[i][0]==0.0)
      {
         lambda=new double[sd_Npeaks[i]];
         gamma=new double_complex[sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15), 0.);
         }
         plot_spectral_density(sd_CoupledSite[i], sd_Npeaks[i], lambda, gamma, fn_output, Ediffmax_invcm);
         // NOTE since 1d spectra get the zero exciton state |0_ex><0_ex|, I need the "+1" in sd_CoupledSite[i]+1
         liouPhon.AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
      }
      else
      {
         lambda=new double[2*sd_Npeaks[i]];
         gamma=new double_complex[2*sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[2*k]=sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            lambda[2*k+1]=sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            gamma[2*k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),+sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
            gamma[2*k+1]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),-sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
         }
         plot_spectral_density(sd_CoupledSite[i], 2*sd_Npeaks[i], lambda, gamma, fn_output, Ediffmax_invcm);
         // NOTE since 1d spectra get the zero exciton state |0_ex><0_ex|, I need the "+1" in sd_CoupledSite[i]+1
         liouPhon.AddSite(sd_CoupledSite[i]+1,2*sd_Npeaks[i], lambda, gamma);
      }
   }
   fprintf(stdout, "# initialize coupling to the phonon-bath ...\n");
   liouPhon.initialize();
   fprintf(stdout, "# ... done\n");

   Matrix *rho0;
   rho0=new Matrix(system->Nsites);
   Matrix *rhoInit;
   rhoInit=new Matrix(system->Nsites); //in Site Basis
   Matrix *Help3;
   Help3=new Matrix(system->Nsites);
   bool muminus,left;
   int Nmod=Nsteps/10; //no output
   // SE-RP path
   if(SERP)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      evolution.Add_Liouville(liouPhon);
      // define the second and third dipole multiplications:
      muminus=false;	// since I need mu_plus
      left=true;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas RP_SE2(*system, *OpenCLinfo, muminus,left,0,Ntuples);
      evolution.Add_ManipulatingObservable(RP_SE2);
      muminus=false;
      left=false;
      Return_Witness witness_SERP(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_SERP);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_SERP_t[i]+=double_complex(witness_SERP.tracer[i], witness_SERP.tracei[i]);
         data_sum_t[i]+=double_complex(witness_SERP.tracer[i], witness_SERP.tracei[i]);
      }
      evolution.freeDeviceMemory();
      fprintf(stdout,"# SERP method=HEOM done\n");
   }
   // SE-NR path
   if(SENR)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, system->mup, rho0, system->Nsites);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      evolution.Add_Liouville(liouPhon);
      // define the second and third dipole multiplications:
      muminus=true;	// since I need muplus
      left=false;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas NR_SE2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(NR_SE2);
      muminus=false;
      left=false;
      Return_Witness witness_SENR(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_SENR);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_SENR_t[i]+=double_complex(witness_SENR.tracer[i], witness_SENR.tracei[i]);
         data_sum_t[i]+=double_complex(witness_SENR.tracer[i], witness_SENR.tracei[i]);
      }
      evolution.freeDeviceMemory();
      fprintf(stdout,"# SENR method=HEOM done\n");
   }
   // GB-RP path
   if (GBRP)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      evolution.Add_Liouville(liouPhon);
      muminus=false;	// since I need muplus
      left=false;	// since I need to multiply from the right
      Multiply_Dipole_Sigmas RP_GB2(*system, *OpenCLinfo, muminus,left,0,Ntuples);
      evolution.Add_ManipulatingObservable(RP_GB2);
      muminus=false;
      left=true;
      Return_Witness witness_GBRP(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_GBRP);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_GBRP_t[i]+=double_complex(witness_GBRP.tracer[i], witness_GBRP.tracei[i]);
         data_sum_t[i]+=double_complex(witness_GBRP.tracer[i], witness_GBRP.tracei[i]);
      }
      evolution.freeDeviceMemory();
      fprintf(stdout,"# GBRP method=HEOM done\n");
   }
   // GB-NR path
   if (GBNR)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, system->mup, rho0, system->Nsites);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      evolution.Add_Liouville(liouPhon);
      muminus=true;	// since I need muplus
      left=true;	// since I need to multiply from the right
      Multiply_Dipole_Sigmas NR_GB2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(NR_GB2);
      muminus=false;
      left=true;
      Return_Witness witness_GBNR(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_GBNR);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_GBNR_t[i]+=double_complex(witness_GBNR.tracer[i], witness_GBNR.tracei[i]);
         data_sum_t[i]+=double_complex(witness_GBNR.tracer[i], witness_GBNR.tracei[i]);
      }
      evolution.freeDeviceMemory();
      fprintf(stdout,"# GBNR method=HEOM done\n");
   }
   liouPhon.freeDeviceMemory();
   // skip ESA stuff if not needed, no need to allocate unnecessary memory of double-excitons
   if (ESANR==0&&ESARP==0)
   {
      return;
   } // no ESA terms

   /* ESA PREPARATION */
   system->double_Excitation();
   LiouvilleHExciton liouHexciDouble(*system, *OpenCLinfo, Ntuples);
   LiouvillePhononDoubleExcitonLowTemp  liouPhonDouble(*system, *OpenCLinfo, *sigma_Host, Nmax);

   for(int i=0; i<sd_CoupledSite.size(); i++)
   {
      //handle original DL spectral density
      if(sd_Npeaks[i]==1 and sd_gammapos_invcm[i][0]==0.0)
      {
         lambda=new double[sd_Npeaks[i]];
         gamma=new double_complex[sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15), 0.);
         }
         // NOTE since 1d spectra get the zero exciton state |0_ex><0_ex|, I need the "+1" in sd_CoupledSite[i]+1
         liouPhonDouble.AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
      }
      else
      {
         lambda=new double[2*sd_Npeaks[i]];
         gamma=new double_complex[2*sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[2*k]=sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            lambda[2*k+1]=sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            gamma[2*k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),+sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
            gamma[2*k+1]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),-sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
         }
         // NOTE since 1d spectra get the zero exciton state |0_ex><0_ex|, I need the "+1" in sd_CoupledSite[i]+1
         liouPhonDouble.AddSite(sd_CoupledSite[i]+1,2*sd_Npeaks[i], lambda, gamma);
      }
   }
   fprintf(stdout, "# initialize coupling to the phonon-bath ...\n");
   liouPhonDouble.initialize();
   fprintf(stdout, "# ... done\n");

   rho0=new Matrix(system->Nsites);
   rhoInit=new Matrix(system->Nsites); //in Site Basis
   Help3=new Matrix(system->Nsites);
   // ESA-RP path
   if (ESARP)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexciDouble);
      evolution.Add_Liouville(liouPhonDouble);
      muminus=false;	// since I need muplus
      left=true;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas RP_ESA2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(RP_ESA2);
      muminus=false;
      left=true;
      Return_Witness witness_ESARP(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_ESARP);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      ///ESA goes with "-"
      for(int i=0; i<Ndata;  i++)
      {
         data_ESARP_t[i]-=double_complex(witness_ESARP.tracer[i], witness_ESARP.tracei[i]);
         data_sum_t[i]-=double_complex(witness_ESARP.tracer[i], witness_ESARP.tracei[i]);
      }
      evolution.freeDeviceMemory();
      fprintf(stdout,"# ESARP method=HEOM done\n");
   }
   // ESA-NR path
   if (ESANR)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit,system->mup, rho0,  system->Nsites);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexciDouble);
      evolution.Add_Liouville(liouPhonDouble);
      muminus=true;	// since I need muplus
      left=false;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas NR_ESA2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(NR_ESA2);
      muminus=false;
      left=true;
      Return_Witness witness_ESANR(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_ESANR);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      /// ESA goes with "-"
      for(int i=0; i<Ndata;  i++)
      {
         data_ESANR_t[i]-=double_complex(witness_ESANR.tracer[i], witness_ESANR.tracei[i]);
         data_sum_t[i]-=double_complex(witness_ESANR.tracer[i], witness_ESANR.tracei[i]);
      }
      evolution.freeDeviceMemory();
      fprintf(stdout,"# ESANR method=HEOM done\n");
   }
   liouPhonDouble.freeDeviceMemory();
}


void calculate_PumpProbeWitness::execute_calculate_PumpProbeWitness_rotav_10shots(const char* method)
{
   for(int shot=1; shot<=10; shot++)
   {
      fprintf(stdout,"#\n# ten shot rotational average: run shot %i\n",shot);
      double Pol[3];

      /// for the 10 shoot model:
      double phi=(1.+sqrt(5.))/2.;
      double norm=1./sqrt(3.);
      if (shot==1)
      {
         Pol[0]=1.*norm;
         Pol[1]=1.*norm;
         Pol[2]=1.*norm;
      }
      else if (shot==2)
      {
         Pol[0]=1.*norm;
         Pol[1]=-1.*norm;
         Pol[2]=1.*norm;
      }
      else if (shot==3)
      {
         Pol[0]=-1.*norm;
         Pol[1]=1.*norm;
         Pol[2]=1.*norm;
      }
      else if (shot==4)
      {
         Pol[0]=-1.*norm;
         Pol[1]=-1.*norm;
         Pol[2]=1.*norm;
      }
      else if (shot==5)
      {
         Pol[0]=0;
         Pol[1]=1./phi*norm;
         Pol[2]=phi*norm;
      }
      else if (shot==6)
      {
         Pol[0]=0.;
         Pol[1]=-1./phi*norm;
         Pol[2]=phi*norm;
      }
      else if (shot==7)
      {
         Pol[0]=1./phi*norm;
         Pol[1]=phi*norm;
         Pol[2]=0.;
      }
      else if (shot==8)
      {
         Pol[0]=-1./phi*norm;
         Pol[1]=phi*norm;
         Pol[2]=0;
      }
      else if (shot==9)
      {
         Pol[0]=phi*norm;
         Pol[1]=0;
         Pol[2]=1./phi*norm;
      }
      else if (shot==10)
      {
         Pol[0]=-phi*norm;
         Pol[1]=0;
         Pol[2]=1./phi*norm;
      }
      // set polarization in system
      for(int i=0; i<system->Nsites1; i++)
      {
         system->DipoleOperator[i+1]=Pol[0]*system->DipoleMoments[i*3+0]+Pol[1]*system->DipoleMoments[i*3+1]+Pol[2]*system->DipoleMoments[i*3+2];
      }
      if(strncmp(method,"HEOM",strlen(method))==0)
      {
         execute_calculate_PumpProbeWitness();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_PumpProbeWitness(method);
      }
   }
   write_data("ten shot rotational average");
}

void calculate_PumpProbeWitness::execute_calculate_PumpProbeWitness_rotav_4shots(const char* method)
{
   for(int shot=1; shot<=4; shot++)
   {
      fprintf(stdout,"#\n# four shot rotational average: run shot %i\n",shot);
      double Pol[3];

      if (shot==1)
      {
         Pol[0]=+1.0/sqrt(3.0);
         Pol[1]=+1.0/sqrt(3.0);
         Pol[2]=+1.0/sqrt(3.0);
      }
      else if (shot==2)
      {
         Pol[0]=+1.0/sqrt(3.0);
         Pol[1]=+1.0/sqrt(3.0);
         Pol[2]=-1.0/sqrt(3.0);
      }
      else if (shot==3)
      {
         Pol[0]=+1.0/sqrt(3.0);
         Pol[1]=-1.0/sqrt(3.0);
         Pol[2]=-1.0/sqrt(3.0);
      }
      else if (shot==4)
      {
         Pol[0]=+1.0/sqrt(3.0);
         Pol[1]=-1.0/sqrt(3.0);
         Pol[2]=+1.0/sqrt(3.0);
      }

      // set polarization in system
      for(int i=0; i<system->Nsites1; i++)
      {
         system->DipoleOperator[i+1]=Pol[0]*system->DipoleMoments[i*3+0]+Pol[1]*system->DipoleMoments[i*3+1]+Pol[2]*system->DipoleMoments[i*3+2];
      }
      if(strncmp(method,"HEOM",strlen(method))==0)
      {
         execute_calculate_PumpProbeWitness();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_PumpProbeWitness(method);
      }
   }
   write_data("four shot rotational average");
}


void calculate_PumpProbeWitness::execute_calculate_Redfield_PumpProbeWitness(const char* method)
{
   system->create_mum_mup();
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      fprintf(stdout, "# initialize propagation method=full Redfield ...\n");
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      fprintf(stdout, "# initialize propagation method=secular Redfield ...\n");
   }
   Nmax=0; // for full Redfield only one denstity operator needs to be propagated -> set Nmax=0, then we can use same propagation class as for HEOM
   sigma_Host=new SigmaHelpMatrices(system->Nsites, Nmax, system->Ncoupled); // kind of dummy, required to use the same propagation class as for HEOM
   Ntuples=sigma_Host->AllsigmaTuples.size(); // now per construcution Ntuples=1
   fprintf(stdout, "# ... done\n");
   LiouvilleHExcitonRedfield liouHexci(*system,  *OpenCLinfo);
   LiouvillePhonon_fullRedfield* liouPhonfull;
   LiouvillePhonon_secularRedfield* liouPhonsecular;
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      liouPhonfull=new LiouvillePhonon_fullRedfield(*system, *OpenCLinfo);
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      liouPhonsecular=new LiouvillePhonon_secularRedfield(*system, *OpenCLinfo);
   }
   for(int i=0; i<sd_CoupledSite.size(); i++)
   {
      //handle original DL spectral density
      if(sd_Npeaks[i]==1 and sd_gammapos_invcm[i][0]==0.0)
      {
         lambda=new double[sd_Npeaks[i]];
         gamma=new double_complex[sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15), 0.);
         }
         plot_spectral_density(sd_CoupledSite[i], sd_Npeaks[i], lambda, gamma, fn_output, Ediffmax_invcm);
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            liouPhonfull->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            liouPhonsecular->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
         }
      }
      else
      {
         lambda=new double[sd_Npeaks[i]];
         gamma=new double_complex[sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),+sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
         }
         plot_spectral_densityRedfield(sd_CoupledSite[i], sd_Npeaks[i], lambda, gamma, fn_output, Ediffmax_invcm);
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            liouPhonfull->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            liouPhonsecular->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
         }
      }
   }
   fprintf(stdout, "# initialize coupling to the phonon-bath ...\n");
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      liouPhonfull->initialize();
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      liouPhonsecular->initialize();
   }
   fprintf(stdout, "# ... done\n");
   fflush(stdout);
   Matrix *rho0;
   rho0=new Matrix(system->Nsites);
   Matrix *rhoInit;
   rhoInit=new Matrix(system->Nsites); //in Site Basis
   Matrix *Help3;
   Help3=new Matrix(system->Nsites);
   bool muminus,left;
   int Nmod=Nsteps/10;
   // SE-RP path
   if(SERP)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
      system->convertEigenbasis(rhoInit);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonfull);
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonsecular);
      }
      // define the second and third dipole multiplications:
      muminus=false;	// since I need mu_plus
      left=true;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas_Redfield RP_SE2(*system, *OpenCLinfo, muminus,left,0,Ntuples);
      evolution.Add_ManipulatingObservable(RP_SE2);
      muminus=false;
      left=false;
      Return_Witness_Redfield witness_SERP(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_SERP);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_SERP_t[i]+=double_complex(witness_SERP.tracer[i], witness_SERP.tracei[i]);
         data_sum_t[i]+=double_complex(witness_SERP.tracer[i], witness_SERP.tracei[i]);
      }
      evolution.freeDeviceMemory();
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# SERP method=full Redfield done\n");
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# SERP method=secular Redfield done\n");
      }
      fflush(stdout);
   }
   // SE-NR path
   if(SENR)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, system->mup, rho0, system->Nsites);
      system->convertEigenbasis(rhoInit);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonfull);
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonsecular);
      }
      // define the second and third dipole multiplications:
      muminus=true;	// since I need muplus
      left=false;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas_Redfield NR_SE2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(NR_SE2);
      muminus=false;
      left=false;
      Return_Witness_Redfield witness_SENR(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_SENR);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_SENR_t[i]+=double_complex(witness_SENR.tracer[i], witness_SENR.tracei[i]);
         data_sum_t[i]+=double_complex(witness_SENR.tracer[i], witness_SENR.tracei[i]);
      }
      evolution.freeDeviceMemory();
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# SENR method=full Redfield done\n");
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# SENR method=secular Redfield done\n");
      }
      fflush(stdout);
   }
   // GB-RP path
   if (GBRP)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
      system->convertEigenbasis(rhoInit);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonfull);
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonsecular);
      }
      muminus=false;	// since I need muplus
      left=false;	// since I need to multiply from the right
      Multiply_Dipole_Sigmas_Redfield RP_GB2(*system, *OpenCLinfo, muminus,left,0,Ntuples);
      evolution.Add_ManipulatingObservable(RP_GB2);
      muminus=false;
      left=true;
      Return_Witness_Redfield witness_GBRP(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_GBRP);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_GBRP_t[i]+=double_complex(witness_GBRP.tracer[i], witness_GBRP.tracei[i]);
         data_sum_t[i]+=double_complex(witness_GBRP.tracer[i], witness_GBRP.tracei[i]);
      }
      evolution.freeDeviceMemory();
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# GBRP method=full Redfield done\n");
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# GBRP method=secular Redfield done\n");
      }
      fflush(stdout);
   }
   // GB-NR path
   if (GBNR)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, system->mup, rho0, system->Nsites);
      system->convertEigenbasis(rhoInit);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexci);
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonfull);
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonsecular);
      }
      muminus=true;	// since I need muplus
      left=true;	// since I need to multiply from the right
      Multiply_Dipole_Sigmas_Redfield NR_GB2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(NR_GB2);
      muminus=false;
      left=true;
      Return_Witness_Redfield witness_GBNR(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_GBNR);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      for(int i=0; i<Ndata;  i++)
      {
         data_GBNR_t[i]+=double_complex(witness_GBNR.tracer[i], witness_GBNR.tracei[i]);
         data_sum_t[i]+=double_complex(witness_GBNR.tracer[i], witness_GBNR.tracei[i]);
      }
      evolution.freeDeviceMemory();
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# GBNR method=full Redfield done\n");
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# GBNR method=secular Redfield done\n");
      }
      fflush(stdout);
   }
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      liouPhonfull->freeDeviceMemory();
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      liouPhonsecular->freeDeviceMemory();
   }
   // skip ESA stuff if not needed, no need to allocate unnecessary memory of double-excitons
   if (ESANR==0&&ESARP==0)
   {
      return;
   } // no ESA terms

   /* ESA PREPARATION */
   system->double_Excitation();
   system->computeEigenvects_Eigenvals_doubleExciton();
   LiouvilleHExcitonRedfield liouHexciDouble(*system, *OpenCLinfo);
   LiouvillePhononDouble_fullRedfield* liouPhonDoublefull;
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      liouPhonDoublefull=new LiouvillePhononDouble_fullRedfield(*system, *OpenCLinfo);
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      fprintf(stdout,"# ERROR secular Redfield can't be used together with double excitons.\n");
      fprintf(stdout,"# ERROR ESANR or ESARP do not support secular Redfield. ABORT\n");
      exit(1);
   }
   for(int i=0; i<sd_CoupledSite.size(); i++)
   {
      //handle original DL spectral density
      if(sd_Npeaks[i]==1 and sd_gammapos_invcm[i][0]==0.0)
      {
         lambda=new double[sd_Npeaks[i]];
         gamma=new double_complex[sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15), 0.);
         }
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            liouPhonDoublefull->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
         }
//         if(strncmp(method,"secular Redfield",strlen(method))==0)
//         {
//            liouPhonDoublesecular->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
//         }
      }
      else
      {
         lambda=new double[sd_Npeaks[i]];
         gamma=new double_complex[sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),+sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
         }
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            liouPhonDoublefull->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
         }
//         if(strncmp(method,"secular Redfield",strlen(method))==0)
//         {
//            liouPhonDoublesecular->AddSite(sd_CoupledSite[i]+1,sd_Npeaks[i], lambda, gamma);
//         }
      }
   }
   fprintf(stdout, "# initialize coupling to the phonon-bath ...\n");
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      liouPhonDoublefull->initialize();
   }
//   if(strncmp(method,"secular Redfield",strlen(method))==0)
//   {
//      liouPhonDoublesecular->initialize();
//   }
   fprintf(stdout, "# ... done\n");
   rho0=new Matrix(system->Nsites);
   rhoInit=new Matrix(system->Nsites); //in Site Basis
   Help3=new Matrix(system->Nsites);
   // ESA-RP path
   if (ESARP)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
      system->convertEigenbasis(rhoInit);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexciDouble);
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonDoublefull);
      }
//         if(strncmp(method,"secular Redfield",strlen(method))==0)
//         {
//            evolution.Add_Liouville(*liouPhonDoublesecular);
//         }
      muminus=false;	// since I need muplus
      left=true;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas_Redfield RP_ESA2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(RP_ESA2);
      muminus=false;
      left=true;
      Return_Witness_Redfield witness_ESARP(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_ESARP);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      ///ESA goes with "-"
      for(int i=0; i<Ndata;  i++)
      {
         data_ESARP_t[i]-=double_complex(witness_ESARP.tracer[i], witness_ESARP.tracei[i]);
         data_sum_t[i]-=double_complex(witness_ESARP.tracer[i], witness_ESARP.tracei[i]);
      }
      evolution.freeDeviceMemory();
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# ESARP method=full Redfield done\n");
      }
//            if(strncmp(method,"secular Redfield",strlen(method))==0)
//            {
//               fprintf(stdout,"# ESARP method=secular Redfield done\n");
//            }
      fflush(stdout);
   }
   // ESA-NR path
   if (ESANR)
   {
      Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
      rho0->assignEntry(0,0,1.0);
      MatrixMuliplication(rhoInit,system->mup, rho0,  system->Nsites);
      system->convertEigenbasis(rhoInit);
      evolution.initialize(*rhoInit);
      evolution.Add_Liouville(liouHexciDouble);
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         evolution.Add_Liouville(*liouPhonDoublefull);
      }
//         if(strncmp(method,"secular Redfield",strlen(method))==0)
//         {
//            evolution.Add_Liouville(*liouPhonDoublesecular);
//         }
      muminus=true;	// since I need muplus
      left=false;	// since I need to multiply from the left
      Multiply_Dipole_Sigmas_Redfield NR_ESA2(*system, *OpenCLinfo, muminus,left,0, Ntuples);
      evolution.Add_ManipulatingObservable(NR_ESA2);
      muminus=false;
      left=true;
      Return_Witness_Redfield witness_ESANR(*system, *OpenCLinfo, DNT,  muminus, left, Ntuples, system->mum);
      evolution.Add_ManipulatingObservable(witness_ESANR);
      evolution.propagate(Nmod); // Nmod defines the time steps for an output
      /// ESA goes with "-"
      for(int i=0; i<Ndata;  i++)
      {
         data_ESANR_t[i]-=double_complex(witness_ESANR.tracer[i], witness_ESANR.tracei[i]);
         data_sum_t[i]-=double_complex(witness_ESANR.tracer[i], witness_ESANR.tracei[i]);
      }
      evolution.freeDeviceMemory();
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# ESANR method=full Redfield done\n");
      }
//            if(strncmp(method,"secular Redfield",strlen(method))==0)
//            {
//               fprintf(stdout,"# ESANR method=secular Redfield done\n");
//            }
      fflush(stdout);
   }
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      liouPhonDoublefull->freeDeviceMemory();
   }
//   if(strncmp(method,"secular Redfield",strlen(method))==0)
//   {
//      liouPhonDoublesecular->freeDeviceMemory();
//   }
}




#endif
