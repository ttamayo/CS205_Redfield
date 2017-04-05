/*
 * Run_2dSpectra.h
 *
 *  Created on: Jan 24, 2014
 *      Author: christoph
 */

#ifndef RUN_2DSPECTRA_H_
#define RUN_2DSPECTRA_H_


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "headers.h"
#include <iostream>

class calculate_2d
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Nmax;
   double dt;
   int DNT;
   int Nsteps3;
   int delayTime;
   double** sd_lambda_invcm;
   double** sd_gamma_fs;
   double** sd_gammapos_invcm;
   vector<int> sd_CoupledSite;
   vector<int> sd_Npeaks;
   int padding;
   const char *fn_output;
   int GBRP;
   int GBNR;
   int SERP;
   int SENR;
   int ESARP;
   int ESANR;
   double Emin;
   double Emax;
   bool specify_Enrange;
private:
   SigmaHelpMatrices *sigma_Host;
   int Ntuples;
   double* lambda;
   double_complex* gamma;
   double Ediffmax_invcm;
   double w1_invcm;
   double w3_invcm;
   int NTP1;
   int NTP3;
   fftw_complex *data_SERP_t;
   fftw_complex *data_SENR_t;
   fftw_complex *data_GBRP_t;
   fftw_complex *data_GBNR_t;
   fftw_complex *data_ESARP_t;
   fftw_complex *data_ESANR_t;
   fftw_complex *data_SERP_f;
   fftw_complex *data_SENR_f;
   fftw_complex *data_GBRP_f;
   fftw_complex *data_GBNR_f;
   fftw_complex *data_ESARP_f;
   fftw_complex *data_ESANR_f;
   fftw_complex *data_sum_f;
   char Fileheader_calculatedPath[255];
public:
   calculate_2d(System& INsystem,
                OpenCL_Init& INOpenCLinfo,
                int INNmax,
                double INdt,
                int INDNT,
                int INNsteps3,
                int INdelayTime,
                double** INsd_lambda_invcm,
                double** INsd_gamma_fs,
                double** INsd_gammapos_invcm,
                vector<int> INsd_CoupledSite,
                vector<int> INsd_Npeaks,
                int INpadding,
                const char* INfn_output,
                int INGBRP,
                int INGBNR,
                int INSERP,
                int INSENR,
                int INESARP,
                int INESANR,
                double INEmin,
                double INEmax,
                bool INspecify_Enrange
               );
   void execute_calculate_2d();
   void execute_calculate_Redfield_2d(const char* method);
   void execute_calculate_2d(const char* dipole_moments,const char* method); //FIXME here depending on method switch to execute_calculate_Redfied_2d
   void execute_calculate_2d_rotav_4shots(const char* method); //FIXME here depending on method switch to execute_calculate_Redfied_2d
   void execute_calculate_2d_rotav_10shots(const char* method); //FIXME here depending on method switch to execute_calculate_Redfied_2d
   void print_2d_gnuplot(int delayTime, const char* fn_output);
   void write_data(const char* info_shotdipole);
   void perform_fftRP(fftw_complex* data_t, fftw_complex* data_f);
   void perform_fftNR(fftw_complex* data_t, fftw_complex* data_f);
};


calculate_2d::calculate_2d(System& INsystem,
                           OpenCL_Init& INOpenCLinfo,
                           int INNmax,
                           double INdt,
                           int INDNT,
                           int INNsteps3,
                           int INdelayTime,
                           double** INsd_lambda_invcm,
                           double** INsd_gamma_fs,
                           double** INsd_gammapos_invcm,
                           vector<int> INsd_CoupledSite,
                           vector<int> INsd_Npeaks,
                           int INpadding,
                           const char* INfn_output,
                           int INGBRP,
                           int INGBNR,
                           int INSERP,
                           int INSENR,
                           int INESARP,
                           int INESANR,
                           double INEmin,
                           double INEmax,
                           bool INspecify_Enrange
                          )  :
   system(&INsystem),
   OpenCLinfo(&INOpenCLinfo),
   Nmax(INNmax),
   dt(INdt),
   DNT(INDNT),
   Nsteps3(INNsteps3),
   delayTime(INdelayTime),
   sd_lambda_invcm(INsd_lambda_invcm),
   sd_gamma_fs(INsd_gamma_fs),
   sd_gammapos_invcm(INsd_gammapos_invcm),
   sd_CoupledSite(INsd_CoupledSite),
   sd_Npeaks(INsd_Npeaks),
   padding(INpadding),
   fn_output(INfn_output),
   GBRP(INGBRP),
   GBNR(INGBNR),
   SERP(INSERP),
   SENR(INSENR),
   ESARP(INESARP),
   ESANR(INESANR),
   Emin(INEmin),
   Emax(INEmax),
   specify_Enrange(INspecify_Enrange)
{
   fprintf(stdout, "#\n# initialize parameter for 2d-echo spectrum\n");
   sprintf(Fileheader_calculatedPath,"# output: sum of all specified pathways:");
   // initialize data for time lines
   NTP1=Nsteps3/DNT*padding; // steps along t1
   NTP3=Nsteps3/DNT*padding; // steps along t3
   int sizeN=NTP1*NTP3*sizeof(fftw_complex);
   if(SERP==1)
   {
      sprintf(Fileheader_calculatedPath,"%s SE-RP,",Fileheader_calculatedPath);
      data_SERP_t=(fftw_complex*)fftw_malloc(sizeN);
      data_SERP_f=(fftw_complex*)fftw_malloc(sizeN);
      for(int i=0; i<NTP1*NTP3; i++)
      {
         data_SERP_t[i][0]=0.0;
         data_SERP_t[i][1]=0.0;
         data_SERP_f[i][0]=0.0;
         data_SERP_f[i][1]=0.0;
      }
   }
   if(SENR==1)
   {
      sprintf(Fileheader_calculatedPath,"%s SE-NR,",Fileheader_calculatedPath);
      data_SENR_t=(fftw_complex*)fftw_malloc(sizeN);
      data_SENR_f=(fftw_complex*)fftw_malloc(sizeN);
      for(int i=0; i<NTP1*NTP3; i++)
      {
         data_SENR_t[i][0]=0.0;
         data_SENR_t[i][1]=0.0;
         data_SENR_f[i][0]=0.0;
         data_SENR_f[i][1]=0.0;
      }
   }
   if(GBRP==1)
   {
      sprintf(Fileheader_calculatedPath,"%s GB-RP,",Fileheader_calculatedPath);
      data_GBRP_t=(fftw_complex*)fftw_malloc(sizeN);
      data_GBRP_f=(fftw_complex*)fftw_malloc(sizeN);
      for(int i=0; i<NTP1*NTP3; i++)
      {
         data_GBRP_t[i][0]=0.0;
         data_GBRP_t[i][1]=0.0;
         data_GBRP_f[i][0]=0.0;
         data_GBRP_f[i][1]=0.0;
      }
   }
   if(GBNR==1)
   {
      sprintf(Fileheader_calculatedPath,"%s GB-NR,",Fileheader_calculatedPath);
      data_GBNR_t=(fftw_complex*)fftw_malloc(sizeN);
      data_GBNR_f=(fftw_complex*)fftw_malloc(sizeN);
      for(int i=0; i<NTP1*NTP3; i++)
      {
         data_GBNR_t[i][0]=0.0;
         data_GBNR_t[i][1]=0.0;
         data_GBNR_f[i][0]=0.0;
         data_GBNR_f[i][1]=0.0;
      }
   }
   if(ESARP==1)
   {
      sprintf(Fileheader_calculatedPath,"%s ESA-RP,",Fileheader_calculatedPath);
      data_ESARP_t=(fftw_complex*)fftw_malloc(sizeN);
      data_ESARP_f=(fftw_complex*)fftw_malloc(sizeN);
      for(int i=0; i<NTP1*NTP3; i++)
      {
         data_ESARP_t[i][0]=0.0;
         data_ESARP_t[i][1]=0.0;
         data_ESARP_f[i][0]=0.0;
         data_ESARP_f[i][1]=0.0;
      }
   }
   if(ESANR==1)
   {
      sprintf(Fileheader_calculatedPath,"%s ESA-NR,",Fileheader_calculatedPath);
      data_ESANR_t=(fftw_complex*)fftw_malloc(sizeN);
      data_ESANR_f=(fftw_complex*)fftw_malloc(sizeN);
      for(int i=0; i<NTP1*NTP3; i++)
      {
         data_ESANR_t[i][0]=0.0;
         data_ESANR_t[i][1]=0.0;
         data_ESANR_f[i][0]=0.0;
         data_ESANR_f[i][1]=0.0;
      }
   }
   data_sum_f=(fftw_complex*)fftw_malloc(sizeN);
   for(int i=0; i<NTP1*NTP3; i++)
   {
      data_sum_f[i][0]=0.0;
      data_sum_f[i][1]=0.0;
   }
   //
   // specify maximal difference between to exction eigen energies to specify plot range of the spectral density (plot energy range where differences in exciton
   // energies live (extend by factor of 1.5 to be sure that important part is plotted))
   system->computeEigenvects_Eigenvals();
   if(specify_Enrange==false)
   {
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
      Emax=system->Eigenvals[system->Nsites1-1]/const_invcmtomeV+Ediffmax_invcm;
      Emin=system->Eigenvals[0]/const_invcmtomeV-Ediffmax_invcm;
   }
   //max energy range
   double Eminall=-NTP1*0.5*const_hbar/const_invcmtomeV*2.*M_PI/((NTP1-1)*dt);
   double Emaxall=+NTP1*0.5*const_hbar/const_invcmtomeV*2.*M_PI/((NTP1-1)*dt);
   if(specify_Enrange==true)
   {
      // check if specified Erange fits into [Eminall, Emaxall]
      if(Emin<Eminall)
      {
         fprintf(stderr,"error: specified En_range Emin=%e must be larger than %e (constraint by FFT grid)\n",Emin, Eminall);
         fprintf(stderr,"change either Emin or reduce the time step dt\n");
         exit(1);
      }
      if(Emax>Emaxall)
      {
         fprintf(stderr,"error: specified En_range Emax=%e must be smaller than %e (constraint by FFT grid)\n",Emax, Emaxall);
         fprintf(stderr,"change either Emax or reduce the time step dt\n");
         exit(1);
      }
   }
}

void calculate_2d::print_2d_gnuplot(int delayTime, const char* fn_output)
{
   char postfix[255];
   char fn[500];
   char submit_command[500];
   // sum of all spectra:

   strcpy(postfix,"sum");
   sprintf(fn,"%s_2d-echo_delay%d_%s_freq",fn_output,delayTime,postfix);
   sprintf(submit_command,"generate_2d_freq_with_gnuplot.sh %s",fn);
   std::system(submit_command);
   fprintf(stdout,"#\n# plot 2d echo spectrum\n");
   fprintf(stdout,"# SUBMIT: %s\n",submit_command);
   int err;
   err=std::system(submit_command);
   if (err!=0)
   {
      fprintf(stdout,"WARNING: std::system(%s) failed (err=%d).\n",submit_command,err);
//      exit(err);
   }
   std::system("rm cont.dat");
}

void calculate_2d::execute_calculate_2d(const char* dipole_moments, const char* method)
{
   fprintf(stdout,"# specified dipole moments: %s\n",dipole_moments);
   if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
   {
      if(strncmp(method,"HEOM",strlen(method))==0)
      {
         execute_calculate_2d();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_2d(method);
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
         execute_calculate_2d();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_2d(method);
      }
   }
   else
   {
      fprintf(stderr,"effective dipole moments not specified correctly. ABORT.\n");
      exit(1);
   }
   write_data(dipole_moments);
}


void calculate_2d::perform_fftRP(fftw_complex* data_t, fftw_complex* data_f)
{
   // RP: exp(-I*t1,+I*t3) -> two 1d transform
   fftw_plan fftplan1stRP; // for rephasing FFT along t1
   fftw_plan fftplan2ndRP; // for rephasing FFT along t3
   // 1st data_t -> data_f
   fftplan1stRP=fftw_plan_many_dft(1,&NTP3,NTP1,
                                   (fftw_complex*)data_t,NULL,1,NTP3,
                                   (fftw_complex*)data_f,NULL,1,NTP3,
                                   FFTW_FORWARD, FFTW_ESTIMATE);
   // 2nd data_f -> data_f
   fftplan2ndRP=fftw_plan_many_dft(1,&NTP1,NTP3,
                                   (fftw_complex*)data_f,NULL,NTP3,1,
                                   (fftw_complex*)data_f,NULL,NTP3,1,
                                   FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(fftplan1stRP);
   fftw_execute(fftplan2ndRP);
}
void calculate_2d::perform_fftNR(fftw_complex* data_t, fftw_complex* data_f)
{
   fftw_plan fftplan2dNR; // for non-rephasing FFT along both t1 and t3
   // NR: exp(+I*t1,+I*t3) -> FFTW_BACKWARD
   fftplan2dNR=fftw_plan_dft_2d(NTP1,NTP3,(fftw_complex*)data_t,(fftw_complex*)data_f,FFTW_BACKWARD,FFTW_ESTIMATE);
   fftw_execute(fftplan2dNR);
}
void calculate_2d::write_data(const char* info_shotdipole)
{
   char fn[500];
   FILE *fo;
   time_t t;
   // FIXME: in the end w1 and w3 seem to be swapped ??? ugly hack to write out transpose: FIXME: same for (t1,t3) timeline!
   if(SERP==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_2d-echo_delay%d_SERP.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway SE-RP, time signal\n#\n");
      fprintf(fo,"#  t1 in s     t2 in s     Re[signal]    Im[signal]\n");
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            fprintf(fo, "%e %e %e %e\n",NstepsOverDNT*DNT*dt, t1*dt*DNT, data_SERP_t[t1*NTP3+NstepsOverDNT][0], data_SERP_t[t1*NTP3+NstepsOverDNT][1]);
         }
         fprintf(fo, "\n");
      }
      fclose(fo);
      // perform FFT
      for(int t1=0; t1<NTP1; t1++)
      {
         data_SERP_t[t1*NTP3+0][0]=0.5*(data_SERP_t[t1*NTP3+0][0]+data_SERP_t[t1*NTP3+(NTP3-1)][0]);
         data_SERP_t[t1*NTP3+0][1]=0.5*(data_SERP_t[t1*NTP3+0][1]+data_SERP_t[t1*NTP3+(NTP3-1)][1]);
      }
      for(int t3=0; t3<NTP3; t3++)
      {
         data_SERP_t[0*NTP3+t3][0]=0.5*(data_SERP_t[0*NTP3+t3][0]+data_SERP_t[(NTP1-1)*NTP3+t3][0]);
         data_SERP_t[0*NTP3+t3][1]=0.5*(data_SERP_t[0*NTP3+t3][1]+data_SERP_t[(NTP1-1)*NTP3+t3][1]);
      }
      perform_fftRP(data_SERP_t, data_SERP_f);
      sprintf(fn,"%s_2d-echo_delay%d_SERP_freq.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway SE-RP, spectrum\n#\n");
      fprintf(fo,"# w1 in cm^-1  w2 in cm^-1  Re[signal]  Im[signal]\n");
      for(int w1=0; w1<NTP1; w1++)
      {
         bool enter=false;
         for(int w3=0; w3<NTP3; w3++)
         {
            double w1_invcm,w3_invcm;
            int i1=w1,i3=w3;
            if (w1<NTP1/2) i1+=NTP1/2; // w1=0..NTP1/2-1 gives the negative freq
            if (w3<NTP3/2) i3+=NTP3/2;
            if (w1>=NTP1/2) i1-=NTP1/2; // w1=NTP1/2..NTP1 gives the positive freq
            if (w3>=NTP3/2) i3-=NTP3/2;
            w1_invcm=double(w1-NTP1/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP1-1)*DNT*dt);
            w3_invcm=double(w3-NTP3/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP3-1)*DNT*dt);
            if(w1_invcm>Emin and w1_invcm<Emax and w3_invcm>Emin and w3_invcm<Emax)
            {
               // FIXME ugly hack to swap w1 <-> w3
               enter=true;
               data_sum_f[i1*NTP3+i3][0]+=data_SERP_f[i1*NTP3+i3][0];
               data_sum_f[i1*NTP3+i3][1]+=data_SERP_f[i1*NTP3+i3][1];
               fprintf(fo, "%e %e %e %e\n",w3_invcm,w1_invcm,data_SERP_f[i1*NTP3+i3][0], data_SERP_f[i1*NTP3+i3][1]);
            }
         }
         if(enter==true)
         {
            fprintf(fo, "\n");
         }
      }
      fclose(fo);
      fftw_free(data_SERP_t);
      fftw_free(data_SERP_f);
   }
   //
   if(GBRP==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_2d-echo_delay%d_GBRP.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway GB-RP, time signal\n#\n");
      fprintf(fo,"#  t1 in s     t2 in s     Re[signal]    Im[signal]\n");
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            fprintf(fo, "%e %e %e %e\n",NstepsOverDNT*DNT*dt, t1*dt*DNT, data_GBRP_t[t1*NTP3+NstepsOverDNT][0], data_GBRP_t[t1*NTP3+NstepsOverDNT][1]);
         }
         fprintf(fo, "\n");
      }
      fclose(fo);
      // perform FFT
      for(int t1=0; t1<NTP1; t1++)
      {
         data_GBRP_t[t1*NTP3+0][0]=0.5*(data_GBRP_t[t1*NTP3+0][0]+data_GBRP_t[t1*NTP3+(NTP3-1)][0]);
         data_GBRP_t[t1*NTP3+0][1]=0.5*(data_GBRP_t[t1*NTP3+0][1]+data_GBRP_t[t1*NTP3+(NTP3-1)][1]);
      }
      for(int t3=0; t3<NTP3; t3++)
      {
         data_GBRP_t[0*NTP3+t3][0]=0.5*(data_GBRP_t[0*NTP3+t3][0]+data_GBRP_t[(NTP1-1)*NTP3+t3][0]);
         data_GBRP_t[0*NTP3+t3][1]=0.5*(data_GBRP_t[0*NTP3+t3][1]+data_GBRP_t[(NTP1-1)*NTP3+t3][1]);
      }
      perform_fftRP(data_GBRP_t, data_GBRP_f);
      sprintf(fn,"%s_2d-echo_delay%d_GBRP_freq.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway GB-RP, spectrum\n#\n");
      fprintf(fo,"# w1 in cm^-1  w2 in cm^-1  Re[signal]  Im[signal]\n");
      for(int w1=0; w1<NTP1; w1++)
      {
         bool enter=false;
         for(int w3=0; w3<NTP3; w3++)
         {
            double w1_invcm,w3_invcm;
            int i1=w1,i3=w3;
            if (w1<NTP1/2) i1+=NTP1/2; // w1=0..NTP1/2-1 gives the negative freq
            if (w3<NTP3/2) i3+=NTP3/2;
            if (w1>=NTP1/2) i1-=NTP1/2; // w1=NTP1/2..NTP1 gives the positive freq
            if (w3>=NTP3/2) i3-=NTP3/2;
            w1_invcm=double(w1-NTP1/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP1-1)*DNT*dt);
            w3_invcm=double(w3-NTP3/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP3-1)*DNT*dt);
            if(w1_invcm>Emin and w1_invcm<Emax and w3_invcm>Emin and w3_invcm<Emax)
            {
               // FIXME ugly hack to swap w1 <-> w3
               enter=true;
               data_sum_f[i1*NTP3+i3][0]+=data_GBRP_f[i1*NTP3+i3][0];
               data_sum_f[i1*NTP3+i3][1]+=data_GBRP_f[i1*NTP3+i3][1];
               fprintf(fo, "%e %e %e %e\n",w3_invcm,w1_invcm,data_GBRP_f[i1*NTP3+i3][0], data_GBRP_f[i1*NTP3+i3][1]);
            }
         }
         if(enter==true)
         {
            fprintf(fo, "\n");
         }
      }
      fclose(fo);
      fftw_free(data_GBRP_t);
      fftw_free(data_GBRP_f);
   }
   //
   if(ESARP==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_2d-echo_delay%d_ESARP.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway ESA-RP, time signal\n#\n");
      fprintf(fo,"#  t1 in s     t2 in s     Re[signal]    Im[signal]\n");
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            fprintf(fo, "%e %e %e %e\n",NstepsOverDNT*DNT*dt, t1*dt*DNT, data_ESARP_t[t1*NTP3+NstepsOverDNT][0], data_ESARP_t[t1*NTP3+NstepsOverDNT][1]);
         }
         fprintf(fo, "\n");
      }
      fclose(fo);
      // perform FFT
      for(int t1=0; t1<NTP1; t1++)
      {
         data_ESARP_t[t1*NTP3+0][0]=0.5*(data_ESARP_t[t1*NTP3+0][0]+data_ESARP_t[t1*NTP3+(NTP3-1)][0]);
         data_ESARP_t[t1*NTP3+0][1]=0.5*(data_ESARP_t[t1*NTP3+0][1]+data_ESARP_t[t1*NTP3+(NTP3-1)][1]);
      }
      for(int t3=0; t3<NTP3; t3++)
      {
         data_ESARP_t[0*NTP3+t3][0]=0.5*(data_ESARP_t[0*NTP3+t3][0]+data_ESARP_t[(NTP1-1)*NTP3+t3][0]);
         data_ESARP_t[0*NTP3+t3][1]=0.5*(data_ESARP_t[0*NTP3+t3][1]+data_ESARP_t[(NTP1-1)*NTP3+t3][1]);
      }
      perform_fftRP(data_ESARP_t, data_ESARP_f);
      sprintf(fn,"%s_2d-echo_delay%d_ESARP_freq.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway ESA-RP, spectrum\n#\n");
      fprintf(fo,"# w1 in cm^-1  w2 in cm^-1  Re[signal]  Im[signal]\n");
      for(int w1=0; w1<NTP1; w1++)
      {
         bool enter=false;
         for(int w3=0; w3<NTP3; w3++)
         {
            double w1_invcm,w3_invcm;
            int i1=w1,i3=w3;
            if (w1<NTP1/2) i1+=NTP1/2; // w1=0..NTP1/2-1 gives the negative freq
            if (w3<NTP3/2) i3+=NTP3/2;
            if (w1>=NTP1/2) i1-=NTP1/2; // w1=NTP1/2..NTP1 gives the positive freq
            if (w3>=NTP3/2) i3-=NTP3/2;
            w1_invcm=double(w1-NTP1/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP1-1)*DNT*dt);
            w3_invcm=double(w3-NTP3/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP3-1)*DNT*dt);
            if(w1_invcm>Emin and w1_invcm<Emax and w3_invcm>Emin and w3_invcm<Emax)
            {
               // FIXME ugly hack to swap w1 <-> w3
               enter=true;
               data_sum_f[i1*NTP3+i3][0]+=data_ESARP_f[i1*NTP3+i3][0];
               data_sum_f[i1*NTP3+i3][1]+=data_ESARP_f[i1*NTP3+i3][1];
               fprintf(fo, "%e %e %e %e\n",w3_invcm,w1_invcm,data_ESARP_f[i1*NTP3+i3][0], data_ESARP_f[i1*NTP3+i3][1]);
            }
         }
         if(enter==true)
         {
            fprintf(fo, "\n");
         }
      }
      fclose(fo);
      fftw_free(data_ESARP_t);
      fftw_free(data_ESARP_f);
   }

   if(SENR==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_2d-echo_delay%d_SENR.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway SE-NR, time signal\n#\n");
      fprintf(fo,"#  t1 in s     t2 in s     Re[signal]    Im[signal]\n");
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            fprintf(fo, "%e %e %e %e\n",NstepsOverDNT*DNT*dt, t1*dt*DNT, data_SENR_t[t1*NTP3+NstepsOverDNT][0], data_SENR_t[t1*NTP3+NstepsOverDNT][1]);
         }
         fprintf(fo, "\n");
      }
      fclose(fo);
      // perform FFT
      for(int t1=0; t1<NTP1; t1++)
      {
         data_SENR_t[t1*NTP3+0][0]=0.5*(data_SENR_t[t1*NTP3+0][0]+data_SENR_t[t1*NTP3+(NTP3-1)][0]);
         data_SENR_t[t1*NTP3+0][1]=0.5*(data_SENR_t[t1*NTP3+0][1]+data_SENR_t[t1*NTP3+(NTP3-1)][1]);
      }
      for(int t3=0; t3<NTP3; t3++)
      {
         data_SENR_t[0*NTP3+t3][0]=0.5*(data_SENR_t[0*NTP3+t3][0]+data_SENR_t[(NTP1-1)*NTP3+t3][0]);
         data_SENR_t[0*NTP3+t3][1]=0.5*(data_SENR_t[0*NTP3+t3][1]+data_SENR_t[(NTP1-1)*NTP3+t3][1]);
      }
      perform_fftNR(data_SENR_t, data_SENR_f);
      sprintf(fn,"%s_2d-echo_delay%d_SENR_freq.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway SE-NR, spectrum\n#\n");
      fprintf(fo,"# w1 in cm^-1  w2 in cm^-1  Re[signal]  Im[signal]\n");
      for(int w1=0; w1<NTP1; w1++)
      {
         bool enter=false;
         for(int w3=0; w3<NTP3; w3++)
         {
            double w1_invcm,w3_invcm;
            int i1=w1,i3=w3;
            if (w1<NTP1/2) i1+=NTP1/2; // w1=0..NTP1/2-1 gives the negative freq
            if (w3<NTP3/2) i3+=NTP3/2;
            if (w1>=NTP1/2) i1-=NTP1/2; // w1=NTP1/2..NTP1 gives the positive freq
            if (w3>=NTP3/2) i3-=NTP3/2;
            w1_invcm=double(w1-NTP1/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP1-1)*DNT*dt);
            w3_invcm=double(w3-NTP3/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP3-1)*DNT*dt);
            if(w1_invcm>Emin and w1_invcm<Emax and w3_invcm>Emin and w3_invcm<Emax)
            {
               // FIXME ugly hack to swap w1 <-> w3
               enter=true;
               data_sum_f[i1*NTP3+i3][0]+=data_SENR_f[i1*NTP3+i3][0];
               data_sum_f[i1*NTP3+i3][1]+=data_SENR_f[i1*NTP3+i3][1];
               fprintf(fo, "%e %e %e %e\n",w3_invcm,w1_invcm,data_SENR_f[i1*NTP3+i3][0], data_SENR_f[i1*NTP3+i3][1]);
            }
         }
         if(enter==true)
         {
            fprintf(fo, "\n");
         }
      }
      fclose(fo);
      fftw_free(data_SENR_t);
      fftw_free(data_SENR_f);
   }
   //
   if(GBNR==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_2d-echo_delay%d_GBNR.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway GB-NR, time signal\n#\n");
      fprintf(fo,"#  t1 in s     t2 in s     Re[signal]    Im[signal]\n");
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            fprintf(fo, "%e %e %e %e\n",NstepsOverDNT*DNT*dt, t1*dt*DNT, data_GBNR_t[t1*NTP3+NstepsOverDNT][0], data_GBNR_t[t1*NTP3+NstepsOverDNT][1]);
         }
         fprintf(fo, "\n");
      }
      fclose(fo);
      // perform FFT
      for(int t1=0; t1<NTP1; t1++)
      {
         data_GBNR_t[t1*NTP3+0][0]=0.5*(data_GBNR_t[t1*NTP3+0][0]+data_GBNR_t[t1*NTP3+(NTP3-1)][0]);
         data_GBNR_t[t1*NTP3+0][1]=0.5*(data_GBNR_t[t1*NTP3+0][1]+data_GBNR_t[t1*NTP3+(NTP3-1)][1]);
      }
      for(int t3=0; t3<NTP3; t3++)
      {
         data_GBNR_t[0*NTP3+t3][0]=0.5*(data_GBNR_t[0*NTP3+t3][0]+data_GBNR_t[(NTP1-1)*NTP3+t3][0]);
         data_GBNR_t[0*NTP3+t3][1]=0.5*(data_GBNR_t[0*NTP3+t3][1]+data_GBNR_t[(NTP1-1)*NTP3+t3][1]);
      }
      perform_fftNR(data_GBNR_t, data_GBNR_f);
      sprintf(fn,"%s_2d-echo_delay%d_GBNR_freq.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway GB-NR, spectrum\n#\n");
      fprintf(fo,"# w1 in cm^-1  w2 in cm^-1  Re[signal]  Im[signal]\n");
      for(int w1=0; w1<NTP1; w1++)
      {
         bool enter=false;
         for(int w3=0; w3<NTP3; w3++)
         {
            double w1_invcm,w3_invcm;
            int i1=w1,i3=w3;
            if (w1<NTP1/2) i1+=NTP1/2; // w1=0..NTP1/2-1 gives the negative freq
            if (w3<NTP3/2) i3+=NTP3/2;
            if (w1>=NTP1/2) i1-=NTP1/2; // w1=NTP1/2..NTP1 gives the positive freq
            if (w3>=NTP3/2) i3-=NTP3/2;
            w1_invcm=double(w1-NTP1/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP1-1)*DNT*dt);
            w3_invcm=double(w3-NTP3/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP3-1)*DNT*dt);
            if(w1_invcm>Emin and w1_invcm<Emax and w3_invcm>Emin and w3_invcm<Emax)
            {
               // FIXME ugly hack to swap w1 <-> w3
               enter=true;
               data_sum_f[i1*NTP3+i3][0]+=data_GBNR_f[i1*NTP3+i3][0];
               data_sum_f[i1*NTP3+i3][1]+=data_GBNR_f[i1*NTP3+i3][1];
               fprintf(fo, "%e %e %e %e\n",w3_invcm,w1_invcm,data_GBNR_f[i1*NTP3+i3][0], data_GBNR_f[i1*NTP3+i3][1]);
            }
         }
         if(enter==true)
         {
            fprintf(fo, "\n");
         }
      }
      fclose(fo);
      fftw_free(data_GBNR_t);
      fftw_free(data_GBNR_f);
   }
   //
   if(ESANR==1)
   {
      /// TIME line output example SERP
      sprintf(fn,"%s_2d-echo_delay%d_ESANR.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway ESA-NR, time signal\n#\n");
      fprintf(fo,"#  t1 in s     t2 in s     Re[signal]    Im[signal]\n");
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            fprintf(fo, "%e %e %e %e\n",NstepsOverDNT*DNT*dt, t1*dt*DNT, data_ESANR_t[t1*NTP3+NstepsOverDNT][0], data_ESANR_t[t1*NTP3+NstepsOverDNT][1]);
         }
         fprintf(fo, "\n");
      }
      fclose(fo);
      // perform FFT
      for(int t1=0; t1<NTP1; t1++)
      {
         data_ESANR_t[t1*NTP3+0][0]=0.5*(data_ESANR_t[t1*NTP3+0][0]+data_ESANR_t[t1*NTP3+(NTP3-1)][0]);
         data_ESANR_t[t1*NTP3+0][1]=0.5*(data_ESANR_t[t1*NTP3+0][1]+data_ESANR_t[t1*NTP3+(NTP3-1)][1]);
      }
      for(int t3=0; t3<NTP3; t3++)
      {
         data_ESANR_t[0*NTP3+t3][0]=0.5*(data_ESANR_t[0*NTP3+t3][0]+data_ESANR_t[(NTP1-1)*NTP3+t3][0]);
         data_ESANR_t[0*NTP3+t3][1]=0.5*(data_ESANR_t[0*NTP3+t3][1]+data_ESANR_t[(NTP1-1)*NTP3+t3][1]);
      }
      perform_fftNR(data_ESANR_t, data_ESANR_f);
      sprintf(fn,"%s_2d-echo_delay%d_ESANR_freq.dat",fn_output,delayTime);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# %s\n",info_shotdipole);
      fprintf(fo,"# output: 2d-echo pathway ESA-NR, spectrum\n#\n");
      fprintf(fo,"# w1 in cm^-1  w2 in cm^-1  Re[signal]  Im[signal]\n");
      for(int w1=0; w1<NTP1; w1++)
      {
         bool enter=false;
         for(int w3=0; w3<NTP3; w3++)
         {
            double w1_invcm,w3_invcm;
            int i1=w1,i3=w3;
            if (w1<NTP1/2) i1+=NTP1/2; // w1=0..NTP1/2-1 gives the negative freq
            if (w3<NTP3/2) i3+=NTP3/2;
            if (w1>=NTP1/2) i1-=NTP1/2; // w1=NTP1/2..NTP1 gives the positive freq
            if (w3>=NTP3/2) i3-=NTP3/2;
            w1_invcm=double(w1-NTP1/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP1-1)*DNT*dt);
            w3_invcm=double(w3-NTP3/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP3-1)*DNT*dt);
            if(w1_invcm>Emin and w1_invcm<Emax and w3_invcm>Emin and w3_invcm<Emax)
            {
               // FIXME ugly hack to swap w1 <-> w3
               enter=true;
               data_sum_f[i1*NTP3+i3][0]+=data_ESANR_f[i1*NTP3+i3][0];
               data_sum_f[i1*NTP3+i3][1]+=data_ESANR_f[i1*NTP3+i3][1];
               fprintf(fo, "%e %e %e %e\n",w3_invcm,w1_invcm,data_ESANR_f[i1*NTP3+i3][0], data_ESANR_f[i1*NTP3+i3][1]);
            }
         }
         if(enter==true)
         {
            fprintf(fo, "\n");
         }
      }
      fclose(fo);
      fftw_free(data_ESANR_t);
      fftw_free(data_ESANR_f);
   }
// TOTAL Spectrum
   sprintf(fn,"%s_2d-echo_delay%d_sum_freq.dat",fn_output,delayTime);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time(&t);
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# %s\n",info_shotdipole);
   fprintf(fo,"%s spectrum\n#\n",Fileheader_calculatedPath);
   fprintf(fo,"# w1 in cm^-1  w2 in cm^-1  Re[signal]  Im[signal]\n");
   for(int w1=0; w1<NTP1; w1++)
   {
      bool enter=false;
      for(int w3=0; w3<NTP3; w3++)
      {
         double w1_invcm,w3_invcm;
         int i1=w1,i3=w3;
         if (w1<NTP1/2) i1+=NTP1/2; // w1=0..NTP1/2-1 gives the negative freq
         if (w3<NTP3/2) i3+=NTP3/2;
         if (w1>=NTP1/2) i1-=NTP1/2; // w1=NTP1/2..NTP1 gives the positive freq
         if (w3>=NTP3/2) i3-=NTP3/2;
         w1_invcm=double(w1-NTP1/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP1-1)*DNT*dt);
         w3_invcm=double(w3-NTP3/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP3-1)*DNT*dt);
         if(w1_invcm>Emin and w1_invcm<Emax and w3_invcm>Emin and w3_invcm<Emax)
         {
            // FIXME ugly hack to swap w1 <-> w3
            enter=true;
            fprintf(fo, "%e %e %e %e\n",w3_invcm,w1_invcm,data_sum_f[i1*NTP3+i3][0], data_sum_f[i1*NTP3+i3][1]);
         }
      }
      if(enter==true)
      {
         fprintf(fo, "\n");
      }
   }
   fclose(fo);
   fftw_free(data_sum_f);
}


void calculate_2d::execute_calculate_2d()
{
   system->create_mum_mup();
   fprintf(stdout, "# initialize propagation method=HEOM ...\n");
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

   int Nmod=2*Nsteps3+delayTime+1; // no output

   // SE-RP path
   if(SERP)
   {
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.0);
         MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
         evolution.initialize(*rhoInit);
         evolution.Add_Liouville(liouHexci);
         evolution.Add_Liouville(liouPhon);
         // define the second and third dipole multiplications:
         muminus=false;	// since I need mu_plus
         left=true;	// since I need to multiply from the left
         Multiply_Dipole_Sigmas RP_SE2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=false;
         Multiply_Dipole_Sigmas RP_SE3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(RP_SE2);
         evolution.Add_ManipulatingObservable(RP_SE3);
         // returns the trace of mu_minus * rho
         ReturnTrace_mumRho RP_SE_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(RP_SE_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_SERP_t[t1*NTP3+NstepsOverDNT][0]+=RP_SE_out.tracer[t1];
            data_SERP_t[t1*NTP3+NstepsOverDNT][1]+=RP_SE_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         fprintf(stdout,"# PROGRESS=>%d percent of SERP method=HEOM done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         fflush(stdout);
      }
      fprintf(stdout,"# SERP method=HEOM done\n");
   }
   if(SENR)
   {
      // SE-NR path
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
         MatrixMuliplication(rhoInit, system->mup, rho0, system->Nsites);
         evolution.initialize(*rhoInit);
         evolution.Add_Liouville(liouHexci);
         evolution.Add_Liouville(liouPhon);
         // define the second and third dipole multiplications:
         muminus=true;	// since I need muplus
         left=false;	// since I need to multiply from the left
         Multiply_Dipole_Sigmas NR_SE2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=false;
         Multiply_Dipole_Sigmas NR_SE3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(NR_SE2);
         evolution.Add_ManipulatingObservable(NR_SE3);
         // returns the trace of mu_minus * rho
         ReturnTrace_mumRho NR_SE_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(NR_SE_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_SENR_t[t1*NTP3+NstepsOverDNT][0]+=NR_SE_out.tracer[t1];
            data_SENR_t[t1*NTP3+NstepsOverDNT][1]+=NR_SE_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         fprintf(stdout,"# PROGRESS=>%d percent of SENR method=HEOM done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         fflush(stdout);
      }
      fprintf(stdout,"# SENR method=HEOM done\n");
   }
   if (GBRP)
   {
      // GB-RP path
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
         MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
         evolution.initialize(*rhoInit);
         evolution.Add_Liouville(liouHexci);
         evolution.Add_Liouville(liouPhon);
         muminus=false;	// since I need muplus
         left=false;	// since I need to multiply from the right
         Multiply_Dipole_Sigmas RP_GB2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas RP_GB3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(RP_GB2);
         evolution.Add_ManipulatingObservable(RP_GB3);
         ReturnTrace_mumRho RP_GB_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(RP_GB_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_GBRP_t[t1*NTP3+NstepsOverDNT][0]+=RP_GB_out.tracer[t1];
            data_GBRP_t[t1*NTP3+NstepsOverDNT][1]+=RP_GB_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         fprintf(stdout,"# PROGRESS=>%d percent of GBRP method=HEOM done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         fflush(stdout);
      }
      fprintf(stdout,"# GBRP method=HEOM done\n");
   }
   if (GBNR)
   {
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
         MatrixMuliplication(rhoInit, system->mup, rho0, system->Nsites);
         evolution.initialize(*rhoInit);
         evolution.Add_Liouville(liouHexci);
         evolution.Add_Liouville(liouPhon);
         muminus=true;	// since I need muplus
         left=true;	// since I need to multiply from the right
         Multiply_Dipole_Sigmas NR_GB2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas NR_GB3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(NR_GB2);
         evolution.Add_ManipulatingObservable(NR_GB3);
         ReturnTrace_mumRho NR_GB_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(NR_GB_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_GBNR_t[t1*NTP3+NstepsOverDNT][0]+=NR_GB_out.tracer[t1];
            data_GBNR_t[t1*NTP3+NstepsOverDNT][1]+=NR_GB_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         fprintf(stdout,"# PROGRESS=>%d percent of GBNR method=HEOM done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         fflush(stdout);
      }
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

   if (ESARP)
   {
      // ESA-RP path
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
         MatrixMuliplication(rhoInit, rho0, system->mum, system->Nsites);
         evolution.initialize(*rhoInit);
         evolution.Add_Liouville(liouHexciDouble);
         evolution.Add_Liouville(liouPhonDouble);
         muminus=false;	// since I need muplus
         left=true;	// since I need to multiply from the left
         Multiply_Dipole_Sigmas RP_ESA2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas RP_ESA3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(RP_ESA2);
         evolution.Add_ManipulatingObservable(RP_ESA3);
         ReturnTrace_mumRho RP_ESA_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(RP_ESA_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         // ESA goes with minus
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_ESARP_t[t1*NTP3+NstepsOverDNT][0]-=RP_ESA_out.tracer[t1];
            data_ESARP_t[t1*NTP3+NstepsOverDNT][1]-=RP_ESA_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         fprintf(stdout,"# PROGRESS=>%d percent of ESARP method=HEOM done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         fflush(stdout);
      }
      fprintf(stdout,"# ESARP method=HEOM done\n");
   }
   if (ESANR)
   {
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);
         rho0->assignEntry(0,0,1.);
         MatrixMuliplication(rhoInit,system->mup, rho0,  system->Nsites);
         evolution.initialize(*rhoInit);
         evolution.Add_Liouville(liouHexciDouble);
         evolution.Add_Liouville(liouPhonDouble);
         muminus=true;	// since I need muplus
         left=false;	// since I need to multiply from the left
         Multiply_Dipole_Sigmas NR_ESA2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas NR_ESA3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(NR_ESA2);
         evolution.Add_ManipulatingObservable(NR_ESA3);
         ReturnTrace_mumRho NR_ESA_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(NR_ESA_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         // ESA goes with minus
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_ESANR_t[t1*NTP3+NstepsOverDNT][0]-=NR_ESA_out.tracer[t1];
            data_ESANR_t[t1*NTP3+NstepsOverDNT][1]-=NR_ESA_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         fprintf(stdout,"# PROGRESS=>%d percent of ESANR method=HEOM done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         fflush(stdout);
      }
      fprintf(stdout,"# ESANR method=HEOM done\n");
   }
   liouPhonDouble.freeDeviceMemory();
}


void calculate_2d::execute_calculate_2d_rotav_10shots(const char* method)
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
         execute_calculate_2d();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_2d(method);
      }
   }
   write_data("ten shot rotational average");
}


void calculate_2d::execute_calculate_2d_rotav_4shots(const char* method)
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
         execute_calculate_2d();
      }
      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
      {
         execute_calculate_Redfield_2d(method);
      }
   }
   write_data("four shot rotational average");
}





void calculate_2d::execute_calculate_Redfield_2d(const char* method)
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
   int Nmod=2*Nsteps3+delayTime+1; // no output
   // SE-RP path
   if(SERP)
   {
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
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
         Multiply_Dipole_Sigmas_Redfield RP_SE2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=false;
         Multiply_Dipole_Sigmas_Redfield RP_SE3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(RP_SE2);
         evolution.Add_ManipulatingObservable(RP_SE3);
         // returns the trace of mu_minus * rho
         ReturnTrace_mumRho_Redfield RP_SE_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(RP_SE_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_SERP_t[t1*NTP3+NstepsOverDNT][0]+=RP_SE_out.tracer[t1];
            data_SERP_t[t1*NTP3+NstepsOverDNT][1]+=RP_SE_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of SERP method=full Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of SERP method=secular Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         fflush(stdout);
      }
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
   if(SENR)
   {
      // SE-NR path
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
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
         Multiply_Dipole_Sigmas_Redfield NR_SE2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=false;
         Multiply_Dipole_Sigmas_Redfield NR_SE3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(NR_SE2);
         evolution.Add_ManipulatingObservable(NR_SE3);
         // returns the trace of mu_minus * rho
         ReturnTrace_mumRho_Redfield NR_SE_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(NR_SE_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_SENR_t[t1*NTP3+NstepsOverDNT][0]+=NR_SE_out.tracer[t1];
            data_SENR_t[t1*NTP3+NstepsOverDNT][1]+=NR_SE_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of SENR method=full Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of SENR method=secular Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         fflush(stdout);
      }
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
   if (GBRP)
   {
      // GB-RP path
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
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
         Multiply_Dipole_Sigmas_Redfield RP_GB2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas_Redfield RP_GB3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(RP_GB2);
         evolution.Add_ManipulatingObservable(RP_GB3);
         ReturnTrace_mumRho_Redfield RP_GB_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(RP_GB_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_GBRP_t[t1*NTP3+NstepsOverDNT][0]+=RP_GB_out.tracer[t1];
            data_GBRP_t[t1*NTP3+NstepsOverDNT][1]+=RP_GB_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of GBRP method=full Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of GBRP method=secular Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         fflush(stdout);
      }
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
   if (GBNR)
   {
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
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
         Multiply_Dipole_Sigmas_Redfield NR_GB2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas_Redfield NR_GB3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(NR_GB2);
         evolution.Add_ManipulatingObservable(NR_GB3);
         ReturnTrace_mumRho_Redfield NR_GB_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(NR_GB_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_GBNR_t[t1*NTP3+NstepsOverDNT][0]+=NR_GB_out.tracer[t1];
            data_GBNR_t[t1*NTP3+NstepsOverDNT][1]+=NR_GB_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of GBNR method=full Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of GBNR method=secular Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
         fflush(stdout);
      }
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

   if (ESARP)
   {
      // ESA-RP path
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);  // Nsteps-1 generates Nsteps output values (from 0 to Nsteps-1)
         rho0->assignEntry(0,0,1.);
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
         Multiply_Dipole_Sigmas_Redfield RP_ESA2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas_Redfield RP_ESA3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(RP_ESA2);
         evolution.Add_ManipulatingObservable(RP_ESA3);
         ReturnTrace_mumRho_Redfield RP_ESA_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(RP_ESA_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         // ESA goes with minus
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_ESARP_t[t1*NTP3+NstepsOverDNT][0]-=RP_ESA_out.tracer[t1];
            data_ESARP_t[t1*NTP3+NstepsOverDNT][1]-=RP_ESA_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of ESARP method=full Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
//         if(strncmp(method,"secular Redfield",strlen(method))==0)
//         {
//            fprintf(stdout,"# PROGRESS=>%d percent of ESARP method=secular Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
//         }
         fflush(stdout);
      }
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout,"# ESARP method=full Redfield done\n");
      }
//            if(strncmp(method,"secular Redfield",strlen(method))==0)
//            {
//               fprintf(stdout,"# ESARP method=secular Redfield done\n");
//            }
//            fflush(stdout);
   }
   if (ESANR)
   {
      for (int NstepsOverDNT=0; NstepsOverDNT<Nsteps3/DNT; NstepsOverDNT++) // Nsteps=NstepsOverDNT * DNT *dt gives the time t1
      {
         Propagation evolution(*system, *OpenCLinfo, *sigma_Host, NstepsOverDNT*DNT+Nsteps3+delayTime, dt, Nmax);
         rho0->assignEntry(0,0,1.);
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
         Multiply_Dipole_Sigmas_Redfield NR_ESA2(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT,Ntuples);
         muminus=false;
         left=true;
         Multiply_Dipole_Sigmas_Redfield NR_ESA3(*system, *OpenCLinfo, muminus,left,NstepsOverDNT*DNT+delayTime,Ntuples);
         evolution.Add_ManipulatingObservable(NR_ESA2);
         evolution.Add_ManipulatingObservable(NR_ESA3);
         ReturnTrace_mumRho_Redfield NR_ESA_out(*system, *OpenCLinfo, DNT, system->mum,Help3,NstepsOverDNT*DNT+delayTime);
         evolution.Add_ManipulatingObservable(NR_ESA_out);
         evolution.propagate(Nmod); // Nmod defines the time steps for an output
         // ESA goes with minus
         for(int t1=0; t1<Nsteps3/DNT; t1++)
         {
            data_ESANR_t[t1*NTP3+NstepsOverDNT][0]-=NR_ESA_out.tracer[t1];
            data_ESANR_t[t1*NTP3+NstepsOverDNT][1]-=NR_ESA_out.tracei[t1];
         }
         evolution.freeDeviceMemory();
         if(strncmp(method,"full Redfield",strlen(method))==0)
         {
            fprintf(stdout,"# PROGRESS=>%d percent of ESANR method=full Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
         }
//         if(strncmp(method,"secular Redfield",strlen(method))==0)
//         {
//            fprintf(stdout,"# PROGRESS=>%d percent of ESANR method=secular Redfield done\n",100*(NstepsOverDNT+1)/(Nsteps3/DNT));
//         }
         fflush(stdout);
      }
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



#endif /* GPU_RUN2DSPECTRA_H_ */
