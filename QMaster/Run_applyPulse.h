/*
 * Run_applyPulse.h
 *
 *  Created on: Feb 17, 2014
 *      Author: christoph
 */

#ifndef RUN_APPLYPULSE_H_
#define RUN_APPLYPULSE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "headers.h"


// Many functions and varibables are not used right now. FIXME clean-up
class calculate_applyPulse
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Nmax;
   double dt;
   int DNT;
   int Nsteps;
   double** sd_lambda_invcm;
   double** sd_gamma_fs;
   double** sd_gammapos_invcm;
   vector<int> sd_CoupledSite;
   vector<int> sd_Npeaks;
   double sigma_pulse;
   double center_pulse;
   double strength_pulse;
   double omega_pulse;
   int padding;
   bool rotav;
   double* Pol;
   const char *fn_output;
   const char *dipole_moments;
   double Emin;
   double Emax;
   bool specify_Enrange;
   bool return_densityMatrix; // if specified to return all entries of the desity matrix in site-basis
   double diff;
private:
   struct timeval t0, t1;
   bool use_specified_laserDirection;
   double* lambda;
   double_complex* gamma;
   double Ediffmax_invcm;
   int NTP;
   SigmaHelpMatrices *sigma_Host;
   int Ntuples;
   double_complex* data_t;
public:
   calculate_applyPulse(System& INsystem,
                        OpenCL_Init& INOpenCLinfo,
                        int INNmax,
                        double INdt,
                        int INDNT,
                        int INNsteps,
                        double** INsd_lambda_invcm,
                        double** INsd_gamma_fs,
                        double** INsd_gammapos_invcm,
                        vector<int> INsd_CoupledSite,
                        vector<int> INsd_Npeaks,
                        double INsigma_pulse,
                        double INcenter_pulse,
                        double INstrength_pulse,
                        double INomega_pulse,
                        int INpadding,
                        bool INrotav,
                        double* INPol,
                        const char* INfn_output,
                        const char* INdipole_moments,
                        double INEmin,
                        double INEmax,
                        bool INspecify_Enrange,
                        bool INreturn_densityMatrix // if specified to return all entries of the desity matrix in site-basis
                       );
   void execute_calculate_applyPulse();
   void execute_calculate_Redfield_applyPulse(const char* method);
};


calculate_applyPulse::calculate_applyPulse(System& INsystem,
      OpenCL_Init& INOpenCLinfo,
      int INNmax,
      double INdt,
      int INDNT,
      int INNsteps,
      double** INsd_lambda_invcm,
      double** INsd_gamma_fs,
      double** INsd_gammapos_invcm,
      vector<int> INsd_CoupledSite,
      vector<int> INsd_Npeaks,
      double INsigma_pulse,
      double INcenter_pulse,
      double INstrength_pulse,
      double INomega_pulse,
      int INpadding,
      bool INrotav,
      double* INPol,
      const char* INfn_output,
      const char* INdipole_moments,
      double INEmin,
      double INEmax,
      bool INspecify_Enrange,
      bool INreturn_densityMatrix
                                          )  :
   system(&INsystem),
   OpenCLinfo(&INOpenCLinfo),
   Nmax(INNmax),
   dt(INdt),
   DNT(INDNT),
   Nsteps(INNsteps),
   sd_lambda_invcm(INsd_lambda_invcm),
   sd_gamma_fs(INsd_gamma_fs),
   sd_gammapos_invcm(INsd_gammapos_invcm),
   sd_CoupledSite(INsd_CoupledSite),
   sd_Npeaks(INsd_Npeaks),
   sigma_pulse(INsigma_pulse),
   center_pulse(INcenter_pulse),
   strength_pulse(INstrength_pulse),
   omega_pulse(INomega_pulse),
   padding(INpadding),
   rotav(INrotav),
   Pol(INPol),
   fn_output(INfn_output),
   dipole_moments(INdipole_moments),
   Emin(INEmin),
   Emax(INEmax),
   specify_Enrange(INspecify_Enrange),
   return_densityMatrix(INreturn_densityMatrix)
{
   use_specified_laserDirection=false;
   fprintf(stdout, "#\n# initialize parameter for absorption spectrum\n");
   // specify maximal difference between to exction eigen energies to specify plot range of the spectral density (plot energy range where differences in exciton
   // energies live (extend by factor of 1.5 to be sure that important part is plotted))
   system->computeEigenvects_Eigenvals();
   if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
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
   }
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
   diff=0;
   NTP=Nsteps*padding; // steps along time including padding
   //max energy range
   double Eminall=-NTP*0.5*const_hbar/const_invcmtomeV*2.*M_PI/((NTP-1)*dt);
   double Emaxall=+NTP*0.5*const_hbar/const_invcmtomeV*2.*M_PI/((NTP-1)*dt);
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

void calculate_applyPulse::execute_calculate_applyPulse()
{
   fftw_plan fftplan1d;
   int sizeN=NTP*sizeof(fftw_complex);
   char fn[500];
   int num_shots;
   FILE *fo;

   fftw_complex *data_f,*data_fsum; //,*data_t
   data_t=new double_complex[Nsteps*system->Nsites];
   for(int i=0; i<Nsteps*system->Nsites; i++)
   {
      data_t[i]=0;
   }

   // arrays to hold time line and its FFT to energy space
   // data_t=(fftw_complex*)fftw_malloc(sizeN);
   data_f=(fftw_complex*)fftw_malloc(sizeN);
   data_fsum=(fftw_complex*)fftw_malloc(sizeN);
   for(int i=0; i<NTP; i++)
   {
      //  data_t[i][0]=0.0;
      //  data_t[i][1]=0.0;
      data_fsum[i][0]=0.0;
      data_fsum[i][1]=0.0;
   }

   fftplan1d=fftw_plan_dft_1d(NTP,(fftw_complex*)data_t,(fftw_complex*)data_f,FFTW_BACKWARD,FFTW_ESTIMATE);

   if(rotav==true)
   {
      num_shots=3;
   }
   else
   {
      num_shots=1;
      system->set_new_laser_direction(Pol[0],Pol[1],Pol[2]); // default Pol=(1,0,0) if not default Pol is specified by keyword laser_polarization
   }
   fprintf(stdout,"# specified dipole moments: %s\n",dipole_moments);

   int Nmod=Nsteps/10;  // every 10% progress report to stdout
   fprintf(stdout, "# initialize propagation method=HEOM ...\n");
   sigma_Host=new SigmaHelpMatrices(system->Nsites, Nmax, system->Ncoupled);
   Ntuples=sigma_Host->AllsigmaTuples.size();
   Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
   fprintf(stdout, "# ... done\n");
   LiouvilleHExciton liouHExc(*system, *OpenCLinfo, Ntuples);
   LiouvilleHPulse liouHPulse(*system, *OpenCLinfo, Ntuples, sigma_pulse, center_pulse, strength_pulse, omega_pulse);
   LiouvillePhononLowTemp liouPhon(*system, *OpenCLinfo,  *sigma_Host, Nmax);
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
   evolution.Add_Liouville(liouHExc);
   evolution.Add_LiouvilleTime(liouHPulse);
   evolution.Add_Liouville(liouPhon);
   returnDensityMatrixSiteBasis *totalrho;
   GetDensityMatrixSiteBasisDiag *printrho;
   for(int shot=0; shot<num_shots; shot++)
   {
      if (rotav==false and (strcmp(dipole_moments,"effective dipole moments eigen")!=0) and (strcmp(dipole_moments,"effective dipole moments site")!=0))
      {
         fprintf(stdout, "# laser polarization along l=(%e,%e,%e)\n",Pol[0],Pol[1],Pol[2]);
         use_specified_laserDirection=true;
      }
      if ((rotav==true)&&(shot==0))
      {
         system->set_new_laser_direction(1.0,0.0,0.0);
      }
      if ((rotav==true)&&(shot==1))
      {
         system->set_new_laser_direction(0.0,1.0,0.0);
      }
      if ((rotav==true)&&(shot==2))
      {
         system->set_new_laser_direction(0.0,0.0,1.0);
      }
      if (rotav==true)
      {
         fprintf(stdout, "#\n# perform rotational average shot=%i\n",shot+1);
      }
      system->create_mum_mup();
      Matrix rhoInit(system->Nsites);
      rhoInit.assignEntry(0,0,1.0);
      evolution.initialize(rhoInit);
      printrho=new GetDensityMatrixSiteBasisDiag(*system,*OpenCLinfo,DNT, dt);
      evolution.Add_Observable(*printrho);
      if(return_densityMatrix==true)
      {
         totalrho=new returnDensityMatrixSiteBasis(*system,*OpenCLinfo,DNT, dt);
         evolution.Add_Observable(*totalrho);
      }
      fprintf(stdout, "#\n# starting propagation method=HEOM ...\n");
      gettimeofday(&t0, 0);
      evolution.propagate(Nmod);
      cout<<"# ... done"<<endl;
      gettimeofday(&t1, 0);
      double diff_local = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
      diff+=diff_local;
      evolution.clear_Sigma();
   }
   cout<<"# free Device-Memory ..."<<endl;
   liouPhon.freeDeviceMemory();
   evolution.freeDeviceMemory();
   liouHExc.freeDeviceMemory();
   liouHPulse.freeDeviceMemory();
   clFlush(OpenCLinfo->queue);
   clFinish(OpenCLinfo->queue);
   cout<<"# ... done"<<endl;


   /*
          *********************************************************
           * Save output values to file
           *********************************************************
          */
   sprintf(fn,"%s_applyPulse.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time_t t;
   time(&t);
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: population dynamics\n#\n");
   fprintf(fo,"# t in s       ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"site %i       ",j);
   }
   fprintf(fo,"\n");
   for(int it=0; it<printrho->Ndata; it++)
   {
      fprintf(fo,"%e ",(double)it*dt*DNT);
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo,"%e ",abs(printrho->data[it*system->Nsites+j]));
      }
      fprintf(fo,"\n");
   }
   fclose(fo);

   if(return_densityMatrix==true)
   {
      for(int i=0; i<system->Nsites; i++)
      {
         sprintf(fn,"%s_applyPulse_densityMatrix_rho_%i_k.dat",fn_output, i);
         FILE* out;
         out=fopen(fn,"w");
         if (out==NULL)
         {
            fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
            exit(1);
         }
         time_t t;
         time(&t);
         fprintf(out,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
         fprintf(out,"# output: applyPulse, all entries of the density Matrix of column %i \n#\n",i);
         fprintf(out,"# t in s  ");
         for(int j=0; j<system->Nsites; j++)
         {
            fprintf(out,"Re rho(%i,%i) Im rho(%i,%i) ",i,j, i ,j);
         }

         fprintf(out,"\n");
         for(int it=0; it<totalrho->Ndata; it++)
         {
            fprintf(out,"%e ",(double)it*dt*DNT);
            for(int j=0; j<system->Nsites; j++)
            {
               int id=i*system->Nsites+j;
               fprintf(out,"%e %e ",real(totalrho->data[it*system->Nsites*system->Nsites+id]),imag(totalrho->data[it*system->Nsites*system->Nsites+id]));
            }
            fprintf(out,"\n");
         }
         fclose(out);
      }
   }
   /*
     *********************************************************
      * Save output values to file
      *********************************************************
     */
//   sprintf(fn,"%s_spectrum1d.dat",fn_output);
//   fo=fopen(fn,"w");
//   if (fo==NULL)
//   {
//      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
//      exit(1);
//   }
//   time_t t;
//   time(&t);
//   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
//   if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(rotav==true)
//   {
//      fprintf(fo,"# xyz-rotational average\n");
//   }
//   if(use_specified_laserDirection==true)
//   {
//      fprintf(fo,"# laser polarization along l=(%e,%e,%e)\n",Pol[0],Pol[1],Pol[2]);
//   }
//   fprintf(fo,"# output: absorption spectrum, time signal\n#\n");
//   fprintf(fo,"#  t in s     Re[signal]   Im[signal] \n");
//   for(int i=0; i<NTP; i++)
//   {
//      fprintf(fo, "%e %e %e\n",i*dt, data_t[i][0],data_t[i][1]);
//   }
//   fclose(fo);
//   // perform FFT
//   data_t[0][0]=0.5*(data_t[0][0]+data_t[NTP-1][0]);
//   data_t[0][1]=0.5*(data_t[0][1]+data_t[NTP-1][1]);
//   fftw_execute(fftplan1d);
////   for(int i=0;i<NTP;i++) { data_fsum[i][0]+=data_f[i][0]; data_fsum[i][1]+=data_f[i][1]; }
//   sprintf(fn,"%s_spectrum1d_freq.dat",fn_output);
//   fo=fopen(fn,"w");
//   if (fo==NULL)
//   {
//      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
//      exit(1);
//   }
//   time(&t);
//   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
//   if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(rotav==true)
//   {
//      fprintf(fo,"# xyz-rotational average\n");
//   }
//   if(use_specified_laserDirection==true)
//   {
//      fprintf(fo,"# laser polarization along l=(%e,%e,%e)\n",Pol[0],Pol[1],Pol[2]);
//   }
//   fprintf(fo,"# output: absorption spectrum\n#\n");
//   fprintf(fo,"#  w in cm^-1   Re[signal]  Im[signal] \n");
//   for(int w=0; w<NTP; w++)
//   {
//      double w_invcm;
//      int i=w;
//      if (w<NTP/2)  i+=NTP/2; // w=0..NTP/2-1 gives the negative freq
//      if (w>=NTP/2) i-=NTP/2; // w=NTP/2..NTP gives the positive freq
//      w_invcm=double(w-NTP/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP-1)*dt);
//      if(w_invcm>Emin and w_invcm<Emax)
//      {
//         fprintf(fo, "%e %e %e\n",w_invcm,data_f[i][0], data_f[i][1]);
//      }
//   }
//   fclose(fo);
}

void calculate_applyPulse::execute_calculate_Redfield_applyPulse(const char* method)
{
   fftw_plan fftplan1d;
   int sizeN=NTP*sizeof(fftw_complex);
   char fn[500];
   int num_shots;
   FILE *fo;

   fftw_complex *data_f,*data_fsum; //,*data_t
   data_t=new double_complex[Nsteps*system->Nsites];
   for(int i=0; i<Nsteps*system->Nsites; i++)
   {
      data_t[i]=0;
   }

   // arrays to hold time line and its FFT to energy space
   // data_t=(fftw_complex*)fftw_malloc(sizeN);
   data_f=(fftw_complex*)fftw_malloc(sizeN);
   data_fsum=(fftw_complex*)fftw_malloc(sizeN);
   for(int i=0; i<NTP; i++)
   {
      //  data_t[i][0]=0.0;
      //  data_t[i][1]=0.0;
      data_fsum[i][0]=0.0;
      data_fsum[i][1]=0.0;
   }

   fftplan1d=fftw_plan_dft_1d(NTP,(fftw_complex*)data_t,(fftw_complex*)data_f,FFTW_BACKWARD,FFTW_ESTIMATE);

   if(rotav==true)
   {
      num_shots=3;
   }
   else
   {
      num_shots=1;
      system->set_new_laser_direction(Pol[0],Pol[1],Pol[2]); // default Pol=(1,0,0) if not default Pol is specified by keyword laser_polarization
   }
   fprintf(stdout,"# specified dipole moments: %s\n",dipole_moments);
   int Nmod=Nsteps/10;  // every 10% progress report to stdout
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
   Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
   fprintf(stdout, "# ... done\n");
   LiouvilleHExcitonRedfield liouHExc(*system,  *OpenCLinfo);
   LiouvilleHPulseRedfield liouHPulse(*system, *OpenCLinfo, Ntuples, sigma_pulse, center_pulse, strength_pulse, omega_pulse);
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
   evolution.Add_Liouville(liouHExc);
   evolution.Add_LiouvilleTime(liouHPulse);
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      evolution.Add_Liouville(*liouPhonfull);
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      evolution.Add_Liouville(*liouPhonsecular);
   }
   returnDensityMatrixSiteBasis_fullRedfield *totalrho;
   GetDensityMatrixSiteBasisDiag_fullRedfield *printrho;
   for(int shot=0; shot<num_shots; shot++)
   {
      if (rotav==false and (strcmp(dipole_moments,"effective dipole moments eigen")!=0) and (strcmp(dipole_moments,"effective dipole moments site")!=0))
      {
         fprintf(stdout, "# laser polarization along l=(%e,%e,%e)\n",Pol[0],Pol[1],Pol[2]);
         use_specified_laserDirection=true;
      }
      if ((rotav==true)&&(shot==0))
      {
         system->set_new_laser_direction(1.0,0.0,0.0);
      }
      if ((rotav==true)&&(shot==1))
      {
         system->set_new_laser_direction(0.0,1.0,0.0);
      }
      if ((rotav==true)&&(shot==2))
      {
         system->set_new_laser_direction(0.0,0.0,1.0);
      }
      if (rotav==true)
      {
         fprintf(stdout, "#\n# perform rotational average shot=%i\n",shot+1);
      }
      system->create_mum_mup();
      Matrix rhoInit(system->Nsites);
      // here no need to rotate in eigenbasis, since we start in ground-state
      rhoInit.assignEntry(0,0,1.0);
      evolution.initialize(rhoInit);
      printrho=new GetDensityMatrixSiteBasisDiag_fullRedfield(*system,*OpenCLinfo,DNT, dt);
      evolution.Add_Observable(*printrho);
      if(return_densityMatrix==true)
      {
         totalrho=new returnDensityMatrixSiteBasis_fullRedfield(*system,*OpenCLinfo,DNT, dt);
         evolution.Add_Observable(*totalrho);
      }
      if(strncmp(method,"full Redfield",strlen(method))==0)
      {
         fprintf(stdout, "#\n# starting propagation method=full Redfield ...\n");
      }
      if(strncmp(method,"secular Redfield",strlen(method))==0)
      {
         fprintf(stdout, "#\n# starting propagation method=secular Redfield ...\n");
      }
      gettimeofday(&t0, 0);
      evolution.propagate(Nmod);
      cout<<"# ... done"<<endl;
      gettimeofday(&t1, 0);
      double diff_local = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
      diff+=diff_local;
      evolution.clear_Sigma();
   }
   cout<<"# free Device-Memory ..."<<endl;
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
      liouPhonfull->freeDeviceMemory();
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      liouPhonsecular->freeDeviceMemory();
   }
   evolution.freeDeviceMemory();
   liouHExc.freeDeviceMemory();
   liouHPulse.freeDeviceMemory();
   clFlush(OpenCLinfo->queue);
   clFinish(OpenCLinfo->queue);
   cout<<"# ... done"<<endl;


   /*
          *********************************************************
           * Save output values to file
           *********************************************************
          */
   sprintf(fn,"%s_applyPulse.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time_t t;
   time(&t);
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: population dynamics\n#\n");
   fprintf(fo,"# t in s       ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"site %i       ",j);
   }
   fprintf(fo,"\n");
   for(int it=0; it<printrho->Ndata; it++)
   {
      fprintf(fo,"%e ",(double)it*dt*DNT);
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo,"%e ",abs(printrho->data[it*system->Nsites+j]));
      }
      fprintf(fo,"\n");
   }
   fclose(fo);

   if(return_densityMatrix==true)
   {
      for(int i=0; i<system->Nsites; i++)
      {
         sprintf(fn,"%s_applyPulse_densityMatrix_rho_%i_k.dat",fn_output, i);
         FILE* out;
         out=fopen(fn,"w");
         if (out==NULL)
         {
            fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
            exit(1);
         }
         time_t t;
         time(&t);
         fprintf(out,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
         fprintf(out,"# output: applyPulse, all entries of the density Matrix of column %i \n#\n",i);
         fprintf(out,"# t in s  ");
         for(int j=0; j<system->Nsites; j++)
         {
            fprintf(out,"Re rho(%i,%i) Im rho(%i,%i) ",i,j, i ,j);
         }

         fprintf(out,"\n");
         for(int it=0; it<totalrho->Ndata; it++)
         {
            fprintf(out,"%e ",(double)it*dt*DNT);
            for(int j=0; j<system->Nsites; j++)
            {
               int id=i*system->Nsites+j;
               fprintf(out,"%e %e ",real(totalrho->data[it*system->Nsites*system->Nsites+id]),imag(totalrho->data[it*system->Nsites*system->Nsites+id]));
            }
            fprintf(out,"\n");
         }
         fclose(out);
      }
   }
   /*
     *********************************************************
      * Save output values to file
      *********************************************************
     */
//   sprintf(fn,"%s_spectrum1d.dat",fn_output);
//   fo=fopen(fn,"w");
//   if (fo==NULL)
//   {
//      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
//      exit(1);
//   }
//   time_t t;
//   time(&t);
//   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
//   if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(rotav==true)
//   {
//      fprintf(fo,"# xyz-rotational average\n");
//   }
//   if(use_specified_laserDirection==true)
//   {
//      fprintf(fo,"# laser polarization along l=(%e,%e,%e)\n",Pol[0],Pol[1],Pol[2]);
//   }
//   fprintf(fo,"# output: absorption spectrum, time signal\n#\n");
//   fprintf(fo,"#  t in s     Re[signal]   Im[signal] \n");
//   for(int i=0; i<NTP; i++)
//   {
//      fprintf(fo, "%e %e %e\n",i*dt, data_t[i][0],data_t[i][1]);
//   }
//   fclose(fo);
//   // perform FFT
//   data_t[0][0]=0.5*(data_t[0][0]+data_t[NTP-1][0]);
//   data_t[0][1]=0.5*(data_t[0][1]+data_t[NTP-1][1]);
//   fftw_execute(fftplan1d);
////   for(int i=0;i<NTP;i++) { data_fsum[i][0]+=data_f[i][0]; data_fsum[i][1]+=data_f[i][1]; }
//   sprintf(fn,"%s_spectrum1d_freq.dat",fn_output);
//   fo=fopen(fn,"w");
//   if (fo==NULL)
//   {
//      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
//      exit(1);
//   }
//   time(&t);
//   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
//   if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
//   {
//      fprintf(fo,"# %s\n", dipole_moments);
//   }
//   if(rotav==true)
//   {
//      fprintf(fo,"# xyz-rotational average\n");
//   }
//   if(use_specified_laserDirection==true)
//   {
//      fprintf(fo,"# laser polarization along l=(%e,%e,%e)\n",Pol[0],Pol[1],Pol[2]);
//   }
//   fprintf(fo,"# output: absorption spectrum\n#\n");
//   fprintf(fo,"#  w in cm^-1   Re[signal]  Im[signal] \n");
//   for(int w=0; w<NTP; w++)
//   {
//      double w_invcm;
//      int i=w;
//      if (w<NTP/2)  i+=NTP/2; // w=0..NTP/2-1 gives the negative freq
//      if (w>=NTP/2) i-=NTP/2; // w=NTP/2..NTP gives the positive freq
//      w_invcm=double(w-NTP/2.0)*const_hbar/const_invcmtomeV*2.0*M_PI/((NTP-1)*dt);
//      if(w_invcm>Emin and w_invcm<Emax)
//      {
//         fprintf(fo, "%e %e %e\n",w_invcm,data_f[i][0], data_f[i][1]);
//      }
//   }
//   fclose(fo);
}

#endif
