/* Copyright (C) 2010-2013 Christoph Kreisbeck
   Additional modifications by Tobias Kramer
   This is closed software.
   If you want to use it, please contact:
   christophkreisbeck@gmail.com
   to obtain permission.
*/
#define DIMBLOCKX 16
#define DIMBLOCKY 16
//#define PROFILING
// Duration for demo license
#define DEMOLICENCETIME 31 //guarantee 30 days, add one additional day


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "headers.h"
#include <string>
#include <time.h>

/*
   calling convention:
   argv[1]=name of input paramter file
   argv[2]=name of output file
*/

int main(int argc, char * argv[])
{
#ifdef DEVELOPMENT
   cout<<"# RUN QMaster IN DEVELOPMENT MODE"<<endl;
   cout<<"# "<<endl;
   // development mode affects in OpenCL_initilaize kernel sources:
   // -> in development mode kernels are from files Kernels/*.cl
   // -> default kernels are hardcoded and loaded from Kernels/KernelSource.h
#endif
#ifdef DEMO_MODE
#ifdef TIMESTARTLICENCE
   cout<<"# RUN QMaster IN DEMO MODE"<<endl;
   cout<<"# "<<endl;
   double lincence_end=TIMESTARTLICENCE+DEMOLICENCETIME*24*3600;
   struct timeval today;
   gettimeofday(&today, 0);
   //in seconds since 00:00:00 UTC (Coordinated Universal Time), January 1, 1970
   double licence_runtime=today.tv_sec-TIMESTARTLICENCE;
   if(today.tv_sec>lincence_end)
   {
      if((today.tv_sec-lincence_end)/3600<24)
      {
         fprintf(stdout,"# Licence Error: licence expired %.1f hours ago.\n",(today.tv_sec-lincence_end)/3600);
      }
      else
      {
         fprintf(stdout,"# Licence Error: licence expired %.1f days ago.\n",(today.tv_sec-lincence_end)/3600/24);
      }
      exit(1);
   }
//	   cout<<"licence in use: "<<licence_runtime/3600<<"h"<<endl;
#else
   fprintf(stdout,"# License Error: demo license has not been activated")
#endif
#endif
   double diff=0;
   struct timeval start_2d, end_2d;
   double compute_time=0;
   struct timeval t0, t1;
   gettimeofday(&t0, 0);
   fprintf(stdout,"# QMaster-%s\n#\n",INPUT_FORMAT_VERSION);
   bool efficiency_non_hermitian=true; //non-hermitian treatment of trapping and losses is defined as default
   bool specify_general_initial_state=false;
   bool specify_valrho_init=false;
   bool specify_version=false;
   bool specify_Enrange=false;
   bool specify_DNT=false;
   bool specify_Nsites=false;
   bool specify_Ham=false;
   bool specify_Temp=false;
   bool specify_Nmax=false;
   bool specify_Nsteps=false;
   bool specify_dt=false;
   bool specify_siteinit=false;
   bool specify_Ei=false;
   bool specify_Ej=false;
   bool specify_calculation_mode=false;
   bool specify_2dpathway=false;
   bool specify_method=false;
   bool specify_method_eigenstate_population=false;
   bool specify_Npeaks=false;
   bool specify_peaks=false;
   bool specify_sigma_pulse=false;
   bool specify_center_pulse=false;
   bool specify_strength_pulse=false;
   bool specify_omega_pulse=false;
   bool specify_Npeaks_all=false;
   bool specify_Nsites_coupled=false;
   bool specify_dipole_moments=false;
   bool specify_threshold_dominant_path=false;
   bool rotav=false;
   bool specify_vendor=false;
   bool specify_device_type=false;
   bool specify_shot=false;
   bool specify_laser_polarization=false;
   bool parse_error=false;
   bool specify_tmax=false;
   bool specify_NumIter=false;
   bool dominant_path=false;
   int readin_peaks=0;
   int readin_all=0;
   bool flux_analysis=false;
   bool return_densityMatrix=false;
   bool specify_Jcutoff=false;
   vector<int> sd_readin_peaks;
   vector<int> set_dipoles;
   vector<int> set_sd_peak_all;
   map<int,int> readin_sites;
   char version[500];
   int use_device_id=-1;
   int padding=4;
   int Nsites=0;  /* for input.group(excitonsystem).integer(Nsites) */
   int Ncoupled=0; /* number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site*/
   double Temp=0.0;  /* for input.group(excitonsystem).number(Temp) */
   int Nmax=0;  /* for input.group(excitonsystem).integer(Nmax) */
   int Nsteps=0;  /* for input.group(excitonsystem).integer(Nsteps) */
   double dt=0.0;  /* for input.group(excitonsystem).number(dt) */
   double tmax;
   int NumIter;
   int Npeaks=0;  /* for input.group(spectraldensity).integer(Npeaks) */
   int siteinit=0; /* for input.group(initialstate).integer(siteinit) */
   int Ei=0; /* for input.group(initialstate).integer(Ei) */
   int Ej=0; /* for input.group(initialstate).integer(Ej) */
   int delayTime=0;
   double Emin=0; /* optional argument sets energy range for spectra */
   double Emax=0; /* optional argument sets energy range for spectra*/
   int DNT=1; // initialize to 1 as default
   double norm_threshold=0.; //for transfer efficiency, propagate until population in system falls below norm_threshold
   double *valHam=NULL;
   double Jcutoff=0;
   double_complex *valrho_init=NULL;
   double_complex *valDipoleMoments=NULL;
   bool allsites_different=false;  // optional input. choose different spectral densities for the different sites
   int Nsites_coupled=0;
   vector<int> sd_CoupledSite;       // if spectral density different for all sites, define container for sd_CoupledSite,
   vector<int> sd_Npeaks;            // sd_Npeaks, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm
   vector<int> site_to_target;      // list of sites that couple to a trapping target state
   vector<int> site_to_loss;        // list of sites that couple to a trapping loss channel
   vector<double> rate_to_target;   // coupling rate to the target state
   vector<double> rate_to_loss;     // coupling rate to the loss channel
   double** sd_lambda_invcm=NULL;
   double** sd_gamma_fs=NULL;
   double** sd_gammapos_invcm=NULL;
   double     *lambda_sum;  // add to Hamiltonian
   FILE *fi=NULL;
   char fn_input[500];
   char fn_output[500];
   char calculation_mode[500];
   char method[500];
   char dipole_moments[500];
   char vendor[500];
   char device_type[500];
   char shot[500];
   sprintf(shot,"no shot specified");
   double Pol[3];
   double flux_threshold=0.0;
   Pol[0]=1.;  // set laser_polarization as default
   Pol[1]=0.;
   Pol[2]=0.;
   int GBRP=0;  /* for input.group(group_calculation_mode).group(init_2d).boolean(GBRP) */
   int SERP=0;  /* for input.group(group_calculation_mode).group(init_2d).boolean(SERP) */
   int ESARP=0;  /* for input.group(group_calculation_mode).group(init_2d).boolean(ESARP) */
   int GBNR=0;  /* for input.group(group_calculation_mode).group(init_2d).boolean(GBNR) */
   int SENR=0;  /* for input.group(group_calculation_mode).group(init_2d).boolean(SENR) */
   int ESANR=0;  /* for input.group(group_calculation_mode).group(init_2d).boolean(ESANR) */
   int HEOM=0;  /* calculation_mode=exciton eigenstate population */
   int full_Redfield=0;  /* calculation_mode=exciton eigenstate population */
   int secular_Redfield=0;  /* calculation_mode=exciton eigenstate population */
   int modified_Redfield=0;  /* calculation_mode=exciton eigenstate population */
   int generalized_Foerster=0;  /* calculation_mode=exciton eigenstate population */
   int modified_Redfield_generalized_Foerster=0;  /* calculation_mode=exciton eigenstate population */
   double sigma_pulse=0; //with of the gaussian pulse in units of seconds
   double center_pulse=0; //gaussian pulse centered in time around center_pulse in s
   double strength_pulse=0;//pulse strengt in cm^1/Debye
   double omega_pulse=0; // carrier frequency of the pulse in cm^{-1}
   if (argc!=3)
   {
      fprintf(stderr,"%s: missing input and output filenames. ABORT.",argv[0]);
      exit(1);
   }
   fi=fopen(argv[1],"r");
   if (fi==NULL)
   {
      fprintf(stderr,"%s: cannot open input parameter file \"%s\". ABORT.",argv[0],argv[1]);
      exit(1);
   }
   fflush(stdout);
   strcpy(fn_input,argv[1]);
   strcpy(fn_output,argv[2]);
   // parsing of input file
   {
      char s[500];
      while (fgets(s,500-1,fi)!=NULL)
      {
         if (strlen(s)==0) continue;
         if (s[0]=='#') continue;
         if (s[0]=='\n') continue;
         char key[500];
         sprintf(key,"QMaster-");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            for(int i=strlen(key); i<strlen(s)-1; i++) version[i-strlen(key)]=s[i];
            if (strncmp(INPUT_FORMAT_VERSION,version,strlen(INPUT_FORMAT_VERSION))!=0)
            {
               fprintf(stderr,"%s: wrong format of input file \"%s\". ABORT.\n",argv[0],argv[1]);
               fprintf(stderr,"expect %s%s but got %s%s. ABORT.\n",key,INPUT_FORMAT_VERSION, key,version);
               exit(1);
            }
            else
            {
               specify_version=true;
            }
            continue;
         }
         sprintf(key,"method=");
         if (strncmp(key,s,strlen(key))==0)
         {
            int i;
            for (i=strlen(key); i<strlen(s)-1; i++)
            {
               method[i-strlen(key)]=s[i];
            }
            method[i-strlen(key)]='\0'; // UNIX LINE ENDS REQUIRED IN INPUTFILE
            specify_method=true;
            continue;
         }
         sprintf(key,"device_id=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            use_device_id=atoi(s);
            continue;
         }
         sprintf(key,"Nsites=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            Nsites=atoi(s);
            if (valHam!=NULL) free(valHam);
            valHam=(double*)malloc(sizeof(double)*Nsites*Nsites);
            for(int k=0; k<Nsites*Nsites; k++)
            {
               valHam[k]=0.; //initialize to zero
            }

            valDipoleMoments= new double_complex[Nsites*3];
            if (valHam==NULL)
            {
               fprintf(stderr,"Malloc ERROR with valHam. ABORT.");
               exit(1);
            }
            for(int ii=0; ii<Nsites*3; ii++)
            {
               valDipoleMoments[ii]=0.;
            }
            specify_Nsites=true;
            continue;
         }
         sprintf(key,"Ham:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Nsites==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Nsites\" needs to be specified before Hamiltonian. ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            int k;
            double x;
            res=sscanf(s,"%d:%d=%lg",&i,&k,&x);
            if ((res!=3)||(i<0)||(i>=Nsites)||(k>=Nsites))
            {
               fprintf(stderr,"parse error parameter Hamiltonian \"%s\"\n", argv[1]);
               if(k>=Nsites || i>=Nsites) fprintf(stderr,"error in %s%.*s, dimension of Hamiltonian does not fit with Nsites=%i. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),Nsites);
               if((res!=3)||(i<0)) fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
            if(valHam==NULL)
            {
               fprintf(stderr,"Ham is NULL pointer. ABORT.");
               exit(1);
            }
            int id=i*Nsites+k;
            valHam[id]=x;
            specify_Ham=true;
            continue;
         }
         sprintf(key,"dipole_moments=");
         if (strncmp(key,s,strlen(key))==0)
         {
            int i;
            for (i=strlen(key); i<strlen(s)-1; i++)
            {
               dipole_moments[i-strlen(key)]=s[i];
            }
            dipole_moments[i-strlen(key)]='\0'; // UNIX LINE ENDS REQUIRED IN INPUTFILE
            specify_dipole_moments=true;
            continue;
         }
         sprintf(key,"vendor=");
         if (strncmp(key,s,strlen(key))==0)
         {
            int i;
            for (i=strlen(key); i<strlen(s)-1; i++)
            {
               vendor[i-strlen(key)]=s[i];
            }
            vendor[i-strlen(key)]='\0'; // UNIX LINE ENDS REQUIRED IN INPUTFILE
            if (strcmp(vendor,"Intel")!=0 and strcmp(vendor,"AMD")!=0 and strcmp(vendor,"NVIDIA")!=0 and strcmp(vendor,"Apple")!=0 and strcmp(vendor,"Portable Computing Language")!=0)
            {
               fprintf(stderr,"error unknown 'vendor=%s' in %s.\n",vendor, argv[1]);
               fprintf(stderr,"possible choices: vendor=Intel/AMD/NVIDIA/Apple/Portable Computing Language. ABORT.\n");
               exit(1);
            }
            specify_vendor=true;
            continue;
         }
         sprintf(key,"device_type=");
         if (strncmp(key,s,strlen(key))==0)
         {
            int i;
            for (i=strlen(key); i<strlen(s)-1; i++)
            {
               device_type[i-strlen(key)]=s[i];
            }
            device_type[i-strlen(key)]='\0'; // UNIX LINE ENDS REQUIRED IN INPUTFILE
            if (strcmp(device_type,"CPU")!=0 and strcmp(device_type,"ACCELERATOR")!=0 and strcmp(device_type,"GPU")!=0)
            {
               fprintf(stderr,"error unknown 'device_type=%s' in %s.\n",device_type, argv[1]);
               fprintf(stderr,"possible choices: device_type=CPU/ACCELERATOR/GPU. ABORT.\n");
               exit(1);
            }
            specify_device_type=true;
            continue;
         }
         sprintf(key,"dominant_path");
         if (strncmp(key,s,strlen(key))==0)
         {
            dominant_path=true;
         }
         sprintf(key,"threshold_dominant_path=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            double x;
            res=sscanf(s,"%lg",&x);
            flux_threshold=x;
            specify_threshold_dominant_path=true;
         }
         sprintf(key,"sigma_pulse=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            double x;
            res=sscanf(s,"%lg",&x);
            sigma_pulse=x;
            specify_sigma_pulse=true;
         }
         sprintf(key,"center_pulse=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            double x;
            res=sscanf(s,"%lg",&x);
            center_pulse=x;
            specify_center_pulse=true;
         }
         sprintf(key,"strength_pulse=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            double x;
            res=sscanf(s,"%lg",&x);
            strength_pulse=x*const_invcmtomeV;
            specify_strength_pulse=true;
         }
         sprintf(key,"omega_pulse=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            double x;
            res=sscanf(s,"%lg",&x);
            omega_pulse=x*const_invcmtomeV;
            specify_omega_pulse=true;
         }
         sprintf(key,"Jcutoff=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            double x;
            res=sscanf(s,"%lg",&x);
            Jcutoff=x;
            specify_Jcutoff=true;
         }
//         // Not supported in this QMaster version
//         sprintf(key,"explicit_trap_loss");
//         if (strncmp(key,s,strlen(key))==0)
//         {
//            efficiency_non_hermitian=false;
//            continue;
//         }
         sprintf(key,"dipole_vector:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_dipole_moments==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"How do you want to specify the dipole moments? Set \"dipole_moments=\" before specifying dipole_vectors. ABORT.\n");
               exit(1);
            }
            if(strncmp(dipole_moments,"dipole vectors",strlen(dipole_moments))!=0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"dipole_moments=%s\" is not compatible with \"dipole_vector:\". ABORT.\n", dipole_moments);
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            double x,y,z;
            double imx,imy, imz;
            x=0;
            y=0;
            z=0;
            imx=0;
            imy=0;
            imz=0;
            // if case 1 complex phases in dipole operators
            res=sscanf(s,"%d=%lg:%lg::%lg:%lg::%lg:%lg",&i,&x,&imx,&y,&imy,&z,&imz);
            if ((res!=7)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               parse_error=true;
            }
            else
            {
               if(valDipoleMoments==NULL)
               {
                  fprintf(stderr,"valDipoleMoments is NULL pointer. ABORT.");
                  exit(1);
               }
               valDipoleMoments[i*3+0]=x+double_complex(0,imx);
               valDipoleMoments[i*3+1]=y+double_complex(0,imy);
               valDipoleMoments[i*3+2]=z+double_complex(0,imz);
               set_dipoles.push_back(i);
               continue;
            }
            // if case 2 no complex phases in dipole operators
            res=sscanf(s,"%d=%lg:%lg:%lg",&i,&x,&y,&z);
            if ((res!=4)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               parse_error=true;
            }
            else
            {
               if(valDipoleMoments==NULL)
               {
                  fprintf(stderr,"valDipoleMoments is NULL pointer. ABORT.");
                  exit(1);
               }
               valDipoleMoments[i*3+0]=x;
               valDipoleMoments[i*3+1]=y;
               valDipoleMoments[i*3+2]=z;
               set_dipoles.push_back(i);
               continue;
            }
            if(parse_error==true)
            {
               fprintf(stderr,"parse error in %s\n",argv[1]);
               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
         }
         sprintf(key,"En_range=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            res=sscanf(s,"%lg:%lg",&Emin,&Emax);
            if ((res!=2))
            {
               fprintf(stderr,"parse error in %s\n",argv[1]);
               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
            if(Emin>Emax)
            {
               fprintf(stderr,"error in %s%.*s. Emin=%e>Emax=%e. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),Emin,Emax);
               exit(1);
            }
            specify_Enrange=true;
            continue;
         }
         sprintf(key,"return_density_matrix");
         if (strncmp(key,s,strlen(key))==0)
         {
            return_densityMatrix=true;
         }
         sprintf(key,"dipole_moment_site:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_dipole_moments==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"How do you want to specify the dipole moments? Set \"dipole_moments=\" before specifying dipole_moments. ABORT.\n");
               exit(1);
            }
            if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))!=0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"dipole_moments=%s\" is not compatible with \"dipole_moment_site:\". ABORT.\n", dipole_moments);
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            double x,y,z;
            double impart;
            y=0.;  // for effective dipole moments, dipole moment is stored in mu_x component -> use shot in x-direction
            z=0.;
            x=0.;
            impart=0.;
            // if case 1, complex phases in dipole operators
            res=sscanf(s,"%d=%lg:%lg",&i,&x,&impart);
            if ((res!=3)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               parse_error=true;
               //               fprintf(stderr,"parse error in %s\n",argv[1]);
               //               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               //               exit(1);
            }
            else
            {
               if(valDipoleMoments==NULL)
               {
                  fprintf(stderr,"valDipoleMoments is NULL pointer. ABORT.");
                  exit(1);
               }
               valDipoleMoments[i*3+0]=x+double_complex(0,impart);
               valDipoleMoments[i*3+1]=y;
               valDipoleMoments[i*3+2]=z;
               set_dipoles.push_back(i);
               continue;
            }
            // if case 2 no complex phases in dipole operators
            res=sscanf(s,"%d=%lg",&i,&x);
            if ((res!=2)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               parse_error=true;
            }
            else
            {
               if(valDipoleMoments==NULL)
               {
                  fprintf(stderr,"valDipoleMoments is NULL pointer. ABORT.");
                  exit(1);
               }
               valDipoleMoments[i*3+0]=x;
               valDipoleMoments[i*3+1]=y;
               valDipoleMoments[i*3+2]=z;
               set_dipoles.push_back(i);
               continue;
            }
            if(parse_error==true)
            {
               fprintf(stderr,"parse error in %s\n",argv[1]);
               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
         }
         sprintf(key,"dipole_moment_eigen:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_dipole_moments==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"How do you want to specify the dipole moments? Set \"dipole_moments=\" before specifying dipole_moments. ABORT.\n");
               exit(1);
            }
            if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))!=0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"dipole_moments=%s\" is not compatible with \"dipole_moment_eigen:\". ABORT.\n", dipole_moments);
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            double x,y,z;
            double impart;
            y=0.;  // for effective dipole moments, dipole moment is stored in mu_x component -> use shot in x-direction
            z=0.;
            x=0;
            impart=0;

            // if case 1, complex phases in dipole operators
            res=sscanf(s,"%d=%lg:%lg",&i,&x,&impart);
            if ((res!=3)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               parse_error=true;
            }
            else
            {
               if(valDipoleMoments==NULL)
               {
                  fprintf(stderr,"valDipoleMoments is NULL pointer. ABORT.");
                  exit(1);
               }
               valDipoleMoments[i*3+0]=x+double_complex(0,impart);
               valDipoleMoments[i*3+1]=y;
               valDipoleMoments[i*3+2]=z;
               set_dipoles.push_back(i);
               continue;
            }
            // if case 2 no complex phases in dipole operators
            res=sscanf(s,"%d=%lg",&i,&x);
            if ((res!=2)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               parse_error=true;
            }
            else
            {
               if(valDipoleMoments==NULL)
               {
                  fprintf(stderr,"valDipoleMoments is NULL pointer. ABORT.");
                  exit(1);
               }
               valDipoleMoments[i*3+0]=x;
               valDipoleMoments[i*3+1]=y;
               valDipoleMoments[i*3+2]=z;
               set_dipoles.push_back(i);
               continue;
            }
            if(parse_error==true)
            {
               fprintf(stderr,"parse error in %s\n",argv[1]);
               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
         }
         sprintf(key,"rotational_average");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_dipole_moments==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"dipole_moments=\" needs to be specified before keyword \"rotational_average\" is set. ABORT.\n");
               exit(1);
            }
            if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"rotational_average\" is not compatible with \"dipole_moments=%s\". ABORT.\n", dipole_moments);
               exit(1);
            }
            if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"rotational_average\" is not compatible with \"dipole_moments=%s\". ABORT.\n", dipole_moments);
               exit(1);
            }
            rotav=true;
            continue;
         }
         sprintf(key,"laser_polarization=");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_dipole_moments==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"dipole_moments=\" needs to be specified before keyword \"laser_polarization=\" is set. ABORT.\n");
               exit(1);
            }
            if(strncmp(dipole_moments,"effective dipole moments eigen",strlen(dipole_moments))==0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"laser_polarization=\" is not compatible with \"dipole_moments=%s\". ABORT.\n", dipole_moments);
               exit(1);
            }
            if(strncmp(dipole_moments,"effective dipole moments site",strlen(dipole_moments))==0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"laser_polarization=\" is not compatible with \"dipole_moments=%s\". ABORT.\n", dipole_moments);
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            double x,y,z;
            res=sscanf(s,"%lg:%lg:%lg",&x,&y,&z);
            if ((res!=3))
            {
               fprintf(stderr,"parse error in %s\n",res, argv[1]);
               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
            Pol[0]=x;
            Pol[1]=y;
            Pol[2]=z;
            specify_laser_polarization=true;
            continue;
         }

         sprintf(key,"Temp=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            Temp=atof(s);
            specify_Temp=true;
            continue;
         }
         sprintf(key,"Nmax=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            Nmax=atoi(s);
            specify_Nmax=true;
            continue;
         }
         sprintf(key,"NumIter=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            NumIter=atoi(s);
            if(NumIter<3)
            {
               fprintf(stderr,"error in %s\n",argv[1]);
               fprintf(stderr,"NumIter needs to be >=3. ABORT.\n");
               exit(1);
            }
            specify_NumIter=true;
            continue;
         }
         sprintf(key,"norm_threshold=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            norm_threshold=atof(s);
            continue;
         }
         sprintf(key,"Nsteps=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            Nsteps=atoi(s);
            specify_Nsteps=true;
            if(Nsteps<10)
            {
               fprintf(stderr,"error in %s\n",argv[1]);
               fprintf(stderr,"Nsteps needs to be >=10. ABORT.\n");
               exit(1);
            }
            continue;
         }
         sprintf(key,"dt=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            dt=atof(s);
            specify_dt=true;
            continue;
         }
         sprintf(key,"tmax=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            tmax=atof(s);
            if(tmax<=2.0e-14)
            {
               fprintf(stderr,"error in %s\n",argv[1]);
               fprintf(stderr,"tmax needs to be >2e-14. ABORT.\n");
               exit(1);
            }
            specify_tmax=true;
            continue;
         }
         sprintf(key,"calculation_mode=");
         if (strncmp(key,s,strlen(key))==0)
         {
            int i;
            for (i=strlen(key); i<strlen(s)-1; i++)
            {
               calculation_mode[i-strlen(key)]=s[i];
            }
            calculation_mode[i-strlen(key)]='\0'; // UNIX LINE ENDS REQUIRED IN INPUTFILE
            specify_calculation_mode=true;
            continue;
         }
         sprintf(key,"flux_analysis");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(strcmp(calculation_mode,"transfer efficiency")!=0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"flux_analysis\" can be used with \"calculation_mode=transfer efficiency\" only. ABORT.\n");
               exit(1);
            }
            flux_analysis=true;
            continue;
         }
         sprintf(key,"flux_threshold=");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(flux_analysis==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"flux_threshold\" can not be used without keyword \"flux_analysis\". ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            flux_threshold=atof(s);
            continue;
         }
         sprintf(key,"padding=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            padding=atoi(s);
            continue;
         }
         sprintf(key,"shot=");
         if (strncmp(key,s,strlen(key))==0)
         {
            int i;
            for (i=strlen(key); i<strlen(s)-1; i++)
            {
               shot[i-strlen(key)]=s[i];
            }
            shot[i-strlen(key)]='\0'; // FIXME UNIX LINE ENDS REQUIRED
            specify_shot=true;
            continue;
         }
         sprintf(key,"general_initial_state");
         if (strncmp(key,s,strlen(key))==0)
         {
            if (valrho_init!=NULL) free(valrho_init);
            valrho_init=(double_complex*)malloc(sizeof(double_complex)*Nsites*Nsites);
            for(int k=0; k<Nsites*Nsites; k++)
            {
               valrho_init[k]=double_complex(0.,0); //initialize to zero
            }
            specify_general_initial_state=true;
         }
         sprintf(key,"rho_init:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Nsites==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Nsites\" needs to be specified before rho_init. ABORT.\n");
               exit(1);
            }
            if(specify_general_initial_state==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"general_initial_state\" needs to be specified before rho_init. ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            int k;
            double Re_x;
            double Im_x;
            res=sscanf(s,"%d:%d=%lg:%lg",&i,&k,&Re_x,&Im_x);
            if ((res!=4)||(i<0)||(i>=Nsites)||(k>=Nsites))
            {
               fprintf(stderr,"parse error parameter rho_init \"%s\"\n", argv[1]);
               if(k>=Nsites || i>=Nsites) fprintf(stderr,"error in %s%.*s, dimension of rho_init does not fit with Nsites=%i. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),Nsites);
               if((res!=3)||(i<0)) fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
            if(valrho_init==NULL)
            {
               fprintf(stderr,"rho_init is NULL pointer. ABORT.");
               exit(1);
            }
            int id=i*Nsites+k;
            valrho_init[id]=double_complex(Re_x,Im_x);
            specify_valrho_init=true;
            continue;
         }
         sprintf(key,"siteinit=");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_general_initial_state==true)
            {
               fprintf(stderr,"error in %s\n",argv[1]);
               fprintf(stderr,"\"siteinit\" cannot be specified toghether with \"general_initial_state\". ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            siteinit=atoi(s);
            specify_siteinit=true;
            continue;
         }
         sprintf(key,"Ei=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            Ei=atoi(s);
            specify_Ei=true;
            continue;
         }
         sprintf(key,"Ej=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            Ej=atoi(s);
            specify_Ej=true;
            continue;
         }
         sprintf(key,"delayTime=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            delayTime=atoi(s);
            continue;
         }
         sprintf(key,"DNT=");
         if (strncmp(key,s,strlen(key))==0)
         {
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            DNT=atoi(s);
            specify_DNT=true;
            continue;
         }
         sprintf(key,"allsites_different");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Npeaks==true)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"allsites_different\" is in conflict with specification \"Npeaks\". ABORT.\n");
               exit(1);
            }
            allsites_different=true;
            continue;
         }
         sprintf(key,"Npeaks=");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Nsites==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Nsites\" needs to be specified before \"Npeaks\". ABORT.\n");
               exit(1);
            }
            if(allsites_different==true)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Npeaks\" is in conflict with specification \"allsites_different\". ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            Npeaks=atoi(s);
            sd_lambda_invcm=new double*[Nsites];
            sd_gamma_fs=new double*[Nsites];
            sd_gammapos_invcm=new double*[Nsites];
            for(int i=0; i<Nsites; i++)
            {
               sd_CoupledSite.push_back(i);
               sd_Npeaks.push_back(Npeaks);
               sd_lambda_invcm[i]=new double[Npeaks];
               sd_gamma_fs[i]=new double[Npeaks];
               sd_gammapos_invcm[i]=new double[Npeaks];
            }
            if ((sd_lambda_invcm==NULL)||(sd_gamma_fs==NULL)||(sd_gammapos_invcm==NULL))
            {
               fprintf(stderr,"Malloc ERROR with Npeaks. ABORT.");
               exit(1);
            }
            specify_Npeaks=true;
            continue;
         }

         sprintf(key,"Nsites_coupled=");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Nsites==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Nsites\" needs to be specified before \"Nsites_coupled\". ABORT.\n");
               exit(1);
            }
            if(allsites_different==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"Nsites_coupled\" cannot be defined without \"allsites_different\". ABORT.\n");
               fprintf(stderr,"\"allsites_different\" needs to be specified before \"Nsites_coupled\". ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            sd_Npeaks.clear();
            sd_CoupledSite.clear();
            string str;
            int site;
            int count=0;
            for(int k=strlen(key); k< strlen(s); k++)
            {
               if(s[k]==':')
               {
                  // cout<<str<<endl;
                  site=atoi(str.c_str());
                  sd_CoupledSite.push_back(site);
                  sd_readin_peaks.push_back(0);
                  if(readin_sites.count(site) != 0)
                  {
                     fprintf(stderr,"parse error %s\n",argv[1]);
                     fprintf(stderr,"Error Nsites_coupled: site %d already couples to the phonon bath. ABORT.\n", site);
                     exit(1);
                  };
                  readin_sites[site]=count;
                  if(site>=Nsites)
                  {
                     fprintf(stderr,"parse error %s\n",argv[1]);
                     fprintf(stderr,"Error Nsites_coupled: specified site %d is not available, max. site index is given by %d. ABORT.\n", site, Nsites-1);
                     exit(1);
                  }
                  count+=1;
                  str.clear();
               }
               else str+=s[k];
            }
            site=atoi(str.c_str());
            if(site>=Nsites)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"Error: specified site %d is not available, max. site index is given by %d. ABORT.\n", site, Nsites-1);
               exit(1);
            }
            if(readin_sites.count(site) != 0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"Error Nsites_coupled: site %d already couples to the phonon bath. ABORT.\n", site);
               exit(1);
            };
            readin_sites[site]=count;
            sd_CoupledSite.push_back(site);
            sd_readin_peaks.push_back(0);
            Nsites_coupled=readin_sites.size();
            if(Nsites_coupled>Nsites)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"number of coupled sites\"=%d exceeds number of sites \"Nsites=%d\".\n", Nsites_coupled, Nsites);
               fprintf(stderr,"You cannot couple more sites to the bath than are available in your system. ABORT.\n", Nsites_coupled, Nsites);
               exit(1);
            }
            for(int i=0; i<Nsites_coupled; i++)
            {
               fprintf(stdout,"# Add site %i to phonon bath\n", sd_CoupledSite[i]);
            }
            specify_Nsites_coupled=true;
            continue;
         }
         sprintf(key,"Npeaks_sd=");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(allsites_different==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"Npeaks_sd\" cannot be defined without \"allsites_different\". ABORT.\n");
               fprintf(stderr,"\"allsites_different\" needs to be specified before \"Npeaks_sd\". ABORT.\n");
               exit(1);
            }

            if(specify_Npeaks_all==true)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"Npeaks_sd\" and \"Npeaks_all\" cannot be defined toghether. ABORT.\n");
               exit(1);
            }

            if(specify_Nsites_coupled==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Nsites_coupled\" needs to be specified before \"Npeaks_sd\". ABORT.\n");
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            string str;
            int NDL;
            int count=0;
            for(int k=strlen(key); k< strlen(s); k++)
            {
               if(s[k]==':')
               {
                  NDL=atoi(str.c_str());
                  sd_Npeaks.push_back(NDL);
                  count+=1;
                  str.clear();
               }
               else str+=s[k];
            }
            NDL=atoi(str.c_str());
            sd_Npeaks.push_back(NDL);
            if(sd_Npeaks.size()!=Nsites_coupled)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: \"number of given \"sd_Npeaks entries\" %d != \"number of coupled sites\" %d\".\n", sd_Npeaks.size(),Nsites_coupled);
               exit(1);
            }
            sd_lambda_invcm=new double*[Nsites_coupled];
            sd_gamma_fs=new double*[Nsites_coupled];
            sd_gammapos_invcm=new double*[Nsites_coupled];
            for(int k=0; k<Nsites_coupled; k++)
            {
               sd_lambda_invcm[k]=new double[sd_Npeaks[k]];
               sd_gamma_fs[k]=new double[sd_Npeaks[k]];
               sd_gammapos_invcm[k]=new double[sd_Npeaks[k]];
            }
            if ((sd_lambda_invcm==NULL)||(sd_gamma_fs==NULL)||(sd_gammapos_invcm==NULL))
            {
               fprintf(stderr,"Malloc ERROR with Npeaks. ABORT.");
               exit(1);
            }
            specify_Npeaks_all=true;
            continue;
         }
         sprintf(key,"sd_peak:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(allsites_different==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"conflict: spectral densities for individual peaks found \"sd_peak:i:j=\", cannot be defined without \"allsites_different\". ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            int site;
            double x,y,z;
            if(specify_Npeaks_all==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Npeaks_all\" or \"Npeaks_sd\" needs to be specified before the parameter of the spectral density. ABORT.\n");
               exit(1);
            }
            res=sscanf(s,"%d:%d=%lg:%lg:%lg",&site, &i,&x,&y,&z);
            if(readin_sites.count(site) == 0)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"Error in sd_peak:%i:%i=%lg:%lg:%lg, site %d does not couple to the phonon bath. ABORT.\n",site,i,x,y,z, site);
               exit(1);
            };
            int NDL=sd_Npeaks[readin_sites[site]];
            if ((res!=5)||(i<0)||(i>=NDL))
            {
               if(i>=NDL)
               {
                  fprintf(stderr,"parse error for parameter of the spectral density \"%s\"\n",argv[1]);
                  fprintf(stderr,"error: for site %d more than %d-peaks are specified in the input file.\nThe number of peaks are defined by \"Npeaks_sd\" or \"Npeaks_all\". ABORT.\n",site, NDL);
                  exit(1);
               }
               else
               {
                  fprintf(stderr,"parse error for parameter of the spectral density \"%s\"\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
                  exit(1);
               }
            }
            sd_lambda_invcm[readin_sites[site]][i]=x;
            sd_gamma_fs[readin_sites[site]][i]=y;
            sd_gammapos_invcm[readin_sites[site]][i]=z;
            sd_readin_peaks[readin_sites[site]]+=1;
            readin_all=0;
            for(int kk=0; kk<Nsites_coupled; kk++)
            {
               if(sd_readin_peaks[kk]==sd_Npeaks[kk])
               {
                  readin_all+=1;
               }
            }
            if(readin_all==Nsites_coupled)
            {
               specify_peaks=true;
            }
            continue;
         }
         sprintf(key,"sd_peak_all:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Npeaks==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Npeaks\" needs to be specified before the parameter of the spectral density. ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            double x,y,z;

            res=sscanf(s,"%d=%lg:%lg:%lg",&i,&x,&y,&z);
            if ((res!=4)||(i<0)||(i>=Npeaks))
            {
               if(i>=Npeaks)
               {
                  fprintf(stderr,"parse error for parameter of the spectral density \"%s\"\n",argv[1]);
                  fprintf(stderr,"error more than %d-peaks are specified in the input file. ABORT.\n",Npeaks);
                  exit(1);
               }
               else
               {
                  fprintf(stderr,"parse error for parameter of the spectral density \"%s\"\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
                  exit(1);
               }
            }
            for(int k=0; k<Nsites; k++)
            {
               sd_lambda_invcm  [k][i]=x;
               sd_gamma_fs      [k][i]=y;
               sd_gammapos_invcm[k][i]=z;
            }
            readin_peaks+=1;
            set_sd_peak_all.push_back(i);
            if(readin_peaks==Npeaks)
            {
               specify_peaks=true;
            }
            continue;
         }
         sprintf(key,"GBRP");
         if (strncmp(key,s,strlen(key))==0)
         {
            GBRP=1;
            specify_2dpathway=true;
            continue;
         }
         sprintf(key,"GBNR");
         if (strncmp(key,s,strlen(key))==0)
         {
            GBNR=1;
            specify_2dpathway=true;
            continue;
         }
         sprintf(key,"SERP");
         if (strncmp(key,s,strlen(key))==0)
         {
            SERP=1;
            specify_2dpathway=true;
            continue;
         }
         sprintf(key,"SENR");
         if (strncmp(key,s,strlen(key))==0)
         {
            SENR=1;
            specify_2dpathway=true;
            continue;
         }
         sprintf(key,"ESARP");
         if (strncmp(key,s,strlen(key))==0)
         {
            ESARP=1;
            specify_2dpathway=true;
            continue;
         }
         sprintf(key,"ESANR");
         if (strncmp(key,s,strlen(key))==0)
         {
            ESANR=1;
            specify_2dpathway=true;
            continue;
         }

         sprintf(key,"HEOM");
         if (strncmp(key,s,strlen(key))==0)
         {
            HEOM=1;
            specify_method_eigenstate_population=true;
            continue;
         }
         sprintf(key,"full_Redfield");
         if (strncmp(key,s,strlen(key))==0)
         {
            full_Redfield=1;
            specify_method_eigenstate_population=true;
            continue;
         }
         sprintf(key,"secular_Redfield");
         if (strncmp(key,s,strlen(key))==0)
         {
            secular_Redfield=1;
            specify_method_eigenstate_population=true;
            continue;
         }
         sprintf(key,"modified_Redfield");
         if (strncmp(key,s,strlen(key))==0)
         {
            modified_Redfield=1;
            specify_method_eigenstate_population=true;
            continue;
         }
         sprintf(key,"generalized_Foerster");
         if (strncmp(key,s,strlen(key))==0)
         {
            generalized_Foerster=1;
            specify_method_eigenstate_population=true;
            continue;
         }

         sprintf(key,"combined_modified_Redfield_generalized_Foerster");
         if (strncmp(key,s,strlen(key))==0)
         {
            modified_Redfield_generalized_Foerster=1;
            specify_method_eigenstate_population=true;
            continue;
         }
         sprintf(key,"link_to_target:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Nsites==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Nsites\" needs to be specified before \"link_to_target\". ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            double Gamma;
            res=sscanf(s,"%d=%lg",&i,&Gamma);
            if ((res!=2)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               fprintf(stderr,"parse error in %s\n",argv[1]);
               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
            site_to_target.push_back(i);
            rate_to_target.push_back(Gamma);
         }
         sprintf(key,"link_to_loss:");
         if (strncmp(key,s,strlen(key))==0)
         {
            if(specify_Nsites==false)
            {
               fprintf(stderr,"parse error %s\n",argv[1]);
               fprintf(stderr,"\"Nsites\" needs to be specified before \"link_to_loss\". ABORT.\n");
               exit(1);
            }
            for (int i=0; i<strlen(key); i++) s[i]=' ';
            int res;
            int i;
            double Gamma;
            res=sscanf(s,"%d=%lg",&i,&Gamma);
            if ((res!=2)||(i<0)||(i>=Nsites))
            {
               if(i>=Nsites)
               {
                  fprintf(stderr,"parse error in %s\n",argv[1]);
                  fprintf(stderr,"syntax error in %s%.*s: %i>=Nsites. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key),i);
                  exit(1);
               }
               fprintf(stderr,"parse error in %s\n",argv[1]);
               fprintf(stderr,"syntax error in %s%.*s. ABORT.\n",key,strlen(s)-strlen(key)-1,s+strlen(key));
               exit(1);
            }
            site_to_loss.push_back(i);
            rate_to_loss.push_back(Gamma);
         }
      }
      fclose(fi);
   }

   //Begin: Do some basic checks if input is correct
   fprintf(stdout,"# check input ....\n");
   if(specify_version==false)
   {
      fprintf(stderr,"error no QMaster- version specified in %s. ABORT.\n",argv[1]);
      exit(1);
   }
   if(specify_vendor==false)
   {
      fprintf(stderr,"error no vendor specified in %s. ABORT.\n",argv[1]);
      fprintf(stderr,"possible choices: vendor=Intel/AMD/NVIDIA. ABORT.\n");
      exit(1);
   }
   if(specify_device_type==false)
   {
      fprintf(stderr,"error no device_type specified in %s\n",argv[1]);
      fprintf(stderr,"possible choices: device_type=CPU/ACCELERATOR/GPU. ABORT.\n");
      exit(1);
   }





   //FIXHME tests

   if((strcmp(calculation_mode,"exciton eigenstate population")==0))
   {
      if(specify_method==true)
      {
         fprintf(stderr,"conflict in %s:\nkeyword \"method=\" cannot be specified in calculation mode \"exciton eigenstate population\".\nPropagation method needs to be set by keywords.\nPossile keywords are: HEOM/full_Redfield/secular_Redfield/modified_Redfield/generalized_Foerster/combined_modified_Redfield_generalized_Foerster.  ABORT.\n",argv[1]);
         exit(1);
      }
      if(specify_method_eigenstate_population==false)
      {
         fprintf(stderr,"error no propagation method specified in %s.\nPropagation method needs to be set by keywords.\nPossile keywords are: HEOM/full_Redfield/secular_Redfield/modified_Redfield/generalized_Foerster/combined_modified_Redfield_generalized_Foerster. ABORT.\n",argv[1]);
         exit(1);
      }
      if(HEOM==0 and specify_Nmax==true)
      {
         fprintf(stderr,"parse error %s\n",argv[1]);
         fprintf(stderr,"conflict: \"Nmax\" is only required for HEOM. ABORT.\n");
         exit(1);
      }
      if(modified_Redfield==1)
      {
         if(specify_tmax==false)
         {
            fprintf(stderr,"parse error %s\n",argv[1]);
            fprintf(stderr,"method modified_Redfield requires \"tmax\" to be specified. ABORT.\n");
            exit(1);
         }
         if(specify_NumIter==false)
         {
            fprintf(stderr,"parse error %s\n",argv[1]);
            fprintf(stderr,"method modified_Redfield requires \"NumIter\" to be specified. ABORT.\n");
            exit(1);
         }
      }
      if(modified_Redfield_generalized_Foerster==1)
      {
         if(specify_tmax==false)
         {
            fprintf(stderr,"parse error %s\n",argv[1]);
            fprintf(stderr,"method combined_modified_Redfield_generalized_Foerster requires \"tmax\" to be specified. ABORT.\n");
            exit(1);
         }
         if(specify_NumIter==false)
         {
            fprintf(stderr,"parse error %s\n",argv[1]);
            fprintf(stderr,"method combined_modified_Redfield_generalized_Foerster requires \"NumIter\" to be specified. ABORT.\n");
            exit(1);
         }
         if(specify_Jcutoff==false)
         {
            fprintf(stderr,"parse error %s\n",argv[1]);
            fprintf(stderr,"method combined_modified_Redfield_generalized_Foerster requires \"Jcutoff\" to be specified. ABORT.\n");
            exit(1);
         }
      }
      if(generalized_Foerster==1)
      {
         if(specify_tmax==false)
         {
            fprintf(stderr,"parse error %s\n",argv[1]);
            fprintf(stderr,"method generalized_Foerster requires \"tmax\" to be specified. ABORT.\n");
            exit(1);
         }
      }
      if(specify_Nmax==false and HEOM==1)
      {
         fprintf(stderr,"error no Nmax specified in %s. ABORT.\n",argv[1]);
         exit(1);
      }
   }
   else
   {
      if(specify_method==false)
      {
         fprintf(stderr,"error no propagation method (for example method=HEOM) specified in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(strcmp(method,"HEOM")!=0 and strcmp(method,"full Redfield")!=0 and strcmp(method,"secular Redfield")!=0)
      {

         fprintf(stderr,"error: method=\"%s\" invalid propagation method (allowed: method=HEOM/full Redfield/secular Redfield) in %s. ABORT.\n",method, argv[1]);
         exit(1);
      }
      if(strcmp(method,"HEOM")!=0 and specify_Nmax==true)
      {
         fprintf(stderr,"parse error %s\n",argv[1]);
         fprintf(stderr,"conflict: \"Nmax\" is only required for \"method=HEOM\". ABORT.\n");
         exit(1);
      }
      if(specify_method_eigenstate_population==true)
      {
         fprintf(stderr,"conflict in %s:\nkeywords HEOM/full_Redfield/secular_Redfield/modified_Redfield/generalized_Foerster/combined_modified_Redfield_generalized_Foerster\ncannot be specified in calculation mode \"%s\".  ABORT.\n",argv[1],calculation_mode);
         exit(1);
      }
      if(specify_Nmax==false and (strcmp(method,"HEOM")==0))
      {
         fprintf(stderr,"error no Nmax specified in %s. ABORT.\n",argv[1]);
         exit(1);
      }
   }
//END FIXHME TESTS





   if(specify_peaks==true)
   {
      for(int i=0; i<Nsites; i++)
      {
         int check=count(set_sd_peak_all.begin(), set_sd_peak_all.end(),i);
         if(check>1)
         {
            fprintf(stderr,"conflict in %s:\nsd_peak_all: chromophore %i is defined multiple times. ABORT.\n",argv[1],i);
            exit(1);
         }
      }
   }
   if(specify_peaks==false)
   {
      if(allsites_different==false)
      {
         fprintf(stderr,"error in \"%s\"\n",argv[1]);
         fprintf(stderr,"error %d of %d-peaks of the spectral density are specified. ABORT. \n", readin_peaks, Npeaks);
         exit(1);
      }
      if(allsites_different==true)
      {
         fprintf(stderr,"error in \"%s\"\n",argv[1]);
         for(int i=0; i<Nsites_coupled; i++)
         {
            if(sd_readin_peaks[i]!=sd_Npeaks[i])
            {
               fprintf(stderr,"error in the spectral density of site %d: %d of %d-peaks of the spectral density are specified. ABORT. \n", sd_CoupledSite[i], sd_readin_peaks[i], sd_Npeaks[i]);
            }
         }
         exit(1);
      }
   }
   if(specify_Nsites==false)
   {
      fprintf(stderr,"error no Nsites specified in %s. ABORT.\n",argv[1]);
      exit(1);
   }
   if(specify_Ham==false)
   {
      fprintf(stderr,"error no Hamiltonian specified in %s. ABORT.\n",argv[1]);
      exit(1);
   }
   if(strcmp(calculation_mode,"population dynamics")==0 or strcmp(calculation_mode,"transfer efficiency")==0)
   {
      if(specify_valrho_init==false && specify_siteinit==false)
      {
         fprintf(stderr,"error no \"siteinit\" or \"rho_init\" specified in %s. ABORT.\n",argv[1]);
         exit(1);
      }
   }
   if(specify_Temp==false)
   {
      fprintf(stderr,"error no Temp specified in %s. ABORT.\n",argv[1]);
      exit(1);
   }
   if(specify_Nsteps==false)
   {
      fprintf(stderr,"error no Nsteps specified in %s. ABORT.\n",argv[1]);
      exit(1);
   }
   if(specify_dt==false)
   {
      fprintf(stderr,"error no dt specified in %s. ABORT.\n",argv[1]);
      exit(1);
   }
   if(specify_calculation_mode==false)
   {
      fprintf(stderr,"error no calculation mode specified in %s. ABORT.\n",argv[1]);
      exit(1);
   }
   if((strcmp(calculation_mode,"applyPulse")==0) )
   {
      if(specify_sigma_pulse==false)
      {
         fprintf(stderr,"error in %s\n",argv[1]);
         fprintf(stderr,"\"sigma_pulse\" is not specified. ABORT.\n");
         exit(1);
      }
      if(specify_center_pulse==false)
      {
         fprintf(stderr,"error in %s\n",argv[1]);
         fprintf(stderr,"\"center_pulse\" is not specified. ABORT.\n");
         exit(1);
      }
      if(specify_strength_pulse==false)
      {
         fprintf(stderr,"error in %s\n",argv[1]);
         fprintf(stderr,"\"strength_pulse\" is not specified. ABORT.\n");
         exit(1);
      }
      if(specify_omega_pulse==false)
      {
         fprintf(stderr,"error in %s\n",argv[1]);
         fprintf(stderr,"\"omega_pulse\" is not specified. ABORT.\n");
         exit(1);
      }
   }


   if((strcmp(calculation_mode,"population dynamics")==0 or strcmp(calculation_mode,"transfer efficiency")==0))
   {
      if(specify_siteinit==true && specify_general_initial_state==true)
      {
         fprintf(stderr,"error in %s\n",argv[1]);
         fprintf(stderr,"\"siteinit\" cannot be specified toghether with \"general_initial_state\". ABORT.\n");
         exit(1);
      }
   }
   if(specify_DNT==true and (strcmp(calculation_mode,"absorption spectra")==0))
   {
      fprintf(stdout,"# warning: ignore DNT=%i; DNT is not supported for calculation mode = %s\n",DNT,calculation_mode);
   }
   if(strcmp(calculation_mode,"exciton coherence")==0)
   {
      if(specify_Ei==false)
      {
         fprintf(stderr,"error no value for Ei specified in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(specify_Ej==false)
      {
         fprintf(stderr,"error no value for Ej specified in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(Ei<0)
      {
         fprintf(stderr,"error value for Ei<0 in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(Ei>=Nsites)
      {
         fprintf(stderr,"error value for Ei>=Nsites in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(Ej<0)
      {
         fprintf(stderr,"error value for Ej<0 in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(Ej>=Nsites)
      {
         fprintf(stderr,"error value for Ej>=Nsites in %s. ABORT.\n",argv[1]);
         exit(1);
      }
   }
   if(strcmp(calculation_mode,"exciton eigenstate population")==0)
   {
      if(specify_Ei==false)
      {
         fprintf(stderr,"error no value for Ei specified in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(specify_Ej==true)
      {
         fprintf(stderr,"conflict in %s:\nEj cannot be specified in calculation mode \"exciton eigenstate population\". ABORT.\n",argv[1]);
         exit(1);
      }
      if(Ei<0)
      {
         fprintf(stderr,"error value for Ei<0 in %s. ABORT.\n",argv[1]);
         exit(1);
      }
      if(Ei>=Nsites)
      {
         fprintf(stderr,"error value for Ei>=Nsites in %s. ABORT.\n",argv[1]);
         exit(1);
      }
   }
   if(     strcmp(calculation_mode,"absorption spectra")==0
           or strcmp(calculation_mode,"two-dimensional spectra")==0
           or strcmp(calculation_mode,"pump-probe witness")==0
           or strcmp(calculation_mode,"absorption spectra finite pulse")==0
     )
   {
      // check if all dipole moments are defined
      for(int i=0; i<Nsites; i++)
      {
         int check=count(set_dipoles.begin(), set_dipoles.end(),i);
         if(check==0)
         {
            fprintf(stderr,"error in %s:\ndipole vector or effective dipole moment for chromophore %i is not defined. ABORT.\n",argv[1],i);
            exit(1);
         }
         if(check>1)
         {
            fprintf(stderr,"conflict in %s:\ndipole vector or effective dipole moment for chromophore %i is defined multiple times. ABORT.\n",argv[1],i);
            exit(1);
         }
      }
      if(specify_laser_polarization==true and rotav==true)
      {
         fprintf(stderr,"conflict in %s:\ncannot specify \"laser_polarization\" together with \"rotational_average\". ABORT.\n",argv[1]);
         exit(1);
      }
   }
   if(strcmp(calculation_mode,"two-dimensional spectra")==0 or strcmp(calculation_mode,"pump-probe witness")==0)
   {
      if(specify_2dpathway==false)
      {
         fprintf(stderr,"error in %s:\nno pathway specified (choices are: SERP/SENR/GBRP/GBNR/ESARP/EASNR).  ABORT.\n",argv[1]);
         exit(1);
      }


      if(specify_laser_polarization==true and specify_shot==true)
      {
         fprintf(stderr,"conflict in %s:\ncannot specify \"laser_polarization\" together with \"shot=\". ABORT.\n",argv[1]);
         exit(1);
      }
      if(strcmp(dipole_moments,"effective dipole moments site")==0 and specify_shot==true)
      {
         fprintf(stderr,"conflict in %s:\ncannot specify \"dipole_moments=effective dipole moments site\" together with \"shot=\". ABORT.\n",argv[1]);
         exit(1);
      }
      if(strcmp(dipole_moments,"effective dipole moments eigen")==0 and specify_shot==true)
      {
         fprintf(stderr,"conflict in %s:\ncannot specify \"dipole_moments=effective dipole moments eigen\" together with \"shot=\". ABORT.\n",argv[1]);
         exit(1);
      }
   }
   if(strcmp(calculation_mode,"absorption spectra finite pulse")==0)
   {
      if(specify_sigma_pulse==false)
      {
         fprintf(stderr,"error in %s: \"sigma_pulse\" is not specified. ABORT.\n",argv[1]);
         exit(1);
      }
   }

   if(dominant_path==true and specify_threshold_dominant_path==false)
   {
      fprintf(stderr,"error in %s: if keyword \"dominant_path\" is set one needs to specify \"threshold_dominant_path=\". ABORT.\n",argv[1]);
      exit(1);
   }


   if(strcmp(calculation_mode,"transfer efficiency")==0)
   {
      for(int i=0; i<Nsites; i++)
      {
         int check=count(site_to_target.begin(), site_to_target.end(),i);
         if(check>1)
         {
            fprintf(stderr,"conflict in %s:\nChromophore %i already couples to the target state. ABORT.\n",argv[1],i);
            exit(1);
         }
      }
      for(int i=0; i<Nsites; i++)
      {
         int check=count(site_to_loss.begin(), site_to_loss.end(),i);
         if(check>1)
         {
            fprintf(stderr,"conflict in %s:\nChromophore %i already couples to the loss channel. ABORT.\n",argv[1],i);
            exit(1);
         }
      }
   }
   if(strcmp(calculation_mode,"absorption spectra")==0
         or strcmp(calculation_mode,"two-dimensional spectra")==0
         or strcmp(calculation_mode,"pump-probe witness")==0
         or strcmp(calculation_mode,"absorption spectra finite pulse")==0
         or strcmp(calculation_mode,"exciton eigenstate population")==0)
   {
      if(return_densityMatrix==true)
      {
         fprintf(stderr,"conflict in %s:\nreturn_density_matrix cannot be specified for calculation_mode='%s'. ABORT.\n",argv[1],calculation_mode);
         exit(1);
      }
   }
   fprintf(stdout,"# .... done\n");
   //End: Do some basic checks if input is correct
   // get number of coupled sites and LD peaks
   for(int i=0; i<sd_Npeaks.size(); i++)
   {
      Ncoupled+=sd_Npeaks[i];
   }
   Ncoupled*=2;
   // cout<<"Ncoupled="<<Ncoupled<<endl;
   // check if only 1 peak for all sites
   // special case: non-shifted single Drude-Lorentz spectral density
   // here only one peak
   for(int i=0; i<sd_Npeaks.size(); i++)
   {
      if(sd_Npeaks[i]==1 and sd_gammapos_invcm[i][0]==0.0)
      {
         Ncoupled-=1;  // correct for double counting
      }
   }
   // calculate total reorganization energy
   lambda_sum=new double[Nsites];
   for(int i=0; i<Nsites; i++)
   {
      lambda_sum[i]=0.;
   }
   for(int i=0; i<sd_CoupledSite.size(); i++)
   {
      // loop over all coupled sites
      int site=sd_CoupledSite[i];
      lambda_sum[site]=0.;
      for(int k=0; k<sd_Npeaks[i]; k++)
      {
         // loop over NDLpeaks per of coupled site=sd_CoupledSite[i]
         lambda_sum[site]+=sd_lambda_invcm[i][k];
      }
      lambda_sum[site]*=const_invcmtomeV;
   }
   // initialize Hamiltonian
   Matrix *HamMat=new Matrix(Nsites);
   {
      for(int i=0; i<Nsites; i++)
      {
         for(int j=0; j<Nsites; j++)
         {
            HamMat->assignEntry(i,j,valHam[i*Nsites+j]*const_invcmtomeV);
         }
      }
   }

/// Initialize openCL platforms
   OpenCL_Init OpenCLinfo;
   OpenCLinfo.get_platforms();
   OpenCLinfo.initialize(vendor,device_type, use_device_id);
   if(strcmp(device_type,"CPU")==0)
   {
      OpenCLinfo.get_kernels("LoopIn");   // For kernel execution of Hexciton use LoopIn kernel version, CPUs give better performance for this version
   }
   else if(strcmp(device_type,"ACCELERATOR")==0)
   {
      if(strcmp(vendor,"AMD")==0)
      {
         OpenCLinfo.get_kernels("default"); //use default configuration if AMD ACCELERATOR (not tested but expected that Radeon etc. behaves like NVIDIA GPU)
      }
      else
      {
         OpenCLinfo.get_kernels("LoopIn");   // For kernel execution of Hexciton use LoopIn kernel version, Intel Xeon Phi gives better performance for this version
      }
   }
   else if(strcmp(device_type,"GPU")==0)
   {
      OpenCLinfo.get_kernels("default");
   }
   // depending on calculation mode, initialize system
//    if(strcmp(calculation_mode,"absorption spectra")==0)
//    {
//       // initialize the system
//       System system;
//       system.initialize_spectra(
//          Nsites,            // total number of sites in Hamiltonian
//          Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//          const_invcmtomeV,  // convert cm^-1 to meV
//          const_kb,          // bolzmann constant in meV/K
//          const_hbar,  	    // hbar in meV
//          Temp, 		    // Temperature in K
//          lambda_sum,        // lambda_sum at each site already in meV
//          HamMat,      	    // Hamiltonian already in meV
//          valDipoleMoments,  // 3-vector of dipole moment at each site
//          Pol                // 3-vector or laser polarization
//       );
//       calculate_1d absorption_spectrum(system, OpenCLinfo, Nmax, dt, Nsteps, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, padding, rotav, Pol, fn_output, dipole_moments, Emin, Emax, specify_Enrange);
//       if(strcmp(method,"HEOM")==0)
//       {
//          absorption_spectrum.execute_calculate_1d();
//       }
//       if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//       {
//          absorption_spectrum.execute_calculate_Redfield_1d(method);
//       }
//       compute_time=absorption_spectrum.diff;
//    }


   // depending on calculation mode, initialize system
//   if(strcmp(calculation_mode,"applyPulse")==0)
//   {
//      // initialize the system
//      System system;
//      system.initialize_spectra(
//         Nsites,            // total number of sites in Hamiltonian
//         Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//         const_invcmtomeV,  // convert cm^-1 to meV
//         const_kb,          // bolzmann constant in meV/K
//         const_hbar,  	    // hbar in meV
//         Temp, 		    // Temperature in K
//         lambda_sum,        // lambda_sum at each site already in meV
//         HamMat,      	    // Hamiltonian already in meV
//         valDipoleMoments,  // 3-vector of dipole moment at each site
//         Pol                // 3-vector or laser polarization
//      );
////      cout<<"sigma"<<sigma_pulse<<endl;
////      cout<<"center"<<center_pulse<<endl;
////      cout<<"strength"<<strength_pulse/const_invcmtomeV<<endl;
////      cout<<"freq"<<omega_pulse/const_invcmtomeV<<endl;
//
//      calculate_applyPulse applyPulse(system, OpenCLinfo, Nmax, dt, DNT, Nsteps, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, sigma_pulse, center_pulse, strength_pulse, omega_pulse, padding, rotav, Pol, fn_output, dipole_moments, Emin, Emax, specify_Enrange, return_densityMatrix);
//      if(strcmp(method,"HEOM")==0)
//      {
//         applyPulse.execute_calculate_applyPulse();
//      }
//      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//      {
//         applyPulse.execute_calculate_Redfield_applyPulse(method);
//      }
//      compute_time=applyPulse.diff;
//   }

   /*
      // depending on calculation mode, initialize system
      else if(strcmp(calculation_mode,"absorption spectra finite pulse")==0)
      {
         int padding=4; // FIXME: should be a user configuration, only used for calculating FFTs to boost resolution
         // initialize the system
         System system;
         system.initialize_spectra(
            Nsites,            // total number of sites in Hamiltonian
            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
            const_invcmtomeV,  // convert cm^-1 to meV
            const_kb,          // bolzmann constant in meV/K
            const_hbar,  	    // hbar in meV
            Temp, 		    // Temperature in K
            lambda_sum,        // lambda_sum at each site already in meV
            HamMat,      	    // Hamiltonian already in meV
            valDipoleMoments,  // 3-vector of dipole moment at each site
            Pol                // 3-vector or laser polarization
         );
         calculate_1d_finite_pulse absorption_spectrum_finite_pulse(system, Nmax, dt, DNT, Nsteps, sigma_pulse, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, padding, rotav, Pol, fn_output, dipole_moments);
         absorption_spectrum_finite_pulse.execute_calculate_1d_finite_pulse();
         GPUtime=absorption_spectrum_finite_pulse.diff;
      }
   */
//   else if(strcmp(calculation_mode,"two-dimensional spectra")==0)
//   {
//      // enforce that Nsteps3 is a multiple of DNT
//      int Nsteps3=int(Nsteps/DNT+0.5)*DNT; // enforce that Nsteps3 is a multiple of DNT
//      if (Nsteps3!=Nsteps)
//      {
//         fprintf(stdout,"# enforce Nsteps3-1 to fit with multiples of DNT=%d\n",DNT);
//         fprintf(stdout,"# set Nsteps3=%d (instead of %d)\n",Nsteps3,Nsteps);
//      }
//      //
//      if(strcmp(dipole_moments,"effective dipole moments site")==0 or strcmp(dipole_moments,"effective dipole moments eigen")==0)
//      {
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_2d two_dim_spectra(system, OpenCLinfo, Nmax, dt, DNT, Nsteps3, delayTime, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, padding, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR, Emin, Emax, specify_Enrange);
//         gettimeofday(&start_2d, 0);
//         two_dim_spectra.execute_calculate_2d(dipole_moments, method);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//         two_dim_spectra.print_2d_gnuplot(delayTime, fn_output);
//      }
//      else if(specify_laser_polarization==true)
//      {
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         fprintf(stdout,"#\n# start laser polarization along l=(%e,%e,%e)\n",Pol[0],Pol[1],Pol[2]);
//         char test[500];
//         sprintf(test,"laser polarization along l=(%e,%e,%e)",Pol[0],Pol[1],Pol[2]);
//         calculate_2d two_dim_spectra(system, OpenCLinfo, Nmax, dt, DNT, Nsteps3, delayTime, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, padding, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR, Emin, Emax, specify_Enrange);
//         gettimeofday(&start_2d, 0);
//         if(strcmp(method,"HEOM")==0)
//         {
//            two_dim_spectra.execute_calculate_2d();
//         }
//         if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//         {
//            two_dim_spectra.execute_calculate_Redfield_2d(method);
//         }
//         two_dim_spectra.write_data(test);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//         two_dim_spectra.print_2d_gnuplot(delayTime, fn_output);
//      }
//
//      else if (strcmp(shot,"four shot rotational average")==0)
//      {
//         fprintf(stdout,"#\n# start four shot rotational average\n");
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_2d two_dim_spectra(system, OpenCLinfo, Nmax, dt, DNT, Nsteps3, delayTime, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, padding, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR, Emin, Emax, specify_Enrange);
//         gettimeofday(&start_2d, 0);
//         two_dim_spectra.execute_calculate_2d_rotav_4shots(method);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//         two_dim_spectra.print_2d_gnuplot(delayTime, fn_output);
//      }
//      else if (strcmp(shot,"ten shot rotational average")==0)
//      {
//         fprintf(stdout,"#\n# start ten shot rotational average\n");
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_2d two_dim_spectra(system, OpenCLinfo, Nmax, dt, DNT, Nsteps3, delayTime, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, padding, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR, Emin, Emax, specify_Enrange);
//         gettimeofday(&start_2d, 0);
//         two_dim_spectra.execute_calculate_2d_rotav_10shots(method);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//         two_dim_spectra.print_2d_gnuplot(delayTime, fn_output);
//      }
//      else
//      {
//         fprintf(stdout,"#\n# start shot %s\n",shot);
//         double Pol[3];
//         double phi=(1.+sqrt(5.))/2.;
//         double norm=1./sqrt(3.);
//         if (strcmp(shot,"x only")==0)
//         {
//            Pol[0]=1.0;
//            Pol[1]=0.0;
//            Pol[2]=0.0;
//         }
//         else if (strcmp(shot,"y only")==0)
//         {
//            Pol[0]=0.0;
//            Pol[1]=1.0;
//            Pol[2]=0.0;
//         }
//         else if (strcmp(shot,"z only")==0)
//         {
//            Pol[0]=0.0;
//            Pol[1]=0.0;
//            Pol[2]=1.0;
//         }
//         else if (strcmp(shot,"#1 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=+1.0/sqrt(3.0);
//            Pol[2]=+1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#2 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=+1.0/sqrt(3.0);
//            Pol[2]=-1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#3 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=-1.0/sqrt(3.0);
//            Pol[2]=-1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#4 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=-1.0/sqrt(3.0);
//            Pol[2]=+1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#1 of 10")==0)
//         {
//            Pol[0]=1.*norm;
//            Pol[1]=1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#2 of 10")==0)
//         {
//            Pol[0]=1.*norm;
//            Pol[1]=-1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#3 of 10")==0)
//         {
//            Pol[0]=-1.*norm;
//            Pol[1]=1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#4 of 10")==0)
//         {
//            Pol[0]=-1.*norm;
//            Pol[1]=-1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#5 of 10")==0)
//         {
//            Pol[0]=0;
//            Pol[1]=1./phi*norm;
//            Pol[2]=phi*norm;
//         }
//         else if (strcmp(shot,"#6 of 10")==0)
//         {
//            Pol[0]=0.;
//            Pol[1]=-1./phi*norm;
//            Pol[2]=phi*norm;
//         }
//         else if (strcmp(shot,"#7 of 10")==0)
//         {
//            Pol[0]=1./phi*norm;
//            Pol[1]=phi*norm;
//            Pol[2]=0.;
//         }
//         else if (strcmp(shot,"#8 of 10")==0)
//         {
//            Pol[0]=-1./phi*norm;
//            Pol[1]=phi*norm;
//            Pol[2]=0;
//         }
//         else if (strcmp(shot,"#9 of 10")==0)
//         {
//            Pol[0]=phi*norm;
//            Pol[1]=0;
//            Pol[2]=1./phi*norm;
//         }
//         else if (strcmp(shot,"#10 of 10")==0)
//         {
//            Pol[0]=-phi*norm;
//            Pol[1]=0;
//            Pol[2]=1./phi*norm;
//         }
//         else
//         {
//            fprintf(stderr,"Unknown shot direction %s. ABORT.\n",shot);
//            exit(1);
//         }
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_2d two_dim_spectra(system, OpenCLinfo, Nmax, dt, DNT, Nsteps3, delayTime, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, padding, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR, Emin, Emax, specify_Enrange);
//         gettimeofday(&start_2d, 0);
//         if(strcmp(method,"HEOM")==0)
//         {
//            two_dim_spectra.execute_calculate_2d();
//         }
//         if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//         {
//            two_dim_spectra.execute_calculate_Redfield_2d(method);
//         }
//         two_dim_spectra.write_data(shot);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//         two_dim_spectra.print_2d_gnuplot(delayTime, fn_output);
//      }
//   }

//   else  if(strcmp(calculation_mode,"population dynamics")==0)
//   {
//      System system;
//      system.initialize(
//         Nsites,            // total number of sites in Hamiltonian
//         Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//         const_invcmtomeV,  // convert cm^-1 to meV
//         const_kb,          // bolzmann constant in meV/K
//         const_hbar,  	 // hbar in meV
//         Temp, 		 // Temperature in K
//         lambda_sum,        // lambda_sum at each site already in meV
//         HamMat      	 // Hamiltonian already in meV
//      );
//      calculate_pop population_dynamics(system, OpenCLinfo, Nmax, dt, DNT, Nsteps, siteinit, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output, return_densityMatrix);
//      if(specify_general_initial_state==false)
//      {
//         if(strcmp(method,"HEOM")==0)
//         {
//            population_dynamics.execute_calculate_pop(specify_general_initial_state);
//         }
//         if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//         {
//            population_dynamics.execute_Redfield_calculate_pop(specify_general_initial_state, method);
//         }
//      }
//      // if different initial condition:
//      if(specify_general_initial_state==true)
//      {
//         if(strcmp(method,"HEOM")==0)
//         {
//            population_dynamics.execute_calculate_pop(specify_general_initial_state, valrho_init);
//         }
//         if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//         {
//            population_dynamics.execute_Redfield_calculate_pop(specify_general_initial_state, method, valrho_init);
//         }
//      }
//      compute_time=population_dynamics.diff;
//   }
//
//

//   else if(strcmp(calculation_mode,"pump-probe witness")==0)
//   {
//      int Nsteps3=int(Nsteps/DNT+0.5)*DNT+1; // enforce that Nsteps3 is a multiple of DNT
//      if (Nsteps3!=Nsteps)
//      {
//         fprintf(stdout,"# enforce Nsteps3-1 to fit with multiples of DNT=%d\n",DNT);
//         fprintf(stdout,"# set Nsteps3=%d (instead of %d)\n",Nsteps3,Nsteps);
//      }
//      if(strcmp(dipole_moments,"effective dipole moments site")==0 or strcmp(dipole_moments,"effective dipole moments eigen")==0)
//      {
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_PumpProbeWitness PumpProbeWitness(system, OpenCLinfo, Nmax, dt, Nsteps3, DNT, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks,fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR);
//         gettimeofday(&start_2d, 0);
//         PumpProbeWitness.execute_calculate_PumpProbeWitness(dipole_moments, method);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//      }
//      else if(specify_laser_polarization==true)
//      {
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         char test[500];
//         sprintf(test,"laser polarization along l=(%e,%e,%e)",Pol[0],Pol[1],Pol[2]);
//         calculate_PumpProbeWitness PumpProbeWitness(system, OpenCLinfo, Nmax, dt, Nsteps3, DNT, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks,fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR);
//         gettimeofday(&start_2d, 0);
//         if(strcmp(method,"HEOM")==0)
//         {
//            PumpProbeWitness.execute_calculate_PumpProbeWitness();
//         }
//         if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//         {
//            PumpProbeWitness.execute_calculate_Redfield_PumpProbeWitness(method);
//         }
//         PumpProbeWitness.write_data(test);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//      }
//
//      else if (strcmp(shot,"four shot rotational average")==0)
//      {
//         fprintf(stdout,"#\n# start four shot rotational average\n");
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_PumpProbeWitness PumpProbeWitness(system, OpenCLinfo, Nmax, dt, Nsteps3, DNT, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR);
//         gettimeofday(&start_2d, 0);
//         PumpProbeWitness.execute_calculate_PumpProbeWitness_rotav_4shots(method);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//      }
//      else if (strcmp(shot,"ten shot rotational average")==0)
//      {
//         fprintf(stdout,"#\n# start ten shot rotational average\n");
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_PumpProbeWitness PumpProbeWitness(system, OpenCLinfo, Nmax, dt, Nsteps3, DNT, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR);
//         gettimeofday(&start_2d, 0);
//         PumpProbeWitness.execute_calculate_PumpProbeWitness_rotav_10shots(method);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//      }
//      else
//      {
//         fprintf(stdout,"#\n# run %s\n",shot);
//         double Pol[3];
//         double phi=(1.+sqrt(5.))/2.;
//         double norm=1./sqrt(3.);
//         if (strcmp(shot,"x only")==0)
//         {
//            Pol[0]=1.0;
//            Pol[1]=0.0;
//            Pol[2]=0.0;
//         }
//         else if (strcmp(shot,"y only")==0)
//         {
//            Pol[0]=0.0;
//            Pol[1]=1.0;
//            Pol[2]=0.0;
//         }
//         else if (strcmp(shot,"z only")==0)
//         {
//            Pol[0]=0.0;
//            Pol[1]=0.0;
//            Pol[2]=1.0;
//         }
//         else if (strcmp(shot,"#1 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=+1.0/sqrt(3.0);
//            Pol[2]=+1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#2 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=+1.0/sqrt(3.0);
//            Pol[2]=-1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#3 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=-1.0/sqrt(3.0);
//            Pol[2]=-1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#4 of 4")==0)
//         {
//            Pol[0]=+1.0/sqrt(3.0);
//            Pol[1]=-1.0/sqrt(3.0);
//            Pol[2]=+1.0/sqrt(3.0);
//         }
//         else if (strcmp(shot,"#1 of 10")==0)
//         {
//            Pol[0]=1.*norm;
//            Pol[1]=1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#2 of 10")==0)
//         {
//            Pol[0]=1.*norm;
//            Pol[1]=-1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#3 of 10")==0)
//         {
//            Pol[0]=-1.*norm;
//            Pol[1]=1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#4 of 10")==0)
//         {
//            Pol[0]=-1.*norm;
//            Pol[1]=-1.*norm;
//            Pol[2]=1.*norm;
//         }
//         else if (strcmp(shot,"#5 of 10")==0)
//         {
//            Pol[0]=0;
//            Pol[1]=1./phi*norm;
//            Pol[2]=phi*norm;
//         }
//         else if (strcmp(shot,"#6 of 10")==0)
//         {
//            Pol[0]=0.;
//            Pol[1]=-1./phi*norm;
//            Pol[2]=phi*norm;
//         }
//         else if (strcmp(shot,"#7 of 10")==0)
//         {
//            Pol[0]=1./phi*norm;
//            Pol[1]=phi*norm;
//            Pol[2]=0.;
//         }
//         else if (strcmp(shot,"#8 of 10")==0)
//         {
//            Pol[0]=-1./phi*norm;
//            Pol[1]=phi*norm;
//            Pol[2]=0;
//         }
//         else if (strcmp(shot,"#9 of 10")==0)
//         {
//            Pol[0]=phi*norm;
//            Pol[1]=0;
//            Pol[2]=1./phi*norm;
//         }
//         else if (strcmp(shot,"#10 of 10")==0)
//         {
//            Pol[0]=-phi*norm;
//            Pol[1]=0;
//            Pol[2]=1./phi*norm;
//         }
//         else
//         {
//            fprintf(stderr,"Unknown shot direction %s. ABORT.\n",shot);
//            exit(1);
//         }
//         System system;
//         system.initialize_spectra(
//            Nsites,            // total number of sites in Hamiltonian
//            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//            const_invcmtomeV,  // convert cm^-1 to meV
//            const_kb,          // bolzmann constant in meV/K
//            const_hbar,  	    // hbar in meV
//            Temp, 		    // Temperature in K
//            lambda_sum,        // lambda_sum at each site already in meV
//            HamMat,      	    // Hamiltonian already in meV
//            valDipoleMoments,  // 3-vector of dipole moment at each site
//            Pol                // 3-vector or laser polarization
//         );
//         calculate_PumpProbeWitness PumpProbeWitness(system, OpenCLinfo, Nmax, dt, Nsteps3, DNT, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output, GBRP, GBNR, SERP, SENR, ESARP, ESANR);
//         gettimeofday(&start_2d, 0);
//         if(strcmp(method,"HEOM")==0)
//         {
//            PumpProbeWitness.execute_calculate_PumpProbeWitness();
//         }
//         if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//         {
//            PumpProbeWitness.execute_calculate_Redfield_PumpProbeWitness(method);
//         }
//         PumpProbeWitness.write_data(shot);
//         gettimeofday(&end_2d, 0);
//         compute_time = (1000000.0*(end_2d.tv_sec-start_2d.tv_sec) + end_2d.tv_usec-start_2d.tv_usec)/1000000.0;
//      }
//   }
//

//   else if(strcmp(calculation_mode,"exciton coherence")==0)
//   {
//      System system;
//      system.initialize(
//         Nsites,            // total number of sites in Hamiltonian
//         Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//         const_invcmtomeV,  // convert cm^-1 to meV
//         const_kb,          // bolzmann constant in meV/K
//         const_hbar,  	 // hbar in meV
//         Temp, 		 // Temperature in K
//         lambda_sum,        // lambda_sum at each site already in meV
//         HamMat      	 // Hamiltonian already in meV
//      );
//      calculate_coh exciton_coherence(system, OpenCLinfo, Nmax, dt, DNT, Nsteps, Ei, Ej, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output );
//      if(strcmp(method,"HEOM")==0)
//      {
//         exciton_coherence.execute_calculate_coh();
//      }
//      if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
//      {
//         exciton_coherence.execute_calculate_Redfield_coh(method);
//      }
//      compute_time=exciton_coherence.diff;
//   }
//   else if(strcmp(calculation_mode,"exciton eigenstate population")==0)
//   {
//      System system;
//      system.initialize(
//         Nsites,            // total number of sites in Hamiltonian
//         Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
//         const_invcmtomeV,  // convert cm^-1 to meV
//         const_kb,          // bolzmann constant in meV/K
//         const_hbar,  	 // hbar in meV
//         Temp, 		 // Temperature in K
//         lambda_sum,        // lambda_sum at each site already in meV
//         HamMat      	 // Hamiltonian already in meV
//      );
//      // exciton population dynamics
//      if(HEOM==1)
//      {
//         // do calculation with HEOM
//         calculate_popeigen exciton_popeigen(system, OpenCLinfo, Nmax, dt, DNT, Nsteps, Ei, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output );
//         exciton_popeigen.execute_calculate_popeigen();
//         compute_time+=exciton_popeigen.diff;
//      }
//      // do the same calculation with full Redfield
//      //FIXME[needs to be implemented]
//
//      if(secular_Redfield==1)
//      {
//         // do the same calculation with secular Redfield
//         calculate_popeigen_secular_redfield exciton_popeigen_secular_redfield(system, dt, DNT, Nsteps, Ei, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output );
//         exciton_popeigen_secular_redfield.execute_calculate_popeigen_secular_redfield();
//         compute_time+=exciton_popeigen_secular_redfield.diff;
//      }
//
//      if(generalized_Foerster==1)
//      {
//         // do the same calculation with generalized Foerster
//         calculate_popeigen_general_foerster exciton_popeigen_general_foerster(system, dt, DNT, Nsteps, Ei, tmax, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output );
//         exciton_popeigen_general_foerster.execute_calculate_popeigen_general_foerster();
//         //FIXME TIMIG
//      }
//      if(modified_Redfield==1)
//      {
//         // do the same calculation with modified Redfield
//         calculate_popeigen_modified_redfield exciton_popeigen_modified_redfield(system, dt, DNT, Nsteps, Ei, tmax, NumIter, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output );
//         exciton_popeigen_modified_redfield.execute_calculate_popeigen_modified_redfield();
//         //FIXME TIMIG
//      }
//      if(modified_Redfield_generalized_Foerster==1)
//      {
//         // do the same calcualtion with combinedn mdodified Redfield/generalized Foerster
//         calculate_popeigen_combined exciton_popeigen_combined(system, dt, DNT, Nsteps, Ei, tmax, NumIter,Jcutoff, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, fn_output );
//         exciton_popeigen_combined.execute_calculate_popeigen_combined();
//         //FIXME TIMIG
//      }

//        compute_time=exciton_coherence.diff;

//   }
   if(strcmp(calculation_mode,"transfer efficiency")==0)
   {
      System system;
      if(efficiency_non_hermitian==true)
      {
         system.initialize(
            Nsites,            // total number of sites in Hamiltonian
            Ncoupled,          //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
            const_invcmtomeV,  // convert cm^-1 to meV
            const_kb,          // bolzmann constant in meV/K
            const_hbar,  	 // hbar in meV
            Temp, 		 // Temperature in K
            lambda_sum,        // lambda_sum at each site already in meV
            HamMat      	 // Hamiltonian already in meV
         );
         calculate_efficiency_antiHermitian transfer_efficiency_antiHerm(system, OpenCLinfo, Nmax, dt, DNT, Nsteps, siteinit, sd_lambda_invcm, sd_gamma_fs, sd_gammapos_invcm, sd_CoupledSite, sd_Npeaks, site_to_target, rate_to_target, site_to_loss,rate_to_loss, norm_threshold, fn_output, flux_analysis, dominant_path, flux_threshold, return_densityMatrix);
         if(specify_general_initial_state==false)
         {
            if(strcmp(method,"HEOM")==0)
            {
               transfer_efficiency_antiHerm.execute_calculate_efficiency_antiHermitian(specify_general_initial_state);
            }
            if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
            {
               transfer_efficiency_antiHerm.execute_Redfield_calculate_efficiency_antiHermitian(specify_general_initial_state, method);
            }
         }
         // if different initial condition:
         if(specify_general_initial_state==true)
         {
            if(strcmp(method,"HEOM")==0)
            {
               transfer_efficiency_antiHerm.execute_calculate_efficiency_antiHermitian(specify_general_initial_state, valrho_init);
            }
            if(strcmp(method,"full Redfield")==0 or strcmp(method,"secular Redfield")==0)
            {
               transfer_efficiency_antiHerm.execute_Redfield_calculate_efficiency_antiHermitian(specify_general_initial_state, method, valrho_init);
            }
         }
         compute_time=transfer_efficiency_antiHerm.diff;
      }
   }
   else
   {
      fprintf(stderr,"Unknown calculation mode \"%s\". ABORT.",calculation_mode);
      exit(1);
   }
   OpenCLinfo.clean_up();
// FIXME TIMINGS
   gettimeofday(&t1, 0);
   diff = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
   fprintf(stdout,"#\n# total processing time: %4.2f seconds.\n", diff );
   fprintf(stdout,"# OpenCL-Device computation time: %4.2f seconds.\n", compute_time );
   fprintf(stdout,"# processing time for initializing and output: %4.2f seconds.\n", diff-compute_time );
   time_t t;
   time(&t);
   fprintf(stdout,"# run: successful on %s",ctime(&t));
   print_cite_as();
   exit(0);
}

