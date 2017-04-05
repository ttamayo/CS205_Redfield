/*
 * Run_PopulationDynamicsEigenbasis.h
 *
 *  Created on: May 21, 2014
 *      Author: Tobias Kramer
 */

#ifndef RUN_POPULATIONDYNAMICSEIGENBASIS_H_
#define RUN_POPULATIONDYNAMICSEIGENBASIS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "headers.h"

class calculate_popeigen
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Nmax;
   double dt;
   int DNT;
   int Nsteps;
   int Ei;
   double** sd_lambda_invcm;
   double** sd_gamma_fs;
   double** sd_gammapos_invcm;
   vector<int> sd_CoupledSite;
   vector<int> sd_Npeaks;
   const char *fn_output;
   double diff;
private:
   struct timeval t0, t1;
   double* lambda;
   double_complex* gamma;
   double Ediffmax_invcm;
   SigmaHelpMatrices *sigma_Host;
   int Ntuples;
public:
   calculate_popeigen(System& INsystem,
                      OpenCL_Init& INOpenCLinfo,
                      int INNmax,
                      double INdt,
                      int INDNT,
                      int INNsteps,
                      int INEi,
                      double** INsd_lambda_invcm,
                      double** INsd_gamma_fs,
                      double** INsd_gammapos_invcm,
                      vector<int> INsd_CoupledSite,
                      vector<int> INsd_Npeaks,
                      const char* INfn_output );
   void execute_calculate_popeigen();
};


calculate_popeigen::calculate_popeigen(System& INsystem,
                                       OpenCL_Init& INOpenCLinfo,
                                       int INNmax,
                                       double INdt,
                                       int INDNT,
                                       int INNsteps,
                                       int INEi,
                                       double** INsd_lambda_invcm,
                                       double** INsd_gamma_fs,
                                       double** INsd_gammapos_invcm,
                                       vector<int> INsd_CoupledSite,
                                       vector<int> INsd_Npeaks,
                                       const char* INfn_output )  :
   system(&INsystem),
   OpenCLinfo(&INOpenCLinfo),
   Nmax(INNmax),
   dt(INdt),
   DNT(INDNT),
   Nsteps(INNsteps),
   Ei(INEi),
   sd_lambda_invcm(INsd_lambda_invcm),
   sd_gamma_fs(INsd_gamma_fs),
   sd_gammapos_invcm(INsd_gammapos_invcm),
   sd_CoupledSite(INsd_CoupledSite),
   sd_Npeaks(INsd_Npeaks),
   fn_output(INfn_output)
{
   fprintf(stdout, "#\n# initialize parameter for exciton eigenstate population method=HEOM\n");
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

void calculate_popeigen::execute_calculate_popeigen()
{
   bool PrintAllPopulationEnergyBasis=true;

   int Nmod=Nsteps/10;  // every 10% progress report to stdout
   fprintf(stdout, "# initialize propagation method=HEOM ...\n");
   sigma_Host=new SigmaHelpMatrices(system->Nsites, Nmax, system->Ncoupled);
   Ntuples=sigma_Host->AllsigmaTuples.size();
   Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
   fprintf(stdout, "# ... done");
   LiouvilleHExciton liouHExc(*system, *OpenCLinfo, Ntuples);
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
         liouPhon.AddSite(sd_CoupledSite[i],sd_Npeaks[i], lambda, gamma);
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
         liouPhon.AddSite(sd_CoupledSite[i],2*sd_Npeaks[i], lambda, gamma);
      }
   }
   fprintf(stdout, "\n# initialize coupling to the phonon-bath ...\n");
   liouPhon.initialize();
   fprintf(stdout, "# ... done\n");
   evolution.Add_Liouville(liouHExc);
   evolution.Add_Liouville(liouPhon);
   Matrix rhoInit(system->Nsites); //in Site Basis
   Matrix rhoInitEigbasis(system->Nsites);
   rhoInitEigbasis.assignEntry(Ei, Ei, 1.); //Ei-1, Ek-1 is required since Ei,Ej=1..Nsites
   //Transform rhoInit_EigBasis to Sitebasis
   Matrix *Help1;
   Help1=new Matrix(system->Nsites1);
   MatrixMuliplication(Help1, &rhoInitEigbasis, system->U, system->Nsites1);
   MatrixMuliplication(&rhoInit,system->UT,Help1, system->Nsites1);
   evolution.initialize(rhoInit);
   // create out-file
   GetDensityMatrixElement_EigenBasis *printCoh;
   printCoh=new GetDensityMatrixElement_EigenBasis(*system, *OpenCLinfo, DNT, Ei, Ei, PrintAllPopulationEnergyBasis);
   evolution.Add_Observable(*printCoh);
   // RUN PROPAGATION
   fprintf(stdout, "#\n# starting propagation method=HEOM ...\n");
   gettimeofday(&t0, 0);
   evolution.propagate(Nmod);
   clFlush(OpenCLinfo->queue);
   clFinish(OpenCLinfo->queue);
   gettimeofday(&t1, 0);
   diff = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
   fprintf(stdout, "# ... done\n");
   /*
        *********************************************************
         * Save output values to file
         *********************************************************
        */
   char fn[500];
   FILE *fo;
   time_t t;
   time(&t);

   sprintf(fn,"%s_exciton_population_HEOM.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in eigenbasis, propagation method=HEOM\n#\n");
   fprintf(fo,"# t in s  ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"Re[rho_E%dE%d]  ",j,j);
   }
   fprintf(fo,"\n");
   fprintf(fo,"%e",dt*0);
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo," %e",(j==Ei)?1.0:0.0); // initial population
   }
   fprintf(fo,"\n");
   for(int it=1; it<printCoh->data.size()/(system->Nsites*system->Nsites); it++)
   {
      fprintf(fo,"%e",(double)it*dt*DNT);
      for(int i=0; i<system->Nsites; i++)
      {
         for(int j=0; j<system->Nsites; j++)
         {
            if (i==j) fprintf(fo," %e",real(printCoh->data[it*(system->Nsites*system->Nsites)+i*system->Nsites+j]));
         }
      }
      fprintf(fo,"\n");
   }
   fclose(fo);


   sprintf(fn,"%s_exciton_population_HEOM_site.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in sitebasis, propagation method=HEOM\n#\n");
   fprintf(fo,"# t in s       ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"site %i       ",j);
   }
   fprintf(fo,"\n");
   for(int it=0; it<printCoh->data.size()/(system->Nsites*system->Nsites); it++)
   {
      fprintf(fo,"%e",(double)it*dt*DNT);

      // FIXME: incorrect assumption only diagonal entries important, I need the full density matrix!
      for(int i=0; i<system->Nsites; i++)
      {
         for(int j=0; j<system->Nsites; j++)
         {
            rhoInitEigbasis.assignEntry(i,j,printCoh->data[it*(system->Nsites*system->Nsites)+i*system->Nsites+j]);
         }
      }
      MatrixMuliplication(Help1, &rhoInitEigbasis, system->U, system->Nsites);
      MatrixMuliplication(&rhoInit,system->UT,Help1, system->Nsites);
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo," %e",real(rhoInit.returnEntry(j,j)));
      }
      fprintf(fo,"\n");
   }
   fclose(fo);

   cout<<"# free Device-Memory ..."<<endl;
   liouPhon.freeDeviceMemory();
   evolution.freeDeviceMemory();
   liouHExc.freeDeviceMemory();
   cout<<"# ... done"<<endl;
}

#endif /* RUN_POPULATIONDYNAMICSEIGENBASIS_H_ */
