/*
 * Run_PopulationDynamicsEigenbasis.h
 *
 *  Created on: May 21, 2014
 *      Author: Tobias Kramer
 */

#ifndef RUN_POPULATIONDYNAMICSEIGENBASIS_SECULARREDFIELD_H_
#define RUN_POPULATIONDYNAMICSEIGENBASIS_SECULARREDFIELD_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "headers.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

class calculate_popeigen_secular_redfield
{
public:
   System *system;
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
   bool PrintAllPopulationEnergyBasis;
public:
   calculate_popeigen_secular_redfield(System& INsystem,
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
   void execute_calculate_popeigen_secular_redfield();
};


calculate_popeigen_secular_redfield::calculate_popeigen_secular_redfield(System& INsystem,
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
   fprintf(stdout, "#\n# initialize parameter for exciton eigenstate population method=secular Redfield\n");
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

void calculate_popeigen_secular_redfield::execute_calculate_popeigen_secular_redfield()
{
   int Nmod=Nsteps/10;  // every 10% progress report to stdout
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
         delete[] lambda;
         delete[] gamma;
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
         delete[] lambda;
         delete[] gamma;
      }
   }
   fprintf(stdout, "# initialize Rate-Tensor method=secular Redfield ...\n");
   // compute Rkkp rate tensor (dimension Nsites x Nsites)
   long double *Rkkp;
   long double *Evk,*Evkp;
   long double Ek,Ekp;
   double *Kkkp,*ExpKkkpDt;

   Evk =new long double[system->Nsites];
   Evkp=new long double[system->Nsites];
   Rkkp=new long double[system->Nsites*system->Nsites];

   Kkkp     =new double[system->Nsites*system->Nsites];
   ExpKkkpDt=new double[system->Nsites*system->Nsites];
   for(int ek=0; ek<system->Nsites; ek++)
   {
      for(int ekp=0; ekp<system->Nsites; ekp++)
      {
         Rkkp[ek*system->Nsites+ekp]=0.0L;

         Ek =system->Eigenvals[ek ];
         Ekp=system->Eigenvals[ekp];

//			fprintf(stderr,"Ek [%d]=%f invcm\n",ek ,(double)Ek /const_invcmtomeV);
//			fprintf(stderr,"Ekp[%d]=%f invcm\n",ekp,(double)Ekp/const_invcmtomeV);

         for(int i=0; i<system->Nsites; i++)
         {
            Evk [i]=(long double)real(system->U->returnEntry(ek ,i));
//				fprintf(stderr,"Evk[%d]=%e\n",i,(double)Evk[i]);
         }
         for(int i=0; i<system->Nsites; i++)
         {
            Evkp[i]=(long double)real(system->U->returnEntry(ekp,i));
//				fprintf(stderr,"Evkp[%d]=%e\n",i,(double)Evkp[i]);
         }

         // this sum runs only over the couples sites
         for(int i=0; i<sd_CoupledSite.size(); i++)
         {
            // this is the real site index 0..Nsites (required to look up the correct eigenvector element)
            int siteindex=sd_CoupledSite[i];
            if(sd_Npeaks[i]==1 and sd_gammapos_invcm[i][0]==0.0)
            {
               lambda=new double[sd_Npeaks[i]];
               gamma=new double_complex[sd_Npeaks[i]];
               for(int k=0; k<sd_Npeaks[i]; k++)
               {
                  lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
                  gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15), 0.);
               }

               long double DeltaFreq=(Ekp-Ek)/const_hbar;
               for(int k=0; k<sd_Npeaks[i]; k++)
               {
                  if(ek==ekp)
                  {
                     Rkkp[ek*system->Nsites+ekp]+=
                        2.0L*const_kb*system->Temp*real(gamma[k])*lambda[k]/(powl(real(gamma[k]),2.0))/const_hbar/const_hbar
                        *2.0L*powl(Evk[siteindex]*Evkp[siteindex],2);
                  }
                  else
                  {
                     Rkkp[ek*system->Nsites+ekp]+=
                        2.0L*real(gamma[k])*lambda[k]*DeltaFreq/(powl(real(gamma[k]),2.0)+powl(DeltaFreq,2.0))/const_hbar
                        *2.0L*powl(Evk[siteindex]*Evkp[siteindex],2)
                        *(1.0L/(expl((Ekp-Ek)/(long double)(const_kb*system->Temp))-1.0L)+1.0L);
                  }
               }
               delete[] gamma;
               delete[] lambda;
            }
            else
            {
               long double DeltaFreq=(Ekp-Ek)/const_hbar;
               for(int k=0; k<sd_Npeaks[i]; k++)
               {
                  lambda=new double[2*sd_Npeaks[i]];
                  gamma=new double_complex[2*sd_Npeaks[i]];

                  lambda[2*k]=sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
                  lambda[2*k+1]=sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
                  gamma[2*k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),+sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
                  gamma[2*k+1]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15),-sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);

                  if(ek==ekp)
                  {
                     Rkkp[ek*system->Nsites+ekp]+=
                        2.0L*const_kb*system->Temp*real(gamma[2*k+0])*lambda[2*k+0]/(powl(real(gamma[2*k+0]),2.0)+powl(imag(gamma[2*k+0]),2.0))/const_hbar/const_hbar
                        *2.0L*powl(Evk[siteindex]*Evkp[siteindex],2);
                     Rkkp[ek*system->Nsites+ekp]+=
                        2.0L*const_kb*system->Temp*real(gamma[2*k+1])*lambda[2*k+1]/(powl(real(gamma[2*k+1]),2.0)+powl(imag(gamma[2*k+1]),2.0))/const_hbar/const_hbar
                        *2.0L*powl(Evk[siteindex]*Evkp[siteindex],2);
                  }
                  else
                  {
                     Rkkp[ek*system->Nsites+ekp]+=
                        2.0L*real(gamma[2*k+0])*lambda[2*k+0]*DeltaFreq/(powl(real(gamma[2*k+0]),2.0)+powl(DeltaFreq+imag(gamma[2*k+0]),2.0))/const_hbar
                        *2.0L*powl(Evk[siteindex]*Evkp[siteindex],2)
                        *(1.0L/(expl((Ekp-Ek)/(long double)(const_kb*system->Temp))-1.0L)+1.0L);
                     Rkkp[ek*system->Nsites+ekp]+=
                        2.0L*real(gamma[2*k+1])*lambda[2*k+1]*DeltaFreq/(powl(real(gamma[2*k+1]),2.0)+powl(DeltaFreq+imag(gamma[2*k+1]),2.0))/const_hbar
                        *2.0L*powl(Evk[siteindex]*Evkp[siteindex],2)
                        *(1.0L/(expl((Ekp-Ek)/(long double)(const_kb*system->Temp))-1.0L)+1.0L);
                  }
                  delete[] lambda;
                  delete[] gamma;
               }
            }
         }
      }
   }

   for(int a=0; a<system->Nsites; a++)
   {
      for(int b=0; b<system->Nsites; b++)
      {
         Kkkp[a*system->Nsites+b]=0.0L;
         if (a!=b)
            Kkkp[a*system->Nsites+b]=Rkkp[b*system->Nsites+a];
         else // a==b
         {
            for(int c=0; c<system->Nsites; c++)
            {
               if(a!=c)
                  Kkkp[a*system->Nsites+b]-=Rkkp[c*system->Nsites+a];
            }
         }
      }
   }
   fprintf(stdout, "# ... done\n");
//   fprintf(stderr,"Rkkp matrix [secular Redfield]\n");
//   for(int ek=0; ek<system->Nsites; ek++)
//   {
//      fprintf(stderr,"k=%2d",ek);
//      for(int ekp=0; ekp<system->Nsites; ekp++)
//      {
//         fprintf(stderr," %e",(double)Rkkp[ek*system->Nsites+ekp]);
//      }
//      fprintf(stderr,"\n");
//   }
//
//   fprintf(stderr,"Kkkp matrix [secular Redfield]\n");
//   for(int ek=0; ek<system->Nsites; ek++)
//   {
//      fprintf(stderr,"k=%2d",ek);
//      for(int ekp=0; ekp<system->Nsites; ekp++)
//      {
//         fprintf(stderr," %e",(double)Kkkp[ek*system->Nsites+ekp]);
//      }
//      fprintf(stderr,"\n");
//   }

   // multiply Kkkp with dt
   for(int i=0; i<system->Nsites*system->Nsites; i++)
   {
      Kkkp[i]*=dt;
   }
   // exponentiate exp[K*dt]
   gsl_matrix_view m =gsl_matrix_view_array(Kkkp,     system->Nsites,system->Nsites);
   gsl_matrix_view em=gsl_matrix_view_array(ExpKkkpDt,system->Nsites,system->Nsites);
   gsl_linalg_exponential_ss(&m.matrix, &em.matrix,0);

   double *population;
   population=new double[system->Nsites];
   gsl_vector_view v =gsl_vector_view_array(population,system->Nsites);

   double *Mnew;
   Mnew=new double[system->Nsites*system->Nsites];
   gsl_matrix_view mnew=gsl_matrix_view_array(Mnew,system->Nsites,system->Nsites);

   // RUN PROPAGATION
   fprintf(stdout, "#\n# starting propagation method=secular Redfield ...\n");
   gettimeofday(&t0, 0);

   char fn[500];
   FILE *fo;
   time_t t;

   time(&t);
   sprintf(fn,"%s_exciton_population_secular_Redfield.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in eigenbasis, propagation method=secular Redfield\n#\n");
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
   gsl_matrix_memcpy(/* dest */&m.matrix,/* source */ &em.matrix);
   for(int i=1; i<Nsteps; i++)
   {
      fprintf(fo,"%e",dt*i);
      // take the row corresponding to the scalar product with the initial eigenstate
      gsl_matrix_get_row(&v.vector,&m.matrix,Ei);
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo," %e",gsl_vector_get(&v.vector,j));
      }
      fprintf(fo,"\n");
      // mnew=em.m;
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0, &m.matrix, &em.matrix, 0.0, &mnew.matrix);
      // m=mnew;
      gsl_matrix_memcpy(/* dest */&m.matrix,/* source */ &mnew.matrix);
   }
   fclose(fo);
   // also output in site basis for easier comparison with other program parts
   Matrix rho(system->Nsites); //in Site Basis
   Matrix rhoEigbasis(system->Nsites);
   Matrix *Help1;
   Help1=new Matrix(system->Nsites1);

   time(&t);
   sprintf(fn,"%s_exciton_population_secular_Redfield_site.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in sitebasis, propagation method=secular Redfield\n#\n");
   fprintf(fo,"# t in s       ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"site %i       ",j);
   }
   fprintf(fo,"\n");
   fprintf(fo,"%e",dt*0);

   //Transform rho_EigBasis to Sitebasis
   rhoEigbasis.assignEntry(Ei, Ei, 1.);
   MatrixMuliplication(Help1, &rhoEigbasis, system->U, system->Nsites);
   MatrixMuliplication(&rho,system->UT,Help1, system->Nsites);
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo," %e",rho.returnEntry(j,j)); // initial population
   }
   fprintf(fo,"\n");
   // m=em;
   gsl_matrix_memcpy(/* dest */&m.matrix,/* source */ &em.matrix);
   for(int i=1; i<Nsteps; i++)
   {
      fprintf(fo,"%e",dt*i);
      // take the row corresponding to the scalar product with the initial eigenstate
      gsl_matrix_get_row(&v.vector,&m.matrix,Ei);
      for(int j=0; j<system->Nsites; j++)
      {
         rhoEigbasis.assignEntry(j,j,gsl_vector_get(&v.vector,j));
      }
      MatrixMuliplication(Help1, &rhoEigbasis, system->U, system->Nsites);
      MatrixMuliplication(&rho,system->UT,Help1, system->Nsites);
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo," %e",rho.returnEntry(j,j)); // initial population
      }
      fprintf(fo,"\n");
      // mnew=em.m;
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0, &m.matrix, &em.matrix, 0.0, &mnew.matrix);
      // m=mnew;
      gsl_matrix_memcpy(/* dest */&m.matrix,/* source */ &mnew.matrix);
   }
   fclose(fo);
   gettimeofday(&t1, 0);
   diff = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
   fprintf(stdout, "# ... done\n");

   delete[] population;
   delete[] Mnew;

   delete[] Evk;
   delete[] Evkp;
   delete[] Rkkp;
   delete[] Kkkp;
   delete[] ExpKkkpDt;
}

#endif /* RUN_POPULATIONDYNAMICSEIGENBASIS_SECULARREDFIELD_H_ */
