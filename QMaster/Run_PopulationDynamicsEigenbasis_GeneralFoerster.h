/*
 * Run_PopulationDynamicsEigenbasis_GeneralFoerster.h
 *
 *  Created on: May 21, 2014
 *      Author: Tobias Kramer
 */

#ifndef RUN_POPULATIONDYNAMICSEIGENBASIS_GENERALFOERSTER_H_
#define RUN_POPULATIONDYNAMICSEIGENBASIS_GENERALFOERSTER_H_

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
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_integration.h>
#include "hyp2f1.h"

class calculate_popeigen_general_foerster
{
public:
   System *system;
   double dt;
   int DNT;
   int Nsteps;
   int Ei;
   double tmax;
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
   calculate_popeigen_general_foerster(System& INsystem,
                                       double INdt,
                                       int INDNT,
                                       int INNsteps,
                                       int INEi,
                                       double INtmax,
                                       double** INsd_lambda_invcm,
                                       double** INsd_gamma_fs,
                                       double** INsd_gammapos_invcm,
                                       vector<int> INsd_CoupledSite,
                                       vector<int> INsd_Npeaks,
                                       const char* INfn_output );
   void execute_calculate_popeigen_general_foerster();
   void print_f_function();
};

calculate_popeigen_general_foerster::calculate_popeigen_general_foerster(System& INsystem,
      double INdt,
      int INDNT,
      int INNsteps,
      int INEi,
      double INtmax,
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
   tmax(INtmax),
   sd_lambda_invcm(INsd_lambda_invcm),
   sd_gamma_fs(INsd_gamma_fs),
   sd_gammapos_invcm(INsd_gammapos_invcm),
   sd_CoupledSite(INsd_CoupledSite),
   sd_Npeaks(INsd_Npeaks),
   fn_output(INfn_output)
{
   fprintf(stdout, "#\n# initialize parameter for exciton eigenstate population method=generalized Foerster\n");
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

void calculate_popeigen_general_foerster::execute_calculate_popeigen_general_foerster()
{
   double *l_list;
   l_list=new double[system->Nsites];
   for(int i=0; i<system->Nsites; i++)
   {
      l_list[i]=0.0;
   }

   for(int i=0; i<sd_CoupledSite.size(); i++)
   {
      double lambda_sum=0.0;
      //handle original DL spectral density
      if(sd_Npeaks[i]==1 and sd_gammapos_invcm[i][0]==0.0)
      {
         lambda=new double[sd_Npeaks[i]];
         gamma=new double_complex[sd_Npeaks[i]];
         for(int k=0; k<sd_Npeaks[i]; k++)
         {
            lambda[k]=sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(sd_gamma_fs[i][k]*1.0e-15), 0.);
            lambda_sum+=lambda[k];
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
            lambda_sum+=(lambda[2*k+0]+lambda[2*k+1]);
         }
         plot_spectral_density(sd_CoupledSite[i], 2*sd_Npeaks[i], lambda, gamma, fn_output, Ediffmax_invcm);
         delete[] lambda;
         delete[] gamma;
      }
      // need to know total lambda at each SITE
      int siteindex=sd_CoupledSite[i];
      l_list[siteindex]=lambda_sum;
   }
   fprintf(stdout,"# initialize rate-matrix method=generalized Foerster ...\n");
   fflush(stdout);
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
   double tol=1e-7;
   for(int ek=0; ek<system->Nsites; ek++)
   {
      for(int ekp=0; ekp<system->Nsites; ekp++)
      {
         Rkkp[ek*system->Nsites+ekp]=0.0L;
         Ek =system->Eigenvals[ek ];
         Ekp=system->Eigenvals[ekp];
         //fprintf(stderr,"Ek [%d]=%f invcm\n",ek ,(double)Ek /const_invcmtomeV);
         //fprintf(stderr,"Ekp[%d]=%f invcm\n",ekp,(double)Ekp/const_invcmtomeV);
         for(int i=0; i<system->Nsites; i++)
         {
            Evk [i]=(long double)real(system->U->returnEntry(ek,i));
            //				fprintf(stderr,"Evk[%d]=%e\n",i,(double)Evk[i]);
         }
         for(int i=0; i<system->Nsites; i++)
         {
            Evkp[i]=(long double)real(system->U->returnEntry(ekp,i));
            //				fprintf(stderr,"Evkp[%d]=%e\n",i,(double)Evkp[i]);
         }
         long double Vkkp;
         Vkkp=0.0L;
         for(int n=0; n<system->Nsites; n++)
         {
            for(int m=0; m<system->Nsites; m++)
            {
               //fprintf(stderr,"H[%d,%d]=%e\n",n,m,real(system->Ham->returnEntry(n,m)));
               if(m!=n) Vkkp+=Evk[n]*Evkp[m]*(long double)real(system->Ham->returnEntry(n,m)/const_hbar);
            }
         }
         Vkkp=powl(Vkkp,2);
         f_params params;
         params.modified_redfield=false;
         params.ek=ek;
         params.ekp=ekp;
         params.Ek =Ek;
         params.Ekp=Ekp;
         params.Evk =Evk;
         params.Evkp=Evkp;
         params.l_list=l_list;
         params.Temp=system->Temp;
         params.sd_lambda_invcm=sd_lambda_invcm;
         params.sd_gamma_fs=sd_gamma_fs;
         params.sd_gammapos_invcm=sd_gammapos_invcm;
         params.sd_CoupledSite=sd_CoupledSite;
         params.sd_Npeaks=sd_Npeaks;
         gsl_function F;
         F.function = &f;
         F.params = reinterpret_cast<void *>(&params);

         //Print f-funtion
//FILE *gfu;
//gfu=fopen("gfu.dat","w");
//         for(int i=0; i<2000; i++)
//         {
//          double time=i*1.e-15;
//          fprintf(gfu, "%e %e\n", time, f(time, &params));
//         }
//fclose(gfu);

         // general foerster is not singular near t=0, thus can take 0..tmax/2
         double result1=0.0,error1,result2=0.0,error2,integral=0;
         fprintf(stdout,"# initialize rate-matrix [%d,%d]:\n",ek,ekp);
         fflush(stdout);
         gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
         double tol2=1e-14;
         gsl_integration_qags(&F,0.0,tmax/2.0,0.0,tol,5000,w,&result1,&error1);
//           printf("result1         = % .18e\n", result1);
//           printf("estimated error1= % .18e\n", error1);
//           printf("estimated relatvive_error1= % .18e\n", error1/result1);
//           printf("intervals =  %d\n", w->size);
         integral+=result1;
         gsl_integration_qags(&F,tmax/2.0,tmax,0.0,tol,5000,w,&result2,&error2);
//           printf("result2         = % .18e\n", result2);
//           printf("estimated error2= % .18e\n", error2);
//           printf("result2/integral= % .18e\n", result2/integral);
//           printf("estimated relatvive_error2= % .18e\n", error2/result2);
//           printf("intervals =  %d\n", w->size);
         integral+=result2;
         int iter=2;
         int maxIter=100000;
         double told, tnew, t_total=0.;
         t_total=tmax;
         told=tmax/2.0;
         tnew=tmax;
         while(abs(result2/integral)>tol2)
         {
            if(iter>maxIter)
            {
               fprintf(stderr,"Error in create Foerster-Rates: Integration does not converge!\nWe recommend to use a different type of spectral density! ABORT.\n");
               exit(1);
            }
            iter+=1;
            told=tnew;
            tnew=told+tmax/2.0;
            result2=0;
            // cout<<"integrate from:"<<told*1e12<<":"<<tnew*1e12<<endl;
            // printf("iteration= %i\n", iter);
            gsl_integration_qags(&F,told,tnew,0.0,tol,5000,w,&result2,&error2);
            integral+=result2;
//        	              printf("result2         = % .18e\n", result2);
//        	              printf("estimated error2= % .18e\n", error2);
//        	              printf("integral        = % .18e\n", integral);
//        	              printf("result2/integral= % .18e\n", result2/integral);
//        	              printf("estimated relatvive_error2= % .18e\n", error2/result2);
//        	              printf("intervals =  %d\n", w->size);
            t_total=tnew;
         }
         // FIXME gsl_integration_qagiu(&F, 1.0e-14, 0.0, 1.0e-5, 1000, w, &result3, &error3);
         // seems not to work
         gsl_integration_workspace_free(w);
         Rkkp[ek*system->Nsites+ekp]=2.0L*Vkkp*(integral);
      }
   }
   fprintf(stdout,"# ... done\n");
   fflush(stdout);
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
   fprintf(stdout,"# Rkkp matrix method=generalized Foerster:\n");
   for(int ek=0; ek<system->Nsites; ek++)
   {
      fprintf(stdout,"# k=%2d",ek);
      for(int ekp=0; ekp<system->Nsites; ekp++)
      {
         fprintf(stdout," %e",(double)Rkkp[ek*system->Nsites+ekp]);
      }
      fprintf(stdout,"\n");
   }
   fflush(stdout);
   fprintf(stdout,"# Kkkp matrix method=generalized Foerster:\n");
   for(int ek=0; ek<system->Nsites; ek++)
   {
      fprintf(stdout,"# k=%2d",ek);
      for(int ekp=0; ekp<system->Nsites; ekp++)
      {
         fprintf(stdout," %e",(double)Kkkp[ek*system->Nsites+ekp]);
      }
      fprintf(stdout,"\n");
   }
   fflush(stdout);
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
   fprintf(stdout, "#\n# starting propagation method=generalized Foerster ...\n");
   gettimeofday(&t0, 0);

   char fn[500];
   FILE *fo;
   time_t t;
   time(&t);
   sprintf(fn,"%s_exciton_population_generalized_Foerster.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in eigenbasis, propagation method=generalized Foerster\n#\n");
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
   // m=em;
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
   gettimeofday(&t1, 0);
   diff = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
   fprintf(stdout, "# ... done\n");

   // also output in site basis for easier comparison with other program parts
   Matrix rho(system->Nsites); //in Site Basis
   Matrix rhoEigbasis(system->Nsites);
   Matrix *Help1;
   Help1=new Matrix(system->Nsites1);

   time(&t);
   sprintf(fn,"%s_exciton_population_generalized_Foerster_site.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in sitebasis, propagation method=generalized Foerster\n#\n");
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

   delete[] population;
   delete[] Mnew;
   delete[] Evk;
   delete[] Evkp;
   delete[] Rkkp;
   delete[] Kkkp;
   delete[] ExpKkkpDt;
}

#endif /* RUN_POPULATIONDYNAMICSEIGENBASIS_GENERALFOERSTER_H_ */
