/*
 * Run_PopulationDynamicsEigenbasis_Combined.h
 *
 *  Created on: May 21, 2014
 *      Author: Tobias Kramer
 */

#ifndef RUN_POPULATIONDYNAMICSEIGENBASIS_COMBINED_H_
#define RUN_POPULATIONDYNAMICSEIGENBASIS_COMBINED_H_

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

class calculate_popeigen_combined
{
public:
   System *system;
   double dt;
   int DNT;
   int Nsteps;
   int Ei;
   double tmax;
   int NumIter;
   double M_crit_invcm;
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
   calculate_popeigen_combined(System& INsystem,
                               double INdt,
                               int INDNT,
                               int INNsteps,
                               int INEi,
                               double INtmax,
                               int INNumIter,
                               double INM_crit_invcm,
                               double** INsd_lambda_invcm,
                               double** INsd_gamma_fs,
                               double** INsd_gammapos_invcm,
                               vector<int> INsd_CoupledSite,
                               vector<int> INsd_Npeaks,
                               const char* INfn_output );
   void execute_calculate_popeigen_combined();
};

calculate_popeigen_combined::calculate_popeigen_combined(System& INsystem,
      double INdt,
      int INDNT,
      int INNsteps,
      int INEi,
      double INtmax,
      int INNumIter,
      double INM_crit_invcm,
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
   NumIter(INNumIter),
   M_crit_invcm(INM_crit_invcm),
   sd_lambda_invcm(INsd_lambda_invcm),
   sd_gamma_fs(INsd_gamma_fs),
   sd_gammapos_invcm(INsd_gammapos_invcm),
   sd_CoupledSite(INsd_CoupledSite),
   sd_Npeaks(INsd_Npeaks),
   fn_output(INfn_output)
{
   fprintf(stderr,"# initialize parameter for exciton eigenstate population method=combined modified Redfield/generalized Foerster\n");
   fflush(stdout);
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

void calculate_popeigen_combined::execute_calculate_popeigen_combined()
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

   // compute Rkkp rate tensor (dimension Nsites x Nsites)
   long double *Rkkp;
   long double *Evk,*Evkp;
   long double Ek,Ekp;
   double *Kkkp,*ExpKkkpDt;

   // we need to diagonalize the matrix E_mm+J_mn with |J_mn|>M_crit
   double *JSiteStrong; // J_mn matrix with all entries > M_crit
   double *JSiteWeak;   // J_mn matrix with all entries < M_crit
   double *JEigenWeak;  // JSiteWeak rotated into  eigenbasis
   double *HelpMatrix;
   double *HexSiteStrong; // E_mm+JStrong_mn
   double *eigenvectors;
   double *eigenvalues;

   JSiteStrong  =new double[system->Nsites*system->Nsites];
   JSiteWeak    =new double[system->Nsites*system->Nsites];
   JEigenWeak   =new double[system->Nsites*system->Nsites];
   HexSiteStrong=new double[system->Nsites*system->Nsites];
   HelpMatrix   =new double[system->Nsites*system->Nsites];
   eigenvectors =new double[system->Nsites*system->Nsites];
   eigenvalues  =new double[system->Nsites];


   for(int m=0; m<system->Nsites; m++)
   {
      for(int n=0; n<system->Nsites; n++)
      {
         JSiteStrong  [m*system->Nsites+n]=0.0;
         JSiteWeak    [m*system->Nsites+n]=0.0;
         HelpMatrix   [m*system->Nsites+n]=0.0;
         JEigenWeak   [m*system->Nsites+n]=0.0;
         HexSiteStrong[m*system->Nsites+n]=0.0;
         double Jmn=real(system->Ham->returnEntry(m,n));

         if(m==n)
            HexSiteStrong[m*system->Nsites+n]=Jmn;

         if(fabs(Jmn)>(M_crit_invcm*const_invcmtomeV))
         {
            if(m!=n)
            {
               JSiteStrong[m*system->Nsites+n]=Jmn;
               HexSiteStrong[m*system->Nsites+n]=Jmn;
            }
         }
         else
         {
            if(m!=n)
            {
               JSiteWeak[m*system->Nsites+n]=Jmn;
            }
         }
      }
   }
   // diagonalize the matrix HexSiteStrong and obtain eigenvalues and vectors

   gsl_matrix_view mHexSiteStrong=gsl_matrix_view_array(HexSiteStrong,system->Nsites,system->Nsites);
   gsl_matrix_view meigenvectors=gsl_matrix_view_array(eigenvectors,system->Nsites,system->Nsites);
   gsl_vector_view meigenvalues=gsl_vector_view_array(eigenvalues,system->Nsites);
   gsl_matrix_view mJSiteWeak=gsl_matrix_view_array(JSiteWeak,system->Nsites,system->Nsites);
   gsl_matrix_view mJEigenWeak=gsl_matrix_view_array(JEigenWeak,system->Nsites,system->Nsites);

   gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc(system->Nsites);

   gsl_eigen_symmv(&mHexSiteStrong.matrix,&meigenvalues.vector,&meigenvectors.matrix,w);
   gsl_eigen_symmv_free(w);
   gsl_eigen_symmv_sort(&meigenvalues.vector,&meigenvectors.matrix,GSL_EIGEN_SORT_VAL_ASC); // sort lowest to highest value


   // rotate JSiteWeak into JEigenWeak
   gsl_matrix_view mHelpMatrix = gsl_matrix_view_array(HelpMatrix,system->Nsites,system->Nsites);

   //
//   fprintf(stderr,"meigenvectors.matrix:\n");
//   for (int i = 0; i < system->Nsites; i++)
//   {
//      for (int j = 0; j < system->Nsites; j++)
//      {
////         fprintf(stderr," %+.5f", i, j, gsl_matrix_get(&meigenvectors.matrix, i, j));
//      }
////      fprintf(stderr,"\n");
//   }

//   fprintf(stderr,"mJSiteWeak.matrix:\n");
//   for (int i = 0; i < system->Nsites; i++)
//   {
//      for (int j = 0; j < system->Nsites; j++)
//      {
////         fprintf(stderr," %+.5f", i, j, gsl_matrix_get(&mJSiteWeak.matrix, i, j));
//      }
////      fprintf(stderr,"\n");
//   }

   // eigenvectors are already transposed, JSiteWeak is symmetric
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                   1.0, &mJSiteWeak.matrix, &meigenvectors.matrix,
                   0.0, &mHelpMatrix.matrix);

//   fprintf(stderr,"mHelpMatrix.matrix:\n");
//   for (int i = 0; i < system->Nsites; i++)
//   {
//      for (int j = 0; j < system->Nsites; j++)
//      {
////         fprintf(stderr," %+.5f", i, j, gsl_matrix_get(&mHelpMatrix.matrix, i, j));
//      }
////      fprintf(stderr,"\n");
//   }

   // need to transpose eigenvectors to get the un-transposed matrix as in MM
   gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                   1.0, &meigenvectors.matrix, &mHelpMatrix.matrix,
                   0.0, &mJEigenWeak.matrix);
//   fprintf(stderr,"mJEigenWeak.matrix: [meV]\n");
//   for (int i = 0; i < system->Nsites; i++)
//   {
//      for (int j = 0; j < system->Nsites; j++)
//      {
//         fprintf(stderr," %+.5f", i, j, gsl_matrix_get(&mJEigenWeak.matrix, i, j));
//      }
//      fprintf(stderr,"\n");
//   }
   fprintf(stdout,"# initialize rate-matrix method=combined modified Redfield/generalized Foerster ...\n");
   fflush(stdout);
   Evk =new long double[system->Nsites];
   Evkp=new long double[system->Nsites];
   Rkkp=new long double[system->Nsites*system->Nsites];
   Kkkp     =new double[system->Nsites*system->Nsites];
   ExpKkkpDt=new double[system->Nsites*system->Nsites];
   double tol=1e-7;
   double tol3=1e-9;
   for(int ek=0; ek<system->Nsites; ek++)
   {
      for(int ekp=0; ekp<system->Nsites; ekp++)
      {
         Rkkp[ek*system->Nsites+ekp]=0.0L;

         Ek =(long double)gsl_vector_get(&meigenvalues.vector,ek);
         Ekp=(long double)gsl_vector_get(&meigenvalues.vector,ekp);

//			fprintf(stderr,"Ek [%d]=%f invcm\n",ek ,(double)Ek /const_invcmtomeV);
//			fprintf(stderr,"Ekp[%d]=%f invcm\n",ekp,(double)Ekp/const_invcmtomeV);

         for(int i=0; i<system->Nsites; i++)
         {
            Evk [i]=(long double)gsl_matrix_get(&meigenvectors.matrix,i,ek);
//				fprintf(stderr,"Evk[%d]=%f\n",i,(double)Evk[i]);
         }
         for(int i=0; i<system->Nsites; i++)
         {
            Evkp[i]=(long double)gsl_matrix_get(&meigenvectors.matrix,i,ekp);
//				fprintf(stderr,"Evkp[%d]=%f\n",i,(double)Evkp[i]);
         }

         f_params params;

         params.modified_redfield=true;
         params.ek=ek;
         params.ekp=ekp;

         params.Ek =Ek +gsl_matrix_get(&mJEigenWeak.matrix,ek ,ek );
         params.Ekp=Ekp+gsl_matrix_get(&mJEigenWeak.matrix,ekp,ekp);

         params.Evk =Evk;
         params.Evkp=Evkp;
         params.l_list=l_list;
         params.Temp=system->Temp;
         params.sd_lambda_invcm=sd_lambda_invcm;
         params.sd_gamma_fs=sd_gamma_fs;
         params.sd_gammapos_invcm=sd_gammapos_invcm;
         params.sd_CoupledSite=sd_CoupledSite;
         params.sd_Npeaks=sd_Npeaks;
         params.HCW=gsl_matrix_get(&mJEigenWeak.matrix,ek,ekp);

         params.part=0; // first part for combined method
         gsl_function F;
         F.function = &f_combined;
         F.params = reinterpret_cast<void *>(&params);

         double result1=0.0,error1,result2=0.0,error2,result3=0.0,error3, integral=0.;
         double analytic_result;
         fprintf(stdout,"# initialize rate-matrix [%d,%d]:\n",ek,ekp);
         fflush(stdout);
         gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);

         // first part near singularity up to t=1e-14 s
         gsl_integration_qags(&F,0.0,1.0e-14,0.0,tol3,5000,w,&result1,&error1);
         integral+=result1;
//         printf("result1         = % .18e\n", result1);
//         printf("estimated error1= % .18e\n", error1);
//         printf("intervals =  %d\n", w->size);
         // then the main part t=1e-14..tmax/2.0 s
         gsl_integration_qags(&F, 1.0e-14, tmax/2.0, 0.0, tol, 5000, w, &result2, &error2);
         integral+=result2;
//         printf("result2         = % .18e\n", result2);
//         printf("estimated error2= % .18e\n", error2);
//         printf("intervals =  %d\n", w->size);
         // then the rest t=tmax/2.0..tmax
         // FIXME gsl_integration_qagiu(&F, 1.0e-14, 0.0, 1.0e-5, 1000, w, &result2, &error2);
         // seems not to work
         gsl_integration_qags(&F, tmax/2.0, tmax, 0.0, tol, 5000, w, &result3, &error3);
         integral+=result3;
//         printf("result3         = % .18e\n", result3);
//         printf("estimated error3= % .18e\n", error3);
//         printf("intervals =  %d\n", w->size);

//         printf("result          = % .18e\n", integral);
         int iter=3;
         double told=0., tnew=0., t_total=0.;
         t_total=tmax;
         told=tmax/2.0;
         tnew=tmax;
         while(iter<NumIter)
         {
            iter+=1;
            told=tnew;
            tnew=told+tmax/2.0;
//     	   cout<<"integrate from:"<<told*1e12<<":"<<tnew*1e12<<endl;
//     	   printf("iteration= %i\n", iter);
            result3=0.0;
            gsl_integration_qags(&F,told,tnew,0.0,tol,5000,w,&result3,&error3);
            integral+=result3;
//           printf("result3         = % .18e\n", result3);
//           printf("seg_old         = % .18e\n", seg_old);
//           printf("estimated error3= % .18e\n", error3);
//           printf("result3-seg_old= % .18e\n", abs(result3-seg_old)/seg_old);
//           seg_old=result3;
         }
         gsl_integration_workspace_free(w);
         integral/=const_hbar*const_hbar;
//         if(ek==ekp)
//         {
//        	fprintf(stdout,"# test convergence of numerical integration\n");
//        	fprintf(stdout,"# result[numeric]  = % .18e\n", integral);
////            f(tmax,&params);
//            f((tmax+(NumIter-3)*tmax/2),&params);
//            analytic_result=params.re_gpkkkkp+4.0*params.lkkkk*params.lkkkk*tmax;
//            fprintf(stdout, "# result[analytic] = % .18e\n", analytic_result);
//            fflush(stdout);
////            fprintf(stderr,"lkkkk=%e l_list[0]=%e\n",(double)params.lkkkk,(double)l_list[0]/const_hbar);
//         }

         if(ek!=ekp)
            Rkkp[ek*system->Nsites+ekp]=2.0L*(integral);
         else /* ek==ekp */
         {
            Rkkp[ek*system->Nsites+ekp]=2.0L*(integral);
         }
      }
   }
   delete[] l_list;

   long double result=0.0L;
   for(int c=0; c<system->Nsites; c++)
   {
      if(c==0) continue;
      result+=Rkkp[c*system->Nsites];
//      fprintf(stderr,"%.12f\n",(double)Rkkp[c*system->Nsites]/1.0E11);
   }
//   fprintf(stderr,"result %.12e\n",(double)result);

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
   fprintf(stdout,"# ... done\n");
   fflush(stdout);
   fprintf(stdout,"# Rkkp method=combined modified Redfield/generalized Foerster:\n");
   for(int ek=0; ek<system->Nsites; ek++)
   {
      fprintf(stdout,"# k=%2d",ek);
      for(int ekp=0; ekp<system->Nsites; ekp++)
      {
         fprintf(stdout," %e",(double)Rkkp[ek*system->Nsites+ekp]);
      }
      fprintf(stdout,"\n");
   }
   fprintf(stdout,"# Kkkp matrix method=combined modified Redfield/generalized Foerster:\n");
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
   fprintf(stdout, "#\n# starting propagation method=combined modified Redfield/generalized Foerster ...\n");
   gettimeofday(&t0, 0);

   char fn[500];
   FILE *fo;
   time_t t;

   time(&t);
   sprintf(fn,"%s_exciton_population_combined_modified_Redfield_generalized_Foerster.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in eigenbasis, propagation method=combined modified Redfield/generalized Foerster\n#\n");
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

   //

   // also output in site basis for easier comparison with other program parts
   double *Msite;
   Msite=new double[system->Nsites*system->Nsites];
   double *Menergy;
   Menergy=new double[system->Nsites*system->Nsites];
   gsl_matrix_view msite=gsl_matrix_view_array(Msite,system->Nsites,system->Nsites);
   gsl_matrix_view menergy=gsl_matrix_view_array(Menergy,system->Nsites,system->Nsites);

   time(&t);
   sprintf(fn,"%s_exciton_population_combined_modified_Redfield_generalized_Foerster_site.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   FILE *frho;
   sprintf(fn,"%s_exciton_population_combined_modified_Redfield_generalized_Foerster_rho_init_site.dat",fn_output);
   frho=fopen(fn,"w");
   if (frho==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }


   fprintf(frho,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(frho,"# output: initial density matrix in site basis, format as QMaster input, propagation method=combined modified Redfield/generalized Foerster\n#\n");
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: exciton population in sitebasis, propagation method=combined modified Redfield/generalized Foerster\n#\n");
   fprintf(fo,"# t in s       ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"site %i       ",j);
   }
   fprintf(fo,"\n");
   fprintf(fo,"%e",dt*0);

   // initial matrix
   gsl_matrix_set_all(&m.matrix,0.0);
   gsl_matrix_set(&m.matrix,Ei,Ei,1.0);

   // transform energy -> site basis
   gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                   1.0, &m.matrix, &meigenvectors.matrix,
                   0.0, &mHelpMatrix.matrix);

   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                   1.0, &meigenvectors.matrix, &mHelpMatrix.matrix,
                   0.0, &msite.matrix);
   for(int i=0; i<system->Nsites; i++)
   {
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(frho,"rho_init:%d:%d=%e:0. \n",i,j,gsl_matrix_get(&msite.matrix,i,j));
      }
   }
   fclose(frho);

   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo," %e",gsl_matrix_get(&msite.matrix,j,j));
   }
   fprintf(fo,"\n");

   // m=em;
   gsl_matrix_memcpy(/* dest */&m.matrix,/* source */ &em.matrix);
   for(int i=1; i<Nsteps; i++)
   {
      fprintf(fo,"%e",dt*i);
      // take the row corresponding to the scalar product with the initial eigenstate
      gsl_matrix_get_row(&v.vector,&m.matrix,Ei);
      // construct basis in energy

      // initial matrix
      gsl_matrix_set_all(&menergy.matrix,0.0);

      for(int j=0; j<system->Nsites; j++)
      {
         gsl_matrix_set(&menergy.matrix,j,j,gsl_vector_get(&v.vector,j));
      }

      // transform energy -> site basis
      gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                      1.0, &menergy.matrix, &meigenvectors.matrix,
                      0.0, &mHelpMatrix.matrix);

      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                      1.0, &meigenvectors.matrix, &mHelpMatrix.matrix,
                      0.0, &msite.matrix);

      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo," %e",gsl_matrix_get(&msite.matrix,j,j));
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

#endif /* RUN_POPULATIONDYNAMICSEIGENBASIS_COMBINED_H_ */
