/*
 * Run_Efficiency_aniHermitian.h
 *
 *  Created on: Dec 18, 2013
 *      Author: christoph
 */

#ifndef RUN_EFFICIENCY_ANTIHERMITIAN_H_
#define RUN_EFFICIENCY_ANTIHERMITIAN_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include "headers.h"
#include <sys/time.h>

class calculate_efficiency_antiHermitian
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Nmax;
   double dt;
   int DNT;
   int Nsteps;
   int siteinit;
   double** sd_lambda_invcm;
   double** sd_gamma_fs;
   double** sd_gammapos_invcm;
   vector<int> sd_CoupledSite;
   vector<int> sd_Npeaks;
   vector<int> site_to_target;
   vector<double> rate_to_target;
   vector<int> site_to_loss;
   vector<double> rate_to_loss;
   double norm_threshold;
   const char *fn_output;
   bool flux_analysis;
   bool dominant_path;
   double flux_threshold;
   bool return_densityMatrix; // if specified to return all entries of the desity matrix in site-basis
   double diff;
private:
   FILE *fo;
   struct timeval t0, t1;
   double* lambda;
   double_complex* gamma;
   double Ediffmax_invcm;
   SigmaHelpMatrices *sigma_Host;
   int Ntuples;
private:
   double GridIntegrate(double dx,vector<double> fGrid);
   double return_efficiency( double dt, vector<double_complex> rho_diag_site, int Ndata);
   double return_average_trapping_time( double dt, vector<double_complex> rho_diag_site, int Ndata, double efficiency);
public:
   calculate_efficiency_antiHermitian(System& INsystem,
                                      OpenCL_Init& INOpenCLinfo,
                                      int INNmax,
                                      double INdt,
                                      int INDNT,
                                      int INNsteps,
                                      int INsiteinit,
                                      double** INsd_lambda_invcm,
                                      double** INsd_gamma_fs,
                                      double** INsd_gammapos_invcm,
                                      vector<int> INsd_CoupledSite,
                                      vector<int> INsd_Npeaks,
                                      vector<int> INsite_to_target,
                                      vector<double> INrate_to_target,
                                      vector<int> INsite_to_loss,
                                      vector<double> INrate_to_loss,
                                      double INnorm_threshold,
                                      const char* INfn_output,
                                      bool INflux_analysis,
                                      bool INdominant_path,
                                      double INflux_threshold,
                                      bool INreturn_densityMatrix);

   void execute_calculate_efficiency_antiHermitian(bool general_rho_init, double_complex *rho_init);
   void execute_Redfield_calculate_efficiency_antiHermitian(bool general_rho_init, const char* method, double_complex *rho_init);
};

calculate_efficiency_antiHermitian::calculate_efficiency_antiHermitian(
   System& INsystem,
   OpenCL_Init& INOpenCLinfo,
   int INNmax,
   double INdt,
   int INDNT,
   int INNsteps,
   int INsiteinit,
   double** INsd_lambda_invcm,
   double** INsd_gamma_fs,
   double** INsd_gammapos_invcm,
   vector<int> INsd_CoupledSite,
   vector<int> INsd_Npeaks,
   vector<int> INsite_to_target,
   vector<double> INrate_to_target,
   vector<int> INsite_to_loss,
   vector<double> INrate_to_loss,
   double INnorm_threshold,
   const char* INfn_output,
   bool INflux_analysis,
   bool INdominant_path,
   double INflux_threshold,
   bool INreturn_densityMatrix)  :
   system(&INsystem),
   OpenCLinfo(&INOpenCLinfo),
   Nmax(INNmax),
   dt(INdt),
   DNT(INDNT),
   Nsteps(INNsteps),
   siteinit(INsiteinit),
   sd_lambda_invcm(INsd_lambda_invcm),
   sd_gamma_fs(INsd_gamma_fs),
   sd_gammapos_invcm(INsd_gammapos_invcm),
   sd_CoupledSite(INsd_CoupledSite),
   sd_Npeaks(INsd_Npeaks),
   site_to_target(INsite_to_target),
   rate_to_target(INrate_to_target),
   site_to_loss(INsite_to_loss),
   rate_to_loss(INrate_to_loss),
   norm_threshold(INnorm_threshold),
   fn_output(INfn_output),
   flux_analysis(INflux_analysis),
   dominant_path(INdominant_path),
   flux_threshold(INflux_threshold),
   return_densityMatrix(INreturn_densityMatrix)
{
   fprintf(stdout, "#\n# initialize parameter for transfer efficiency\n");
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

double calculate_efficiency_antiHermitian::GridIntegrate(double dx,vector<double> fGrid)
{
   //Assume aequidistant Grid with discretization dx
   //number of grid points
   int N=fGrid.size();
   //Integrate with trapez rulee
   double Integral;
   Integral=(fGrid[0]+fGrid[N-1])/2.;
   for(int i=0; i<N-1; i++)
   {
      Integral=Integral+fGrid[i+1];
   }
   Integral=dx*Integral;
   return Integral;
}

double calculate_efficiency_antiHermitian::return_efficiency( double dt, vector<double_complex> rho_diag_site, int Ndata)
{
   double efficiency=0.;
   for(int j=0; j<site_to_target.size(); j++)
   {
      vector<double> Re_part; // only realpart relevant
      for(int it=0; it<Ndata; it++)
      {
         Re_part.push_back(real(rho_diag_site[it*system->Nsites+site_to_target[j]]));
      }
      efficiency+=GridIntegrate(dt, Re_part)/rate_to_target[j];
   }
   return efficiency;
}

double calculate_efficiency_antiHermitian::return_average_trapping_time( double dt, vector<double_complex> rho_diag_site, int Ndata, double efficiency)
{
   double average=0.;
   for(int j=0; j<site_to_target.size(); j++)
   {
      vector<double> Re_part; // only realpart relevant
      for(int it=0; it<Ndata; it++)
      {
         Re_part.push_back(dt*it*real(rho_diag_site[it*system->Nsites+site_to_target[j]]));
      }
      average+=GridIntegrate(dt, Re_part)/rate_to_target[j];
   }
   return average/efficiency;
}


void calculate_efficiency_antiHermitian::execute_calculate_efficiency_antiHermitian(bool general_rho_init, double_complex *rho_init=NULL)
{
   int Nmod=Nsteps/10;  // every 10% progress report to stdout
   fprintf(stdout, "# initialize propagation method=HEOM ...\n");
   sigma_Host=new SigmaHelpMatrices(system->Nsites, Nmax, system->Ncoupled);
   Ntuples=sigma_Host->AllsigmaTuples.size();
   Propagation evolution(*system, *OpenCLinfo, *sigma_Host, Nsteps, dt, Nmax);
   fprintf(stdout, "# ... done");
   LiouvilleHExciton liouHExc(*system,  *OpenCLinfo, Ntuples);
   Liouville_antiHermitian liou_antiHerm(*system,  *OpenCLinfo, Ntuples);
   LiouvillePhononLowTemp liouPhon(*system, *OpenCLinfo,  *sigma_Host, Nmax);
   //
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
         liouPhon.AddSite(sd_CoupledSite[i], sd_Npeaks[i], lambda, gamma);
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
   // only if there is a coupling to the target trap
   if(site_to_target.size()>0)
   {
      fprintf(stdout, "# initialize coupling to the trapping target state ...\n");
      liou_antiHerm.add_sink(site_to_target, rate_to_target);
      fprintf(stdout, "# ... done\n");
   }
   else
   {
      fprintf(stdout, "# warning: no coupling to the trapping target state specified\n");
   }
   // only if there is a coupling to the loss channel
   if(site_to_loss.size()>0)
   {
      // only if there is a coupling
      fprintf(stdout, "# initialize coupling to the loss channel ...\n");
      liou_antiHerm.add_sink(site_to_loss, rate_to_loss);
      fprintf(stdout, "# ... done\n");
   }
   else
   {
      fprintf(stdout, "# warning: no coupling to the loss channel specified\n");
   }
   liou_antiHerm.initialize();
   evolution.Add_Liouville(liou_antiHerm);
   // introduce break condition, if population in loss channel |loss> or target trap |trap>  is lager than 1-norm_threshold
   int breakid=0;  //break id assigns each condition with certain id number, beginning form 0
   BreakByNorm_antiHerm  NormCondition(*system, *OpenCLinfo, DNT*100, breakid);
   evolution.Add_BreakRule(NormCondition, norm_threshold, breakid);
   //
   Matrix rhoInit(system->Nsites); //in Site Basis
   if (general_rho_init==true)
   {
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rhoInit.entry[i]=rho_init[i];
      }
   }
   else
   {
      rhoInit.assignEntry(siteinit,siteinit,1.0); //rhoInit=|siteinit><siteinit|
   }
   evolution.initialize(rhoInit);
   //
   GetDensityMatrixSiteBasisDiag *printrho;
   printrho=new GetDensityMatrixSiteBasisDiag(*system, *OpenCLinfo, DNT, dt);
   evolution.Add_Observable(*printrho);
   returnDensityMatrixSiteBasis *totalrho;
   if(return_densityMatrix==true)
   {
      totalrho=new returnDensityMatrixSiteBasis(*system,*OpenCLinfo,DNT, dt);
      evolution.Add_Observable(*totalrho);
   }
   //
   Population_Flux_antiHerm *population_flux;
   if(flux_analysis==true)
   {
      population_flux=new Population_Flux_antiHerm(*system, *OpenCLinfo, DNT, dt);
      evolution.Add_Observable(*population_flux);
   }
   //
   // RUN PROPAGATION
   fprintf(stdout, "#\n# starting propagation method=HEOM ...\n");
   gettimeofday(&t0, 0);
   evolution.propagate(Nmod);
   gettimeofday(&t1, 0);
   diff = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
   fprintf(stdout, "# ... done\n");
   /*
        *********************************************************
         * Save output values to file
         *********************************************************
      */
   // Number of propagation steps and data points
   int Ndata=printrho->Ndata;
   char fn[500];
   FILE *fo;
   time_t t;
   //
   // output of transfer dynamics
   sprintf(fn,"%s_transfer_dynamics.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time(&t);
   fprintf(fo,"# output generated by z-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: transfer dynamics\n#\n");
   fprintf(fo,"# t in s       ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"site %i       ",j);
   }
   fprintf(fo,"\n");
   for(int it=0; it<Ndata; it++)
   {
      fprintf(fo,"%e ",(double)it*dt*DNT);
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo,"%e ",real(printrho->data[it*system->Nsites+j]));
      }
      fprintf(fo,"\n");
   }
   fclose(fo);
   //
   if(return_densityMatrix==true)
   {
      for(int i=0; i<system->Nsites; i++)
      {
         sprintf(fn,"%s_transfer_dynamics_densityMatrix_rho_%i_k.dat",fn_output, i);
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
         fprintf(out,"# output: transfer dynamics, all entries of the density Matrix of column %i \n#\n",i);
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
   clFlush(OpenCLinfo->queue);
   clFinish(OpenCLinfo->queue);
   cout<<"# free Device-Memory ..."<<endl;
   liouPhon.freeDeviceMemory();
   evolution.freeDeviceMemory();
   liouHExc.freeDeviceMemory();
   liou_antiHerm.freeDeviceMemory();
   cout<<"# ... done"<<endl;
   // Number of propagation steps and data points
   //
   if(flux_analysis==true)
   {
      sprintf(fn,"%s_flux_analysis.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time_t t;
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# output: flux analysis and branching probability\n#\n");
      if(evolution.breakProp==false)
      {
         fprintf(fo,"# warning: break conditions are not yet fulfilled:\n");
         fprintf(fo,"# warning: flux analysis might be inaccurate\n");
         fprintf(fo,"# \n");
      }
      fprintf(fo,"# flux analysis:\n");
      population_flux->compute_flux();
      for(int n=0; n<system->Nsites1; n++)
      {
         for(int m=0; m<system->Nsites1; m++)
         {
            if(m!=n)
            {
               fprintf(fo,"F_%i,%i: %e \n",n,m,imag(population_flux->Fnm->returnEntry(n,m)));
            }
         }
      }
      fprintf(fo,"#\n# branching probability:\n");
      if(evolution.breakProp==false)
      {
         fprintf(fo,"# warning: break conditions are not yet fulfilled:\n");
         fprintf(fo,"# warning: branching probability might be inaccurate\n");
      }
      population_flux->compute_branching_probability(siteinit);
      for(int n=0; n<system->Nsites1; n++)
      {
         for(int m=0; m<system->Nsites1; m++)
         {
            if(m!=n)
            {
               fprintf(fo,"q_%i,%i: %e\n",n,m,imag(population_flux->qnm->returnEntry(n,m)));
            }
         }
      }
      fclose(fo);
      //
      // Prepare for in-flux out-flux diagrams of each site, visualized by pyx, needs to be installed on the system
      // Disable in-flux out-flux diagrams  in the current version, needs pyx (python application) to be installed and adds not much value
//         cout<<"#\n# create in-flux out-flux diagrams ... ";
//         for(int i=0; i<system->Nsites1; i++)
//         {
//            double total_flux_in=0.;
//            double total_flux_out=0;
//            for(int n=0; n<system->Nsites1; n++)
//            {
//               double branch_probability=imag(population_flux->qnm->returnEntry(i,n));
//               if(branch_probability>0)
//               {
//                  total_flux_in+=branch_probability;
//               }
//               if(branch_probability<0)
//               {
//                  total_flux_out+=-branch_probability;
//               }
//            }
//            char parse[5500];
//            sprintf(parse,"%i %i ", system->Nsites1, i);
//            for(int n=0; n<system->Nsites1; n++)
//            {
//               sprintf(parse,"%s %5.3f",parse, imag(population_flux->qnm->returnEntry(i,n)));
//            }
//            sprintf(parse,"%s %i",parse, site_to_target.size());
//            for (int j=0; j< site_to_target.size(); j++)
//            {
//               sprintf(parse,"%s %i",parse, site_to_target[j]);
//            }
//            sprintf(parse,"%s %5.4f",parse, total_flux_in);
//            sprintf(parse,"%s %5.4f",parse, total_flux_out);
//            sprintf(parse,"%s %i",parse, siteinit);
//            sprintf(fn,"%s_InOut",fn_output);
//            sprintf(parse,"%s %s",parse,fn);
//            char submit_command[6000];
//      //      cout<<parse<<endl;
//            sprintf(submit_command,"python ~/QMaster_ExcitonicsCenter/workspace_qmaster-v0.2/qmaster-v0.2/Release/Population_Flux_visualize_InOut_at_Nsite.py %s",parse);
//            std::system(submit_command);
//         }
//         cout<<" done"<<endl;
      //
      // dominant path analyis, cannot guarantee that it works correct for all networks, do not allow this functionality for now
      if(dominant_path==true)
      {
         cout<<"#\n# analyze pathways of population flow:"<<endl;
         if(evolution.breakProp==false)
         {
            cout<<"# warning: break conditions are not yet fulfilled:"<<endl;
            cout<<"# warning: pathways of population flow might be inaccurate"<<endl;
         }
         // show dominant pathways and their relative contribution to the total flux from the initial site to the target sites
         // show pathways if contribution to the total flux is larger than flux_threshold
         for(int i=0; i<site_to_target.size(); i++)
         {
            population_flux->get_dominant_flux_pathways(siteinit, site_to_target[i], flux_threshold);
         }
      }
   }
   // Final ressult efficiency and average trapping time
//         double averageTime;
   double efficiency=return_efficiency(dt*DNT, printrho->data, printrho->Ndata);
   double averageTime=return_average_trapping_time(dt*DNT, printrho->data, printrho->Ndata, efficiency);
   cout<<"#"<<endl;
   cout<<"# transfer characteristics:"<<endl;
   if(evolution.breakProp==false)
   {
      cout<<"# warning: break conditions are not yet fulfilled:"<<endl;
      cout<<"# warning: efficiency and average transfer time might be inaccurate"<<endl;
   }
   cout<<"# transfer efficiency:   "<<efficiency<<endl;
   cout<<"# average transfer time in s: "<<averageTime<<endl;
   cout<<"# "<<endl;


   sprintf(fn,"%s_transfer_characteristics.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time(&t);
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: transfer characteristics\n#\n");
   if(evolution.breakProp==false)
   {
      fprintf(fo,"# warning: break conditions are not yet fulfilled:\n");
      fprintf(fo,"# warning: efficiency and average transfer time might be inaccurate\n");
      fprintf(fo,"# \n");
   }
   fprintf(fo,"transfer efficiency: %e \n",efficiency);
   fprintf(fo,"average transfer time in s: %e \n",averageTime);
   fclose(fo);
}

void calculate_efficiency_antiHermitian::execute_Redfield_calculate_efficiency_antiHermitian(bool general_rho_init, const char* method, double_complex *rho_init=NULL)
{
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
   fprintf(stdout, "# ... done");
   LiouvilleHExcitonRedfield liouHExc(*system,  *OpenCLinfo);
   Liouville_antiHermitianRedfield liou_antiHerm(*system,  *OpenCLinfo, Ntuples);
//   LiouvillePhonon_fullRedfield* liouPhonfull;
   LiouvillePhonon_secularRedfield* liouPhonsecular;
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
//      liouPhonfull=new LiouvillePhonon_fullRedfield(*system, *OpenCLinfo);
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
//            liouPhonfull->AddSite(sd_CoupledSite[i],sd_Npeaks[i], lambda, gamma);
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            liouPhonsecular->AddSite(sd_CoupledSite[i],sd_Npeaks[i], lambda, gamma);
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
//            liouPhonfull->AddSite(sd_CoupledSite[i],sd_Npeaks[i], lambda, gamma);
         }
         if(strncmp(method,"secular Redfield",strlen(method))==0)
         {
            liouPhonsecular->AddSite(sd_CoupledSite[i],sd_Npeaks[i], lambda, gamma);
         }
      }
   }
   fprintf(stdout, "\n# initialize coupling to the phonon-bath ...\n");
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
//      liouPhonfull->initialize();
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      liouPhonsecular->initialize();
   }
   fprintf(stdout, "# ... done\n");
   evolution.Add_Liouville(liouHExc);
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
//      evolution.Add_Liouville(*liouPhonfull);
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      evolution.Add_Liouville(*liouPhonsecular);
   }
   // only if there is a coupling to the target trap
   if(site_to_target.size()>0)
   {
      fprintf(stdout, "# initialize coupling to the trapping target state ...\n");
      liou_antiHerm.add_sink(site_to_target, rate_to_target);
      fprintf(stdout, "# ... done\n");
   }
   else
   {
      fprintf(stdout, "# warning: no coupling to the trapping target state specified\n");
   }
   // only if there is a coupling to the loss channel
   if(site_to_loss.size()>0)
   {
      // only if there is a coupling
      fprintf(stdout, "# initialize coupling to the loss channel ...\n");
      liou_antiHerm.add_sink(site_to_loss, rate_to_loss);
      fprintf(stdout, "# ... done\n");
   }
   else
   {
      fprintf(stdout, "# warning: no coupling to the loss channel specified\n");
   }
   liou_antiHerm.initialize();
   evolution.Add_Liouville(liou_antiHerm);
   // introduce break condition, if population in loss channel |loss> or target trap |trap>  is lager than 1-norm_threshold
   int breakid=0;  //break id assigns each condition with certain id number, beginning form 0
   BreakByNorm_antiHerm  NormCondition(*system, *OpenCLinfo, DNT*100, breakid);
   evolution.Add_BreakRule(NormCondition, norm_threshold, breakid); // Trace is the same in site and eigenbasis. That is I can use the traceoperation defined in BreakbyNorm
   //
   Matrix rhoInit(system->Nsites); //in Site Basis
   if (general_rho_init==true)
   {
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rhoInit.entry[i]=rho_init[i];
      }
   }
   else
   {
      rhoInit.assignEntry(siteinit,siteinit,1.0); //rhoInit=|siteinit><siteinit|
   }
   system->convertEigenbasis(&rhoInit);
   evolution.initialize(rhoInit);
   //
   GetDensityMatrixSiteBasisDiag_fullRedfield *printrho;
   printrho=new GetDensityMatrixSiteBasisDiag_fullRedfield(*system,*OpenCLinfo,DNT, dt);
   evolution.Add_Observable(*printrho);

   returnDensityMatrixSiteBasis_fullRedfield *totalrho;
   if(return_densityMatrix==true)
   {
      totalrho=new returnDensityMatrixSiteBasis_fullRedfield(*system,*OpenCLinfo,DNT, dt);
      evolution.Add_Observable(*totalrho);
   }
   //
   Population_Flux_antiHermRedfield *population_flux; // Note: although Redfield propagates in eigenbasis, the population_flux is evaluated in site-basis
   if(flux_analysis==true)
   {
      population_flux=new Population_Flux_antiHermRedfield(*system, *OpenCLinfo, DNT, dt);
      evolution.Add_Observable(*population_flux);
   }
   //
   // RUN PROPAGATION
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
   gettimeofday(&t1, 0);
   diff = (1000000.0*(t1.tv_sec-t0.tv_sec) + t1.tv_usec-t0.tv_usec)/1000000.0;
   fprintf(stdout, "# ... done\n");
   /*
        *********************************************************
         * Save output values to file
         *********************************************************
      */
   // Number of propagation steps and data points
   int Ndata=printrho->Ndata;
   char fn[500];
   FILE *fo;
   time_t t;
   //
   // output of transfer dynamics
   sprintf(fn,"%s_transfer_dynamics.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time(&t);
   fprintf(fo,"# output generated by z-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: transfer dynamics\n#\n");
   fprintf(fo,"# t in s       ");
   for(int j=0; j<system->Nsites; j++)
   {
      fprintf(fo,"site %i       ",j);
   }
   fprintf(fo,"\n");
   for(int it=0; it<Ndata; it++)
   {
      fprintf(fo,"%e ",(double)it*dt*DNT);
      for(int j=0; j<system->Nsites; j++)
      {
         fprintf(fo,"%e ",real(printrho->data[it*system->Nsites+j]));
      }
      fprintf(fo,"\n");
   }
   fclose(fo);
   //
   if(return_densityMatrix==true)
   {
      for(int i=0; i<system->Nsites; i++)
      {
         sprintf(fn,"%s_transfer_dynamics_densityMatrix_rho_%i_k.dat",fn_output, i);
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
         fprintf(out,"# output: transfer dynamics, all entries of the density Matrix of column %i \n#\n",i);
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
   clFlush(OpenCLinfo->queue);
   clFinish(OpenCLinfo->queue);
   cout<<"# free Device-Memory ..."<<endl;
   if(strncmp(method,"full Redfield",strlen(method))==0)
   {
//      liouPhonfull->freeDeviceMemory();
   }
   if(strncmp(method,"secular Redfield",strlen(method))==0)
   {
      liouPhonsecular->freeDeviceMemory();
   }
   evolution.freeDeviceMemory();
   liouHExc.freeDeviceMemory();
   liou_antiHerm.freeDeviceMemory();
   cout<<"# ... done"<<endl;
   // Number of propagation steps and data points
   //


//
   if(flux_analysis==true)
   {
      sprintf(fn,"%s_flux_analysis.dat",fn_output);
      fo=fopen(fn,"w");
      if (fo==NULL)
      {
         fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
         exit(1);
      }
      time_t t;
      time(&t);
      fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
      fprintf(fo,"# output: flux analysis and branching probability\n#\n");
      if(evolution.breakProp==false)
      {
         fprintf(fo,"# warning: break conditions are not yet fulfilled:\n");
         fprintf(fo,"# warning: flux analysis might be inaccurate\n");
         fprintf(fo,"# \n");
      }
      fprintf(fo,"# flux analysis:\n");
      population_flux->compute_flux();
      for(int n=0; n<system->Nsites1; n++)
      {
         for(int m=0; m<system->Nsites1; m++)
         {
            if(m!=n)
            {
               fprintf(fo,"F_%i,%i: %e \n",n,m,imag(population_flux->Fnm->returnEntry(n,m)));
            }
         }
      }
      fprintf(fo,"#\n# branching probability:\n");
      if(evolution.breakProp==false)
      {
         fprintf(fo,"# warning: break conditions are not yet fulfilled:\n");
         fprintf(fo,"# warning: branching probability might be inaccurate\n");
      }
      population_flux->compute_branching_probability(siteinit);
      for(int n=0; n<system->Nsites1; n++)
      {
         for(int m=0; m<system->Nsites1; m++)
         {
            if(m!=n)
            {
               fprintf(fo,"q_%i,%i: %e\n",n,m,imag(population_flux->qnm->returnEntry(n,m)));
            }
         }
      }
      fclose(fo);
      //
      // Prepare for in-flux out-flux diagrams of each site, visualized by pyx, needs to be installed on the system
      // Disable in-flux out-flux diagrams  in the current version, needs pyx (python application) to be installed and adds not much value
         cout<<"#\n# create in-flux out-flux diagrams ... ";
         for(int i=0; i<system->Nsites1; i++)
         {
            double total_flux_in=0.;
            double total_flux_out=0;
            for(int n=0; n<system->Nsites1; n++)
            {
               double branch_probability=imag(population_flux->qnm->returnEntry(i,n));
               if(branch_probability>0)
               {
                  total_flux_in+=branch_probability;
               }
               if(branch_probability<0)
               {
                  total_flux_out+=-branch_probability;
               }
            }
            char parse[5500];
            sprintf(parse,"%i %i ", system->Nsites1, i);
            for(int n=0; n<system->Nsites1; n++)
            {
               sprintf(parse,"%s %5.3f",parse, imag(population_flux->qnm->returnEntry(i,n)));
            }
            sprintf(parse,"%s %i",parse, site_to_target.size());
            for (int j=0; j< site_to_target.size(); j++)
            {
               sprintf(parse,"%s %i",parse, site_to_target[j]);
            }
            sprintf(parse,"%s %5.4f",parse, total_flux_in);
            sprintf(parse,"%s %5.4f",parse, total_flux_out);
            sprintf(parse,"%s %i",parse, siteinit);
            sprintf(fn,"%s_InOut",fn_output);
            sprintf(parse,"%s %s",parse,fn);
            char submit_command[6000];
      //      cout<<parse<<endl;
            sprintf(submit_command,"python ~/QMaster_ExcitonicsCenter/workspace_qmaster-v0.2/qmaster-v0.2/Release/Population_Flux_visualize_InOut_at_Nsite.py %s",parse);
            std::system(submit_command);
         }
         cout<<" done"<<endl;
      //
//      // dominant path analyis, cannot guarantee that it works correct for all networks, do not allow this functionality for now
      if(dominant_path==true)
      {
         cout<<"#\n# analyze pathways of population flow:"<<endl;
         if(evolution.breakProp==false)
         {
            cout<<"# warning: break conditions are not yet fulfilled:"<<endl;
            cout<<"# warning: pathways of population flow might be inaccurate"<<endl;
         }
         // show dominant pathways and their relative contribution to the total flux from the initial site to the target sites
         // show pathways if contribution to the total flux is larger than flux_threshold
         for(int i=0; i<site_to_target.size(); i++)
         {
            population_flux->get_dominant_flux_pathways(siteinit, site_to_target[i], flux_threshold);
         }
      }
   }
   // Final ressult efficiency and average trapping time
//         double averageTime;
   double efficiency=return_efficiency(dt*DNT, printrho->data, printrho->Ndata);
   double averageTime=return_average_trapping_time(dt*DNT, printrho->data, printrho->Ndata, efficiency);
   cout<<"#"<<endl;
   cout<<"# transfer characteristics:"<<endl;
   if(evolution.breakProp==false)
   {
      cout<<"# warning: break conditions are not yet fulfilled:"<<endl;
      cout<<"# warning: efficiency and average transfer time might be inaccurate"<<endl;
   }
   cout<<"# transfer efficiency:   "<<efficiency<<endl;
   cout<<"# average transfer time in s: "<<averageTime<<endl;
   cout<<"# "<<endl;

   cout<<"# ERROR!!! AntiHERMITION DOES NOT WORK ...\n";
   sprintf(fn,"%s_transfer_characteristics.dat",fn_output);
   fo=fopen(fn,"w");
   if (fo==NULL)
   {
      fprintf(stderr,"Cannot open output parameter file \"%s\". ABORT.",fn);
      exit(1);
   }
   time(&t);
   fprintf(fo,"# output generated by QMaster-%s on %s",INPUT_FORMAT_VERSION,ctime(&t));
   fprintf(fo,"# output: transfer characteristics\n#\n");
   if(evolution.breakProp==false)
   {
      fprintf(fo,"# warning: break conditions are not yet fulfilled:\n");
      fprintf(fo,"# warning: efficiency and average transfer time might be inaccurate\n");
      fprintf(fo,"# \n");
   }
   fprintf(fo,"transfer efficiency: %e \n",efficiency);
   fprintf(fo,"average transfer time in s: %e \n",averageTime);
   fclose(fo);
}



#endif /* GPU_RUNEFFICIENCY_H_ */
