/*
 * Population_Flux_antiHermitian_Redfield.h
 *
 *  Created on: May 03, 2014
 *      Author: christoph
 */

#ifndef POPULATION_FLUX_ANTIHERMITIAN_REDFIELD_H_
#define POPULATION_FLUX_ANTIHERMITIAN_REDFIELD_H_

#include "headers.h"
using namespace::std;

class Population_Flux_antiHermRedfield: public Observable
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int t_mod;
   double dt;
   Matrix* Fnm; //Matrix of flux components
   Matrix* qnm; //Matrix for branching probability
private:
   Matrix *rho;
   double *h_rhoEntries;
   vector< vector<int> > path;
   vector< vector<int> > path_old;
   vector<double> flux_per_path;
   vector<double> flux_per_path_old;
public:
   Population_Flux_antiHermRedfield(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, double INdt); // computes the time derivative of rho_00.
   void execute(double time, int it, cl_mem d_sigma_System);
   void compute_flux();
   void compute_branching_probability(int siteinit);
   void get_dominant_flux_pathways(int siteinit, int target, double threshold);
};



Population_Flux_antiHermRedfield::Population_Flux_antiHermRedfield(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, double INdt) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), t_mod(INt_mod), dt(INdt)
{
   h_rhoEntries=new double[2*system->Nsites*system->Nsites];
   rho=new Matrix(system->Nsites);
   Fnm=new Matrix(system->Nsites);
   qnm=new Matrix(system->Nsites);
}


void Population_Flux_antiHermRedfield::compute_flux()
{
   for(int n=0; n<system->Nsites; n++)
   {
      for(int m=0; m<system->Nsites; m++)
      {
         double_complex F;
         F=2./system->hbar*system->Ham->returnEntry(n,m)*Fnm->returnEntry(n, m)*dt*(double)(t_mod); // F_nm=2*J_nm*tau_mn, tau_mn is stored in Fnm(n,m)
         Fnm->assignEntry(n,m,F);
      }
   }
}

void Population_Flux_antiHermRedfield::compute_branching_probability(int siteinit)
{
   for(int n=0; n<system->Nsites; n++)
   {
      for(int m=0; m<system->Nsites; m++)
      {
         double Norm=0;
         for(int nn=0; nn<system->Nsites; nn++)
         {
            if(n==siteinit)
            {
               //no incomming contributions since initial condition. normalize to all outgoing components
               if(imag(Fnm->returnEntry(n,nn))<0)
               {
                  Norm+=abs(imag(Fnm->returnEntry(n,nn)));
               }
            }
            else
            {
               if(imag(Fnm->returnEntry(n,nn))>0)
               {
                  Norm+=imag(Fnm->returnEntry(n,nn));
               }
            }
         }
         double_complex q;
         q=Fnm->returnEntry(n,m)/abs(Norm);
         qnm->assignEntry(n,m,q);
      }
   }
}


void Population_Flux_antiHermRedfield::get_dominant_flux_pathways(int siteinit, int target, double threshold)
{
// initialize path_old and path, flux_per_path, flux_per_path_old
   path_old.erase(path_old.begin(), path_old.end());
   path.erase(path.begin(), path.end());
   flux_per_path_old.erase(flux_per_path_old.begin(), flux_per_path_old.end());
   flux_per_path.erase(flux_per_path.begin(), flux_per_path.end());

   vector<int> start_path;
   start_path.push_back(siteinit);
   path_old.push_back(start_path);
// from site init to network
   for (int node=0; node<system->Nsites; node++)
   {
      double q=imag(qnm->returnEntry(siteinit,node));
      if(q<0 and abs(q)>threshold and node!=siteinit)
      {
         for(int i=0; i<path_old.size(); i++)
         {
            vector<int> flux_flow=path_old[i];
            flux_flow.push_back(node);
            path.push_back(flux_flow);
            flux_per_path.push_back(abs(q));
         }
      }
   }

   for(int k=0; k<system->Nsites-2; k++)
   {
// extend the path by the next site
// step 1: reinitialize path_old and path
      path_old.erase(path_old.begin(), path_old.end());
      for(int i=0; i<path.size(); i++)
      {
         path_old.push_back(path[i]);
      }
      path.erase(path.begin(), path.end());
// reinitialize flux_per_path_old and flux_per_path
      flux_per_path_old.erase(flux_per_path_old.begin(), flux_per_path_old.end());
      for(int i=0; i<flux_per_path.size(); i++)
      {
         flux_per_path_old.push_back(flux_per_path[i]);
      }
      flux_per_path.erase(flux_per_path.begin(), flux_per_path.end());
// step 2: add new site to the path
      for(int i=0; i< path_old.size(); i++)
      {
         // if path already terminates at target, nothing to do keep path
         if(path_old[i].back()==target)
         {
            path.push_back(path_old[i]);
            flux_per_path.push_back(flux_per_path_old[i]);
            continue;
         }
         for(int j=0; j<system->Nsites; j++)
         {
            vector<int> flux_flow =path_old[i];
            // if site j is not yet part of the path
            if(find(flux_flow.begin(), flux_flow.end(), j)==flux_flow.end())
            {
               double q=imag(qnm->returnEntry(flux_flow.back(),j));
               double prob=0;
               if(q<0 and abs(q)>threshold)
               {
                  prob=flux_per_path_old[i]*abs(q);
               }
               if(prob>threshold)
               {
                  flux_flow.push_back(j);
                  path.push_back(flux_flow);
                  flux_per_path.push_back(prob);
               }
            }
         }
      }
   }
   cout<<"# dominant pathways connecting initial site "<<siteinit<<" with the target site "<<target<<endl;
   cout<<"# flux_threshold = "<<threshold<<endl;
// output of the pathways
   for(int i=0; i<path.size(); i++)
   {
      cout<<"# ";
      for(int j=0; j< path[i].size(); j++)
      {
         cout<<path[i][j]<<"->";
      }
      cout<<"relative flux="<<flux_per_path[i]<<endl;
   }
   double total_dominant_target=0.;
   for(int i=0; i<flux_per_path.size(); i++)
   {
      total_dominant_target+=flux_per_path[i];
   }
//   cout<<"all dominant pathways to target "<<target<<": relative contribution to the total flux:"<<total_dominant_target<<endl;
   cout<<"# relative contribution to the total flux of all dominant pathways to target "<<target<<": "<<total_dominant_target<<endl;
}


void Population_Flux_antiHermRedfield::execute(double time, int it, cl_mem d_sigma_System)
{
   // cacluataes on-fly tau_mn=int dt rho_mn  [uses Trapez rule and assumes that rho_mn goes to zero for large times, rho_mn(N_laststep*dt)=0
   if(it%t_mod==0)
   {
      clEnqueueReadBuffer(OpenCLinfo->queue, d_sigma_System, CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries,0, NULL, NULL);
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rho->entry[i]=h_rhoEntries[2*i]+double_complex(0.,h_rhoEntries[2*i+1]);
      }
      system->convertSitebasis(rho);
      if(it==0) // start with 0 at t_0=0
      {
         for(int n=0; n<system->Nsites; n++)
         {
            for(int m=0; m<system->Nsites; m++)
            {
               if(n!=m)
               {
                  // note Fnm goes with tau_mn
                  Fnm->EntryAdd(n, m, 0.5*rho->returnEntry(m,n)); //Trapez-rule -> factor 0.5*
               }
            }
         }
      }
      else
      {
         for(int n=0; n<system->Nsites; n++)
         {
            for(int m=0; m<system->Nsites; m++)
            {
               if(n!=m)
               {
                  // note Fnm goes with tau_mn
                  Fnm->EntryAdd(n, m, rho->returnEntry(m,n));
               }
            }
         }
      }
   }
}






#endif
