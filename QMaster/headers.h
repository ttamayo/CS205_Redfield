#ifndef HEADERS_H
#define HEADERS_H

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

typedef double myreal;
typedef cl_double2 cl_myreal2;

#include "Kernels/KernelSource.h"
#include "DefineQMaster-Version.h"
#include "print_cite_as.h"
#include "Constants_hbarandCo.h"
#include "opencl_error.h"
#include "OpenCL_initialize.h"
#include "PlotSpectralDensity.h"
#include "SigmaHelpMatrices.h"
#include "Tuple.h"
#include "Matrix.h"
#include "System.h"
#include "Liouville.h"
#include "Observable.h"
#include "BreakRule.h"
#include "BreakByNorm_antiHerm.h"
#include "returnTrace_mumRho.h"
#include "returnTrace_mumRho_Redfield.h"
#include "Kernel_InitSigma.h"
#include "RungeKutta.h"
#include "Propagation.h"
#include "GetDensityMatrixSiteBasisDiag.h"
#include "GetDensityMatrixSiteBasisDiag_fullRedfield.h"
#include "GetDensityMatrixElement_EigenBasis.h"
#include "GetDensityMatrixElement_EigenBasis_Redfield.h"
#include "returnDensityMatrixSiteBasis.h"
#include "returnDensityMatrixSiteBasis_fullRedfield.h"
#include "Population_Flux_antiHermitian.h"
#include "Population_Flux_antiHermitian_Redfield.h"
#include "Kernel_HPhononLowTemp.h"
#include "LiouvillePhononLowTemp.h"
#include "Kernel_HPhononDoubleExcitonLowTemp.h"
#include "LiouvillePhononDoubleExcitonLowTemp.h"
//#include "LiouvillePhonon_fullRedfield.h"
//#include "LiouvillePhononDoubleExciton_fullRedfield.h"
#include "LiouvillePhonon_secularRedfield.h"
#include "Kernel_Hexciton.h"
#include "LiouvilleHExciton.h"
#include "Kernel_HexcitonRedfield.h"
#include "LiouvilleHExcitonRedfield.h"
#include "Kernel_antiHermitian.h"
#include "Liouville_antiHermitian.h"
#include "Kernel_antiHermitianRedfield.h"
#include "Liouville_antiHermitianRedfield.h"
#include "LiouvilleHPulse.h"
//#include "LiouvilleHPulseRedfield.h"
#include "Kernel_multiply_dipole_sigmas.h"
#include "multiply_dipole_sigmas.h"
#include "multiply_dipole_sigmas_Redfield.h"
#include "Kernel_return_witness.h"
#include "return_witness.h"
//#include "return_witness_Redfield.h"
//#include "Run_PopulationDynamics.h"
//#include "Run_Coherences.h"
// //#include "Run_1dspectra.h"
//#include "Run_2dspectra.h"
#include "Run_Efficiency_antiHermitian.h"
//#include "Run_PumpProbeWitness.h"
//#include "Run_applyPulse.h"
//#include "Run_PopulationDynamicsEigenbasis.h"
//#include "Run_PopulationDynamicsEigenbasis_SecularRedfield.h"
//#include "Run_PopulationDynamicsEigenbasis_GeneralFoerster.h"
//#include "Run_PopulationDynamicsEigenbasis_ModifiedRedfield.h"
//#include "Run_PopulationDynamicsEigenbasis_Combined.h"


#endif



