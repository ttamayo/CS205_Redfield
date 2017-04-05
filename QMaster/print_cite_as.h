/*
 * print_cite_as.h
 *
 *  Created on: Feb 13, 2014
 *      Author: christoph
 */

#ifndef PRINT_CITE_AS_H_
#define PRINT_CITE_AS_H_

#include <string>
#include <iostream>
using namespace std;

void print_cite_as()
{
   cout<<"#\n# in all publications reporting results obtained with QMaster cite:"<<endl;
   cout<<"# [1] C. Kreisbeck, T. Kramer, A. Aspuru-Guzik, in preparation (2014)"<<endl;
   cout<<"# [2] C. Kreisbeck and T. Kramer, J. Phys. Chem. Lett., 3, 2828 (2012)"<<endl;
   cout<<"# [3] C. Kreisbeck, T. Kramer, M. Rodriguez and B. Hein, J. Chem. Theory Comput. 7, 2166 (2011) "<<endl;
   cout.flush();
}


#endif /* PRINT_CITE_AS_H_ */
