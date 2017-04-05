#ifndef FILEDATA_H
#define FILEDATA_H

#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;
// Einlesen von Spalten aus einer Datei
class Filedata
{
public:
   const char *INfilename;
   int INNCol;
   vector<double> ColVals;
   int Nlines;
private:
   ifstream datei;

public:
   Filedata(const char *filename,int NCol);
   ~Filedata();
   void getcolumn();
   void closeFile(void);
};

Filedata::Filedata(const char *filename, int NCol):
   INfilename(filename), INNCol(NCol)
{
   // Teste ob Datei existiert
   datei.open(INfilename, ios_base::in);
   if(datei.good())
   {
   }
   else
   {
      cout<<"# Datei:"<<endl;
      cout<<"# "<<INfilename<<endl;
      cout<<"# existiert nicht"<<endl;
   }
}

Filedata::~Filedata()
{
   datei.close();
}
void Filedata::closeFile(void)
{
   datei.close();
}

void Filedata::getcolumn()
{
   int breakcolumn=INNCol;
   bool breakvalue=false;
   int colactual=0;
   int end=0;
   vector<int> endlist;
   vector<int> start;
   start.push_back(0);
// #Skip Header
   char buffer[500];
   char s;
   char strdummy='#';
   while(strdummy == '#')
   {
//     getline(datei,buffer,'\n');
      datei.getline(buffer,500);
      strdummy=buffer[0];
   }
// ## lese ersten Datenpunkt ein
   if(! datei.eof())
   {
      while(colactual != INNCol and breakvalue==false)
      {
         colactual=colactual+1;
         s='s';
         while(s!= ' ')
         {
            if(end<strlen(buffer))
            {
               s=buffer[end];
               end=end+1;
            }
            else
            {
               s=' ';
               breakvalue=true;
               breakcolumn=colactual;
            }
         }
         start.push_back(end);
         endlist.push_back(end);
      }
   }
   if(breakcolumn<INNCol)
   {
      cout << "ACHTUNG FEHLER:";
      cout << "In File " << INfilename <<" maximal " << breakcolumn << " Spalten vorhanden" << endl;
   }
   else
   {
      unsigned int k=endlist[INNCol-1]-start[INNCol-1];
      string str;
      for(int i=start[INNCol-1]; i<endlist[INNCol-1]; i++)
      {
         str=str+buffer[i];
      }
      double test2=atof(str.c_str());
      ColVals.push_back(test2);
   }

   // # lese restliche Datenpunkte ein
   //
   datei.getline(buffer,500);
   while(! datei.eof())
   {
      colactual=0;
      end=0;
      breakvalue=false;
      start.erase(start.begin(),start.end());  // loesche Alles in vector start und endlist
      start.push_back(0);
      endlist.erase(endlist.begin(),endlist.end());
      while(colactual != INNCol and breakvalue==false)
      {
         colactual=colactual+1;
         s='s';
         while(s!= ' ')
         {
            if(end<strlen(buffer))
            {
               s=buffer[end];
               end=end+1;
            }
            else
            {
               s=' ';
               breakvalue=true;
               breakcolumn=colactual;
            }
         }
         start.push_back(end);
         endlist.push_back(end);
      }
      if(breakcolumn<INNCol)
      {
      }
      else
      {
         unsigned int k=endlist[INNCol-1]-start[INNCol-1];
         string str;
         for(int i=start[INNCol-1]; i<endlist[INNCol-1]; i++)
         {
            str=str+buffer[i];
         }
         double test2=atof(str.c_str());
         ColVals.push_back(test2);
      }
      datei.getline(buffer,500);
   }
   Nlines=ColVals.size();
}

#endif
