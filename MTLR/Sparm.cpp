#include "Sparm.h"
#include <cstdlib>
#include <iostream>

#include <cassert>
#include <fstream>
#include <stdio.h>

Sparm::Sparm(void)
:m_inputFilename(""),
m_modelFilename(""),
m_maxMonth(60),
m_sizePsi(0),
m_original_sizePsi(0),
m_C1(1),
m_C2(1),
m_lossType("l1"),
m_printProb(false),
m_threshold(30), 
m_bundleSize(20), 
m_trainUncensored(1),
m_initialWeightFile(""),
m_quantTime(std::vector<double>()),
m_regularizer(2),
m_xqueries(std::vector<double>()),
m_yqueries(std::vector<double>())
{
}


Sparm::~Sparm(void)
{
}

void Sparm::ReadParam(int argc, char ** argv)
{
  for(int i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) {
    case 'c': i++; m_C1=atof(argv[i]); break;                   // regularization constant C1
    case 'd': i++; m_C2=atof(argv[i]); break;                   // regularization constant C2
	case 'i': i++; m_inputFilename=std::string(argv[i]); break; // input filename
	case 'o': i++; m_modelFilename=std::string(argv[i]); break; // output model filename
    case 'm': i++; m_maxMonth=atoi(argv[i]); break;             // number of time points
	case 'l': i++; m_lossType=std::string(argv[i]); break;      // loss function: available types are: 'l1', 'l2', 'l1-log', 'l2-log', 'rae'
	case 'p': m_printProb=true; break;                          // output string of survival probability for each month
	case 't': i++; m_threshold=atof(argv[i]); break; 
	case 'b': i++; m_bundleSize = atoi(argv[i]); break;                   // bundle size for L-BFGS
	case 'u': i++; m_trainUncensored=atoi(argv[i]); break;                // treat all as 'uncensored' 
	case 'w': i++; m_initialWeightFile=std::string(argv[i]); break;       // initialization weight file for EM-training
	case 'q': i++; ReadQuantTime(std::string(argv[i])); break; // mar 16
	case 'r': i++; m_regularizer = atoi(argv[i]); break; 
	case 'x': i++; m_xqueries.push_back(atof(argv[i])); break; 
	case 'y': i++; m_yqueries.push_back(atof(argv[i])); break; 
	default: std::cerr << "\nUnrecognized option " << argv[i] << "!" << std::endl;
      exit(0);
    }

  }
  // default for m_quantTime; 1 up to m_maxMonth
  // may 9: change default q-file
  /*
  if (m_quantTime.size()==0)
  {
	for (int i=0; i<m_maxMonth; i++)
	{
		m_quantTime.push_back(i+1);
	}
  }
  */

}




void Sparm::ReadQuantTime(std::string filename) 
{
	std::ifstream infile(filename.c_str());
	if (infile) {
		std::string line; 

		while(getline(infile,line))
		{
			double v = atof(line.c_str());
			m_quantTime.push_back(v); 
		}
		infile.close();
		//assert(m_quantTime.size()==m_maxMonth); 
		if ((int) m_quantTime.size()!=m_maxMonth)
		{
			std::cerr << "Number of time intervals in " << filename << " doe not match number of time intervals specified (=" << m_maxMonth << ")!" << std::endl; 
			exit(1); 
		} 
	} else {
		std::cerr << "Unable to open file " << filename << " for input!" << std::endl; 
		exit(1); 
	}
}

