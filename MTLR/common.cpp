#include "common.h"
#include <cassert>
#include <algorithm>
#include <cmath>
#include <stdio.h>

DenseVector* LoadWeights(std::string filename, Sparm &sparm)
{

	ifstream infile(filename.c_str());
	if (infile) {
		string line; 
		// may 9: new default for q-file
		getline(infile,line);
		vector<string> my_pair1 = Tokenize(line, ":"); 
		size_t numTimePoints = Scan<size_t>(my_pair1[1]); 
		getline(infile,line);
		my_pair1 = Tokenize(line,",");
		vector<double> q_vector;
		for (size_t i=0; i<my_pair1.size(); i++)
		{
			q_vector.push_back(Scan<double>(my_pair1[i])); 
		}
		assert(q_vector.size()==numTimePoints); 
		sparm.SetQuantTime(q_vector); 
		// end may 9
		// aug 2: for reading regularization setting L1/L2
		getline(infile,line); 
		vector<string> my_pair2 = Tokenize(line, ":"); 
		int r = Scan<int>(my_pair2[1]); 
		sparm.SetRegularizer(r); 
		// end aug 2 
		getline(infile,line); 
		vector<string> my_pair = Tokenize(line, ":");
		size_t dim = Scan<size_t>(my_pair[1]);
		
		DenseVector *w = new DenseVector(dim); 

		while(getline(infile,line))
		{
			my_pair = Tokenize(line, ":");
			size_t i = Scan<size_t>(my_pair[0]);
			double v = Scan<double>(my_pair[1]);
			(*w)[i] = v; 
		}
		infile.close();
		return w; 
	} else {
		cerr << "Unable to open model file " << filename << " for input!" << endl; 
		return NULL; 
	}
}



// may 9: construct q vector from input file
void construct_intervals(string input_filename, Sparm &sparm)
{
	std::string line;
	ifstream inputStream(input_filename.c_str());

	if (inputStream.fail())
	{
		//cout << "Cannot read from input file " << input_filename << "!" << endl; 
		cerr << "Cannot read from input file " << input_filename << "!" << endl; 
		exit(1);
	}
	vector<double> event_time; 
	while (!getline(inputStream, line, '\n').eof())
	{
		// process line
		std::string::size_type lastPos = line.find_first_of(" \n",0); 
		double survival_time = 0; 
		survival_time = atof((line.substr(0, lastPos).c_str()));

		event_time.push_back(survival_time); 
	}
	// sort event time
	sort(event_time.begin(), event_time.end());
	vector<double> q_vector;
	// oct 4: avoid duplicates
	double last_event_time = -1;
	for (size_t i=0; i<event_time.size(); i++)
	{
		if ((i % 10)==0) 
		{
			// oct 4: avoid duplicates
			if (fabs(last_event_time-event_time[i])>1E-10)
			{
				q_vector.push_back(event_time[i]); 
			}
			last_event_time = event_time[i]; 
		}
	}
	sparm.SetQuantTime(q_vector); 

	inputStream.close(); 
}



// read input examples
vector<Example> read_input_examples(string input_filename, Sparm &sparm)
{
	vector<Example> sample; 
	std::string line;
	ifstream inputStream(input_filename.c_str());
	
	size_t maxFeatureNum = 0; 

	if (inputStream.fail())
	{
		//cout << "Cannot read from input file " << input_filename << "!" << endl; 
		cerr << "Cannot read from input file " << input_filename << "!" << endl; 
		exit(1);
	}

	const vector<double>& qt_time = sparm.GetQuantTime(); 

	while (!getline(inputStream, line, '\n').eof())
	{
		// process line
		std::string::size_type lastPos = line.find_first_of(" \n",0); 
		double survival_time = 0; 
		survival_time = atof((line.substr(0, lastPos).c_str()));

		// censoring status
		std::string::size_type censoredPos = line.find_first_of(" \n", lastPos+1);
		int censoring_status = atoi(line.substr(lastPos+1,censoredPos-lastPos).c_str()); 
		bool c; 
		if (censoring_status==1) 
		{
			c = true; 
		} else {
			c= false;
		}
		lastPos = censoredPos; 

		vector<pair<size_t,double> > feature_vec; 

		std::string::size_type pos = line.find_first_of(':', lastPos);
		while (std::string::npos != pos || std::string::npos != lastPos)
		{
			size_t i = (size_t) atoi((line.substr(lastPos, pos - lastPos).c_str())); 
			lastPos = line.find_first_of(" \n", pos);
			double v = atof((line.substr(pos+1, lastPos - pos).c_str())); 
			pos = line.find_first_of(':', lastPos);

			if (i>maxFeatureNum) maxFeatureNum = i; 
			feature_vec.push_back(make_pair(i,v)); 
		}

		SparseVector fvec(feature_vec); 

		int n = 0;
		while ((n<sparm.GetMaxMonth())&&(survival_time>qt_time[n]))
		{
			n++; 
		}

		sample.push_back(Example(fvec, n, c, survival_time)); 
	}
	sparm.SetSizePsi(maxFeatureNum+1); 
	inputStream.close(); 

	return(sample); 

}


