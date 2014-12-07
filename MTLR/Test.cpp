

#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "common.h"



vector<double> compute_survival_pdf(const DenseVector &w, const SparseVector &x, const Sparm &sparm)
{
	int maxMonth = sparm.GetMaxMonth();
	double b_constant = 100;
	if (sparm.GetRegularizer()==1) b_constant=1;

	vector<double> scores(maxMonth+1,0); 

	double sum=0;
	double max_score = 0;
	for (int i=1; i<maxMonth+1; i++)
	{
		size_t bias_offset = maxMonth*sparm.GetOriginalSizePsi() + (i-1); 
		size_t offset = (i-1)*sparm.GetOriginalSizePsi(); 
		scores[i] = sprod_ns(w,x, offset) + w[bias_offset]*b_constant;
		// sum stores the scores of the sequence of all 0's
		sum += scores[i];
		if (scores[i]>max_score) max_score = scores[i]; 
	}

	vector<double> probs(maxMonth+1,0); 
	probs[0] = sum; 
	double z=1;


	double running_score=sum; 
	for (int i=1; i<maxMonth+1; i++) 
	{
		running_score -= scores[i]; // flipping signs
		probs[i] = running_score; 
		z += exp(running_score-sum); 
	}
	double log_z = log(z) + sum; 
	for (int i=0; i<maxMonth+1; i++) 
	{
		probs[i] = exp(probs[i] - log_z); 
	}

	return(probs); 
}


bool my_comp(const vector<double> &a, const vector<double> &b)
{
	return (a[0]<b[0]); 
}


// z is vector of 3-tupes (t,p,c)
// assume sorted in t
double concordanceIndex(vector<vector<double> > &z)
{
	size_t n = z.size();
	sort(z.begin(), z.end(), my_comp); 

	int agree = 0;
	int disagree = 0;
	int tie = 0; 
	for (int i=0; i<n; i++)
	{
		if (z[i][2]<1-1E-10) // uncensored
		{
			for (int j=i+1; j<n; j++)
			{
				if (z[j][0]>z[i][0]+1E-10) // only consider non-ties in event time
				{
					if (z[j][1]>z[i][1]+1E-10)
					{
						agree++;
					} else {
						if (z[j][1]<z[i][1]-1E-10)
						{
							disagree++; 
						} else {
							tie++; 
						}
					}
				}
			}
		}
	}

	return ((double) agree + 0.5*tie)/(agree+disagree+tie); 
}


double rae_loss(double prediction, double truth)
{
	double ans; 
	if (prediction>0) {
		ans = fabs(((double) (prediction-truth))/prediction); 
	} else {
		ans = 1; 
	}
	if (ans>1) ans = 1; 

	return(ans);
}


double l2_loss(double prediction, double truth)
{
	return((double) (prediction-truth)*(prediction-truth)); 
}


double l1_loss(double prediction, double truth)
{
	return(fabs((double) (prediction-truth))); 
}


double l2_log_loss(double prediction, double truth)
{
	double log_prediction, log_truth;
	if (prediction==0) {
		// substitute e^-1 for 0; acceptable for survival time prediction
		log_prediction=-1; 
	} else {
		log_prediction=log((double) prediction); 
	}
	if (truth==0) {
		log_truth=-1; 
	} else {
		log_truth=log((double) truth); 
	}
	return((log_prediction-log_truth)*(log_prediction-log_truth)); 
}


double l1_log_loss(double prediction, double truth)
{
	double log_prediction, log_truth;
	if (prediction==0) {
		log_prediction=-1; 
	} else {
		log_prediction=log((double) prediction); 
	}
	if (truth==0) {
		log_truth=-1; 
	} else {
		log_truth=log((double) truth); 
	}

	return(fabs(log_prediction-log_truth)); 
}

/*
double class_loss(double prediction, double truth, double threshold)
{
	if (((prediction>threshold)&&(truth>threshold))||((prediction<=threshold)&&(truth<=threshold)))
	{
		return(0); 
	} else {
		return(1); 
	}
}
*/


DenseVector* remap_weights(DenseVector *w, Sparm &sparm)
{
	const int maxMonth = sparm.GetMaxMonth(); 
	const size_t sizePsi = sparm.GetSizePsi(); 
	const size_t oriSizePsi = sparm.GetOriginalSizePsi(); 

	size_t trainSizePsi = (w->dim()-maxMonth)/maxMonth; 	

	if (trainSizePsi>oriSizePsi)
	{
		sparm.SetSizePsi(trainSizePsi); 
		return w;
	}

	DenseVector *new_w = new DenseVector(sizePsi); 
	for (int i=0; i<maxMonth; i++) 
	{
		size_t train_bias_offset = maxMonth*trainSizePsi + i; 
		size_t train_offset = i*trainSizePsi; 

		size_t bias_offset = maxMonth*oriSizePsi + i; 
		size_t offset = i*oriSizePsi; 

		for (size_t j=0; j<trainSizePsi; j++)
		{
			(*new_w)[offset+j] = (*w)[train_offset+j]; 
		}
		(*new_w)[bias_offset] = (*w)[train_bias_offset]; 
	}

	delete w; 

	return new_w; 
}


double answer_x_query(const vector<double> &x_points, const vector<double> &y_points, double x)
{
	if ((x<0)||(x>x_points[x_points.size()-1]))
	{
		// out of range
		return(-1); // negative value as null
	} else {
		int i=0; 
		while (x_points[i]<x) i++; 
		if (i==0)
		{
			return(1.0);
		} else {
			// interpolate
			double alpha = (x-x_points[i-1])/(x_points[i]-x_points[i-1]); 
			double y = y_points[i-1] + alpha*(y_points[i]-y_points[i-1]);
			return(y); 
		}
	}
}

double answer_y_query(const vector<double> &x_points, const vector<double> &y_points, double y)
{
	if ((y>1)||(y<y_points[y_points.size()-1]))
	{
		// out of range
		return(-1); // negative value as null
	} else {
		int i=0; 
		while (y_points[i]>y) i++;
		if (i==0)
		{
			return(0);
		} else {
			// interpolate
			double alpha = (y-y_points[i-1])/(y_points[i]-y_points[i-1]); 
			double x = x_points[i-1] + alpha*(x_points[i]-x_points[i-1]);
			return(x);
		}
	}
}


int main(int argc, char* argv[])
{

	vector<Example> test_sample; 
	Sparm sparm; 

	sparm.ReadParam(argc, argv); 
	// may 9: load weights and construct q-files before reading examples
	// load weights
	DenseVector *w = LoadWeights(sparm.GetModelFile(), sparm); 

	test_sample = read_input_examples(sparm.GetInputFile(), sparm); 

	
	const int maxMonth = sparm.GetMaxMonth(); 

	// re-map weights
	if (w->dim()!=sparm.GetSizePsi())
	{
		w = remap_weights(w, sparm); 
	}

	// evaluation
	double avg_l1_loss = 0;
	double avg_l2_loss = 0;
	double avg_rae_loss = 0;
	double avg_l1_log_loss = 0;
	double avg_l2_log_loss = 0; 
	double avg_log_loss = 0; 
	// ci
	vector<vector<double> > z; 

	const vector<double>& qt_time = sparm.GetQuantTime(); 

	// pre-compute contingency table
	vector<vector<double> > contingency_table(maxMonth+1, vector<double>(maxMonth+1,0)); 
	for (int j=1; j<maxMonth+1; j++)
	{
		for (int k=1; k<maxMonth+1; k++)
		{
			if (sparm.GetLossType()=="rae")	{
				contingency_table[j][k] = rae_loss(qt_time[j-1],qt_time[k-1]); 
			} else if (sparm.GetLossType()=="l1") {
				contingency_table[j][k] = l1_loss(qt_time[j-1],qt_time[k-1]); 
			} else if (sparm.GetLossType()=="l2") {
				contingency_table[j][k] = l2_loss(qt_time[j-1],qt_time[k-1]); 
			} else if (sparm.GetLossType()=="l1_log") {
				contingency_table[j][k] = l1_log_loss(qt_time[j-1],qt_time[k-1]); 
			} else if (sparm.GetLossType()=="l2_log") {
				contingency_table[j][k] = l2_log_loss(qt_time[j-1],qt_time[k-1]); 
			} else if (sparm.GetLossType()=="class") {
				double threshold = sparm.GetThreshold(); // classification
				// don't use table
				//contingency_table[j][k] = class_loss(qt_time[j-1],qt_time[k-1],threshold); 
			} else {
				//cout << "Unknown loss type '" << sparm.GetLossType() << "'!" << endl; 
				cerr << "Unknown loss type '" << sparm.GetLossType() << "'!" << endl; 
				exit(1); 
			}
		}
	}

	// accuracy and calibration
	double acc=0; 
	double mse=0; 
	long effective_sample_size=0; 

	// now classify the instances
	// first compute the probability
	for (size_t i=0; i<test_sample.size(); i++)
	{
		// return a vector a prob
		vector<double> survival_pdf = compute_survival_pdf(*w, test_sample[i].m_x.features, sparm); 

		// then use the distribution to optimize for different loss
		// can use this to fill in a (maxMonth+1)*(maxMonth+1) table
		vector<double> expected_loss(maxMonth+1, 0); 
		for (int j=1; j<maxMonth+1; j++)
		{
			for (int k=1; k<maxMonth+1; k++)
			{
				expected_loss[j] += survival_pdf[k]*contingency_table[j][k]; 
			}
		}
		// take min as prediction
		int best_prediction = 1; 
		for (int j=2; j<maxMonth+1; j++)
		{
			if (expected_loss[j]<expected_loss[best_prediction])
			{
				best_prediction=j; 
			}
		}

		double a = qt_time[best_prediction-1]; 
		// june 18: interpolation for L1-loss
		if (sparm.GetLossType()=="l1")
		{
			double p = survival_pdf[0];
			int k=1;
			while ((k<sparm.GetMaxMonth()+1)&&(p<0.5))
			{
				p += survival_pdf[k]; 
				k++;
			}
			if (k<sparm.GetMaxMonth()+1)
			{
				// linear interpolation
				double p0 = p-survival_pdf[k-1]; 
				double alpha = (0.5-p0)/(p-p0);
				double t1 = qt_time[k-1]; 
				double t0 = (k>1) ? qt_time[k-2] : 0;
				a = t0 + alpha*(t1-t0); 
			} else {
				a = qt_time[sparm.GetMaxMonth()-1]; 
			}
		}
		// end june 18
		// july 27
		// use expected survival for l2 prediction to avoid discretization effects
		if (sparm.GetLossType()=="l2")
		{
			double expected_survival = 0; 
			for (int j=1; j<sparm.GetMaxMonth()+1; j++)
			{
				expected_survival += survival_pdf[j]*qt_time[j-1]; 
			}
			a = expected_survival; 
		}
		// end july 27

		// output prediction to std out
		if (sparm.GetLossType()!="class")
		{
			//cout << qt_time[best_prediction-1]; 
			cout << a;
		} else { // classification
			int k=1; 
			double p = survival_pdf[0]; 
			while ((k<sparm.GetMaxMonth()+1)&&(qt_time[k-1]<sparm.GetThreshold()))
			{
				p += survival_pdf[k]; 
				k++; 
			}
			if (p>0.5)
			{
				cout << 1 << ":" << p; // event (death) has occurred
			} else {
				cout << 0 << ":" << p; // event has not yet occurred
			}
			/*
			if (qt_time[best_prediction-1]<sparm.GetThreshold())
			{
				cout << 0 << ":" << 1-expected_loss[best_prediction]; // event (death) has not occurred yet
			} else {
				cout << 1 << ":" << expected_loss[best_prediction];  // event has occurred
			}
			*/
		}
		
		// print probabilities for plotting survival cdf
		if (sparm.GetPrintProb()) 
		{
			double surv_prob = 1; 
			for (int j=0; j<maxMonth+1; j++) 
			{
				surv_prob -= survival_pdf[j];
				cout << ", " << surv_prob;
			}
		}
		//cout << endl; 


		// finally, answer queries for survival time and survival prob
		vector<double> x_points; 
		vector<double> y_points;
		x_points.push_back(0);
		x_points.insert(x_points.end(), qt_time.begin(), qt_time.end()); 
		double surv_prob = 1.0;
		y_points.push_back(surv_prob);
		for (int j=0; j<maxMonth; j++)
		{
			surv_prob -= survival_pdf[j]; 
			y_points.push_back(surv_prob); 
		}
		vector<double> x_queries = sparm.GetXQueries(); 
		vector<double> y_queries = sparm.GetYQueries(); 
		for (size_t j=0; j<x_queries.size(); j++)
		{
			double y = answer_x_query(x_points, y_points, x_queries[j]); 
			if (y>=0)
			{
				cout << ", x:" << x_queries[j] << ":" << y; 
			} else {
				cout << ", x:" << x_queries[j] << ":OUT_OF_RANGE"; 
			}
		}
		for (size_t j=0; j<y_queries.size(); j++)
		{
			double x = answer_y_query(x_points, y_points, y_queries[j]); 
			if (x>=0)
			{
				cout << ", y:" << y_queries[j] << ":" << x; 
			} else {
				cout << ", y:" << y_queries[j] << ":OUT_OF_RANGE";
			}
		}
		// print probability at event time, if '-p' option is specified
		if (sparm.GetPrintProb()) 
		{
			double t = test_sample[i].m_y.original_survival_time; 
			double y = answer_x_query(x_points, y_points, t);
			if (y>=0)
			{
				cout << ", t:" << t << ":" << y; 
			} else {
				cout << ", t:" << t << ":OUT_OF_RANGE"; 
			}

		}
		cout << endl; 


		// accuracy and calibration
		if (sparm.GetLossType()=="class")
		{
			int k=1; 
			double p = survival_pdf[0]; 
			while ((k<sparm.GetMaxMonth()+1)&&(qt_time[k-1]<sparm.GetThreshold()))
			{
				p += survival_pdf[k]; 
				k++; 
			}
			// survival prob
			p = 1-p; 
			if (!(((test_sample[i].m_y.censored)&&(sparm.GetTrainUncensored()==0))&&(test_sample[i].m_y.original_survival_time<sparm.GetThreshold())))
			{
				effective_sample_size++; 
				if ((test_sample[i].m_y.original_survival_time-sparm.GetThreshold())*(p-0.5)>0)
				{
					acc+=1.0; 
				}
				int label = (test_sample[i].m_y.original_survival_time>sparm.GetThreshold()) ? 1 : 0;
				mse += (label-p)*(label-p); 
			} 
		}

		double b = test_sample[i].m_y.original_survival_time;
		// ci
		vector<double> temp;
		temp.push_back(b); 
		temp.push_back(a); 
		temp.push_back(test_sample[i].m_y.censored ? 1.0 : 0.0); 
		z.push_back(temp); 

		if ((a<=b)||(!((test_sample[i].m_y.censored)&&(sparm.GetTrainUncensored()==0)))) {
			avg_l1_loss += l1_loss(a,b); 
			avg_l2_loss += l2_loss(a,b); 
			avg_rae_loss += rae_loss(a,b); 
			avg_l1_log_loss += l1_log_loss(a,b); 
			avg_l2_log_loss += l2_log_loss(a,b); 
		} // else the loss is 0; do nothing

		if (test_sample[i].m_y.censored)
		{
			double sum=0; 
			for (int j=test_sample[i].m_y.event_time; j<maxMonth+1; j++)
			{
				sum += survival_pdf[j];
			}
			avg_log_loss -= log(sum); 
		} else {
			avg_log_loss -= log(survival_pdf[test_sample[i].m_y.event_time]); 
		}
	}

	if (sparm.GetLossType()=="class")
	{
		cout << "#avg acc at threshold " << sparm.GetThreshold() << ": " << acc/effective_sample_size << endl; 
		cout << "#avg mse at threshold " << sparm.GetThreshold() << ": " << mse/effective_sample_size << endl; 
	}
	// ci
	double ci = concordanceIndex(z); 

	// output summary 
	cout << "#concordance index: " << ci << endl; 
	cout << "#avg l1-loss: " << avg_l1_loss/test_sample.size() << endl; 
	//cout << "#avg l2-loss: " << avg_l2_loss/test_sample.size() << endl; 
	cout << "#avg l2-loss: " << sqrt(avg_l2_loss/test_sample.size()) << endl; 
	cout << "#avg rae-loss: " << avg_rae_loss/test_sample.size() << endl; 
	cout << "#avg l1-log-loss: " << avg_l1_log_loss/test_sample.size() << endl; 
	cout << "#avg l2-log-loss: " << avg_l2_log_loss/test_sample.size() << endl; 
	cout << "#avg log-likelihood loss: " << avg_log_loss/test_sample.size() << endl; 

	delete w;

	return(0);
}
