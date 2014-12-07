
#include <vector>
#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <stdio.h>

#include "common.h"

#define LINE_SEARCH_C1 1E-4
#define LINE_SEARCH_C2 0.9
#define ALPHA_MAX 1024
#define MAX_ITER 1000
#define L1_MAX_ITER 10000
#define ETA 1.25



using namespace std;

// l2
double compute_gradient_obj_l2(DenseVector &w, const vector<Example> &sample, DenseVector &g, const Sparm &sparm); 
double wolfe_line_search(const DenseVector &w, double v, const DenseVector &g, const DenseVector &r, const vector<Example> &sample, const Sparm &sparm);
void train_l2(const DenseVector& w_initial, const vector<Example> &sample, const Sparm &sparm);

// l1
double compute_gradient_obj_l1(DenseVector &w, const vector<Example> &sample, DenseVector &g, const vector<vector<double > > &h_marginals, const Sparm &sparm); 
double backtrack_ista(DenseVector&w, const DenseVector &r, const vector<Example> &sample, double v, double &L, const vector<vector<double > > &h_marginals, const Sparm &sparm);
vector<vector<double> > compute_hidden_marginals(const DenseVector &w, const vector<Example> &sample, const Sparm &sparm);
void train_l1(const DenseVector& w_initial, const vector<Example> &sample, const Sparm &sparm);



// compute gradient for the smoothing term
void smoothing_grad(DenseVector &w, DenseVector &g, const Sparm &sparm)
{
	size_t oriSizePsi = sparm.GetOriginalSizePsi();
	double C2 = sparm.GetC2(); 

	for (int i=1; i<sparm.GetMaxMonth()-1; i++)
	{
		size_t offset1 = i*oriSizePsi; 
		size_t offset2 = (i+1)*oriSizePsi; 
		size_t offset3 = (i-1)*oriSizePsi; 

		for (size_t j=0; j<oriSizePsi; j++)
		{
			g[offset1+j] += C2*(w[offset1+j] - w[offset2+j]);
			g[offset1+j] += C2*(w[offset1+j] - w[offset3+j]);
		}
	}
	// boundary cases
	size_t offset1 = 0;
	size_t offset2 = oriSizePsi;
	for (size_t j=0; j<oriSizePsi; j++)
	{
		g[offset1+j] += C2*(w[offset1+j] - w[offset2+j]);
	}

	offset1 = (sparm.GetMaxMonth()-1)*oriSizePsi;
	offset2 = (sparm.GetMaxMonth()-2)*oriSizePsi;
	for (size_t j=0; j<oriSizePsi; j++)
	{
		g[offset1+j] += C2*(w[offset1+j] - w[offset2+j]);
	}

}


// multi-task smoothing term 
double smoothing_obj(DenseVector &w, const Sparm &sparm)
{
	double ans=0; 
	size_t oriSizePsi = sparm.GetOriginalSizePsi(); 

	for (int i=0; i<sparm.GetMaxMonth()-1; i++)
	{
		size_t offset1 = i*oriSizePsi; 
		size_t offset2 = (i+1)*oriSizePsi; 
		for (size_t j=0; j<oriSizePsi; j++)
		{
			ans += (w[offset1+j] - w[offset2+j])*(w[offset1+j] - w[offset2+j]);
		}
	}
	return(ans); 
}



int WriteWeights(const DenseVector &w, string filename, const Sparm &sparm)
{
	ofstream outfile(filename.c_str());
	if (outfile) {
		// may 9: new default for q-file
		outfile << "m:" << sparm.GetMaxMonth() << endl; 
		const vector<double> quantTime = sparm.GetQuantTime(); 
		outfile << quantTime[0]; 
		for (size_t i=1; i<quantTime.size(); i++)
		{
			outfile << "," << quantTime[i]; 
		}
		outfile << endl; 
		// end may 9
		// aug 2: regularization parameter, L1/L2
		outfile << "r:" << sparm.GetRegularizer() << endl; 
		// end aug 2
		outfile << "DIM:" << w.dim() << std::endl; 
		for (size_t i=1; i<w.dim(); i++)
		{
			outfile << i << ":" << setprecision(8) << w[i] << endl; 
		}
		outfile.close(); 
		return 0;
	} else {
		cerr << "Unable to open model file " << filename << " for output!" << endl;  
		return -1;
	}
}





int main(int argc, char* argv[])
{
	vector<Example> sample; 
	Sparm sparm; 

	sparm.ReadParam(argc, argv); 
	
	DenseVector *w=NULL;
	size_t sizePsi = 0; 
	// july 20: allow uncensored training to use a weight file for initialization, to speed up cross validation
	//if ((sparm.GetTrainUncensored()==0)&&(sparm.GetInitialWeightFile()!=""))
	if (sparm.GetInitialWeightFile()!="")
	{
		w = LoadWeights(sparm.GetInitialWeightFile(), sparm); // q intervals included in weight file
		sample = read_input_examples(sparm.GetInputFile(), sparm); 
		sizePsi = sparm.GetSizePsi(); 
		assert(w->dim()==sizePsi); 
	} else {
		// may 9: new default for q file
		if (sparm.GetQuantTime().size()==0)
		{
			construct_intervals(sparm.GetInputFile(), sparm); 
		}
		// end may 9
		sample = read_input_examples(sparm.GetInputFile(), sparm); 
		sizePsi = sparm.GetSizePsi(); 
		w = new DenseVector(sizePsi); 
	}

	if (sparm.GetRegularizer()==2)
	{
		train_l2(*w, sample, sparm); 
	} else {
		train_l1(*w, sample, sparm); 
	}

	delete w; 
	return(0); 
}


void train_l2(const DenseVector& w_initial, const vector<Example> &sample, const Sparm &sparm)
{
	// variables for L-BFGS
	const size_t bundle_size = sparm.GetBundleSize(); 
	list<DenseVector*> s_list;
	list<DenseVector*> y_list;
	list<double> rho_list; 
	DenseVector *last_g=NULL, *last_w=NULL; 
	DenseVector *g=NULL;
	DenseVector *w = new DenseVector(w_initial); 

	bool use_lbfgs = true; 
	size_t sizePsi = sparm.GetSizePsi(); 


	int outer_iter=0; 
	double v=1;
	double last_v=1;
	double relative_decrease=1;
	while ((outer_iter<MAX_ITER)&&((outer_iter<1)||(relative_decrease>1E-6)))
	{
		cout << "Iteration: " << outer_iter << endl; 

		delete last_g;
		last_g = g;
		last_v = v; 

		g = new DenseVector(sizePsi); 

		v = compute_gradient_obj_l2(*w, sample, *g, sparm);
		cout << "objective (with regularizer): " << v << endl; 

		if (outer_iter>0) {
			relative_decrease = fabs((v-last_v)/last_v);
		} else {
			relative_decrease = 1; 
		}


		DenseVector *r = new DenseVector(sizePsi);

		if ((!use_lbfgs)||(outer_iter==0))
		{
			// search direction: r = -g, steepest descent
			multadd_nn(*r, *g, -1); 
		} else {
			// update s, y, rho
			s_list.push_back(new DenseVector(*w));
			y_list.push_back(new DenseVector(*g));
			multadd_nn(*(s_list.back()), *last_w, -1.0f); 
			multadd_nn(*(y_list.back()), *last_g, -1.0f);
			double sTy = sprod_nn(*(s_list.back()), *(y_list.back())); 
			double yTy = sprod_nn(*(y_list.back()), *(y_list.back())); 
			rho_list.push_back(1.0f/sTy); 
			double gamma_k = sTy/yTy;  
			if (rho_list.size()>bundle_size)
			{
				rho_list.erase(rho_list.begin());
				delete *(s_list.begin()); 
				s_list.erase(s_list.begin());
				delete *(y_list.begin()); 
				y_list.erase(y_list.begin()); 
			}

			list<double> a; 
			list<double>::const_iterator iter_rho = rho_list.end(); 
			list<DenseVector*>::const_iterator iter_s = s_list.end(); 
			list<DenseVector*>::const_iterator iter_y = y_list.end(); 

			multadd_nn(*r, *g, 1); 
			while (iter_rho!=rho_list.begin())
			{
				iter_rho--;
				iter_s--;
				iter_y--;
  
				double a_i = (*iter_rho)*sprod_nn(*(*iter_s), *r); 
				a.push_front(a_i); 
				multadd_nn(*r, *(*iter_y), -a_i);
			}

			smult_n(*r, 1.0f/gamma_k); 

			iter_rho = rho_list.begin(); 
			iter_s = s_list.begin(); 
			iter_y = y_list.begin(); 
			list<double>::const_iterator iter_a = a.begin(); 

			while (iter_rho!=rho_list.end())
			{
				double b = (*iter_rho)*sprod_nn(*r, *(*iter_y)); 
				multadd_nn(*r, *(*iter_s), (*iter_a)-b); 

				iter_rho++;
				iter_s++;
				iter_y++;
				iter_a++;  
			}

			smult_n(*r, -1.0f); 
		}

		double alpha = wolfe_line_search((const DenseVector &)*w, v, *g, *r, sample, sparm); 

		delete last_w;
		last_w = new DenseVector(*w); // copy before increment
		multadd_nn(*w, *r, alpha); 

		delete r; 

		outer_iter++; 
	}

	// write weights to file
	WriteWeights(*w, sparm.GetModelFile(),sparm); 

	// CLEANUP
	for (list<DenseVector*>::iterator iter_s = s_list.begin(); iter_s!=s_list.end(); ++iter_s)
	{
		delete *iter_s;
	}
	for (list<DenseVector*>::iterator iter_y = y_list.begin(); iter_y!=y_list.end(); ++iter_y)
	{
		delete *iter_y;
	}
	delete last_w; 
	delete last_g;
	delete g; 
	delete w; 

}



void train_l1(const DenseVector& w_initial, const vector<Example> &sample, const Sparm &sparm)
{
	// ISTA
	// deal with the no censoring case first
	size_t sizePsi = sparm.GetSizePsi(); 
	DenseVector *w = new DenseVector(w_initial); 
	DenseVector *last_w = new DenseVector(sizePsi); 
	DenseVector *g = new DenseVector(sizePsi); 

	// L0; or use an estimate based on norm of x_i?
	// figuring out a good Lipschitz constant could be the key here
	//double L = 1.0;

	// mar 26: EM-wrapper
	int outer_iter = 0; 
	double outer_v = 0;
	double outer_last_v = 0;
	while ((outer_iter<MAX_ITER)&&((outer_iter<2)||(abs(outer_last_v/outer_v-1)>1E-5)))
	{
		outer_iter++;

		// E-step
		vector<vector<double> > h_marginals;
		if (sparm.GetTrainUncensored()==0)
		{
			h_marginals = compute_hidden_marginals(*w, sample, sparm); 
		}


		double v = compute_gradient_obj_l1(*w, sample, *g, h_marginals, sparm);
		double last_v = v;

		// initialize L with norm of initial gradient
		double L = sqrt(sprod_nn(*g,*g)); 
		int iter=0;

		// try an initial stopping criterion based on decrease in objective
		bool use_fista=true; 
		if (!use_fista)
		{
			//while ((iter<L1_MAX_ITER)&&((abs((last_v/v-1)>1E-6)||(iter==0)) ))
			while ((iter<L1_MAX_ITER)&&((abs((last_v/v-1)>1E-5)||(iter==0)) ))
			{
				cout << "iter =  " << iter << ", L = " << L << endl;
				cout << "v = " << v << endl; 
				last_v = v;
				g->clear(); 
				v = backtrack_ista(*w, *g, sample, v, L, h_marginals, sparm);
				double temp_v = compute_gradient_obj_l1(*w, sample, *g, h_marginals, sparm);
				iter++;
			}
			// END ISTA
		} else {
			// NOW TRY FISTA
			double t=1.0;

			//DenseVector *y = new DenseVector(w->dim());
			DenseVector *y = new DenseVector(*w);
			//while ((iter<L1_MAX_ITER)&&((abs((last_v/v-1)>1E-6)||(iter==0)) ))
			while ((iter<L1_MAX_ITER)&&((abs((last_v/v-1)>1E-5)||(iter==0)) ))
			{
				cout << "iter: " << iter << ", L: " << L << ", v: " << v << endl; 
				g->clear(); 
				double temp_v = compute_gradient_obj_l1(*y, sample, *g, h_marginals, sparm);
				last_v = v;
				v = backtrack_ista(*y, *g, sample, v, L, h_marginals, sparm); 
				delete last_w;
				last_w = w; 
				w = y;

				double last_t = t;
				t = (1.0+sqrt(1 + 4*t*t))/2;

				y = new DenseVector(w->dim());
				multadd_nn(*y, *w, 1.0+(last_t-1)/t);
				multadd_nn(*y, *last_w, -(last_t-1)/t); 

				iter++;
			}
			delete y;
			//delete last_w;
		}

		outer_last_v = outer_v;
		outer_v = v; 

		// break if not training with censored values
		if (sparm.GetTrainUncensored()!=0) break;

		// END FISTA
	} // end outer loop
	delete last_w; 

	// write weights to file
	WriteWeights(*w, sparm.GetModelFile(), sparm); 

	delete g;
	delete w; 

}




double compute_gradient_obj_l2(DenseVector &w, const vector<Example> &sample, DenseVector &g, const Sparm &sparm)
{
	double v=0;
	int maxMonth = sparm.GetMaxMonth(); 
	double b_constant=100;
	if (sparm.GetRegularizer()==1) b_constant=1;

	for (size_t example_id=0; example_id<sample.size(); ++example_id)
	{
		int event_time = sample[example_id].m_y.event_time;  // 0 <= t <= maxMonth; t==maxMonth means no event during the observed period

		// there are maxMonth+1 possible sequences
		vector<double> scores(maxMonth+1, 0);
		double sum=0; 
		for (int i=1; i<maxMonth+1; i++)
		{
			size_t bias_offset = maxMonth*sparm.GetOriginalSizePsi() + (i-1); 
			size_t offset = (i-1)*sparm.GetOriginalSizePsi(); 
			scores[i] = sprod_ns(w,sample[example_id].m_x.features, offset) + w[bias_offset]*b_constant; // need offset here
			// sum stores the scores of the sequence of all 1's
			sum += scores[i]; 
		}

		vector<double> running_scores(maxMonth+1,0); 
		running_scores[0] = sum; 

		double max_score = running_scores[0]; 
		for (int i=1; i<maxMonth+1; i++)
		{
			running_scores[i] = running_scores[i-1] - scores[i]; // flipping signs
			if (running_scores[i]>max_score) max_score = running_scores[i]; 
		}

		double z=0; 
		vector<double> marginals(maxMonth+1, 0);
		for (int i=0; i<maxMonth+1; i++)
		{
			z += exp(running_scores[i]-max_score); 
			marginals[i] = z; // or cdf; marginals[i] stores p(t<i)
		}
		double log_z = log(z) + max_score; 
		v += log_z; 

		// for censored targets
		double z_h=0; 
		double log_z_h=0; 
		vector<double> marginals_h(maxMonth+1, 0);
		double max_score_h=0; 
		if ((sample[example_id].m_y.censored)&&(sparm.GetTrainUncensored()==0))
		{
			max_score_h = running_scores[event_time]; 
			for (int i=event_time+1; i<maxMonth+1; i++)
			{
				if (running_scores[i]>max_score_h) max_score_h = running_scores[i]; 
			}
			for (int i=event_time; i<maxMonth+1; i++)
			{
				z_h += exp(running_scores[i]-max_score_h); 
				marginals_h[i] = z_h; // or cdf; marginals[i] stores p(t<i)
			}
			log_z_h = log(z_h) + max_score_h; 
			v -= log_z_h; 
		} else {
			v -= running_scores[event_time]; 
		}

		for (int i=1; i<maxMonth+1; i++) 
		{
			size_t offset = (i-1)*sparm.GetOriginalSizePsi(); 
			double p = exp(log(marginals[i-1])+max_score-log_z);
			double factor = p; 
			if ((sample[example_id].m_y.censored)&&(sparm.GetTrainUncensored()==0))
			{
				if (event_time<i) {
					double p_h = exp(log(marginals_h[i-1])+max_score_h-log_z_h);
					factor -= p_h; 
				}
			} else {
				if (event_time<i) factor -= 1; 
			}
			multadd_ns(g, sample[example_id].m_x.features, factor, offset); // need offset here
			size_t bias_offset = maxMonth*sparm.GetOriginalSizePsi() + (i-1); 
			g[bias_offset] += factor*b_constant; 
		}
	}
	v = v*sparm.GetC1()/sample.size(); 
	smult_n(g, sparm.GetC1()/sample.size()); 

	v += 0.5*sprod_nn(w,w); 
	multadd_nn(g, w, 1.0); 

	v += 0.5*sparm.GetC2()*smoothing_obj(w, sparm);
	smoothing_grad(w, g, sparm); 

	return(v); 

}


double cubic_interpolation(double a1, double a2, double f1, double f2, double g1, double g2)
{
	double d1 = g1 + g2 - 3*(f1-f2)/(a1-a2); 
	double d2 = sqrt(d1*d1-g1*g2); 
	if (a1>a2) {
		d2*=-1;
	}
	return (a2 - (a2-a1)*(g2+d2-d1)/(g2-g1+2*d2));
}


double zoom(double a_lo, double a_hi, double f_lo, double f_hi, double g_lo, double g_hi, const DenseVector &w, const DenseVector &r,
		   const vector<Example> &sample, double v, double suff_decrease_value, const Sparm &sparm)
{
	double alpha_lo = a_lo;
	double alpha_hi = a_hi; 

	int iter=0;
	double trial_obj=0;
	double alpha_j=0;
	while (iter<100)
	{
		iter++;

		alpha_j = cubic_interpolation(alpha_lo, alpha_hi, f_lo, f_hi, g_lo, g_hi); 

		DenseVector trial_pt(w); 
		multadd_nn(trial_pt, r, alpha_j);

		DenseVector trial_gradient(w.dim());
		trial_obj = compute_gradient_obj_l2(trial_pt, sample, trial_gradient, sparm); 
		double trial_decrease_value = sprod_nn(r, trial_gradient); 

		cout << "alpha_lo: " << alpha_lo << ", alpha_hi: " << alpha_hi << ", alpha_j: " << alpha_j << ", trial obj: " << trial_obj << endl; 

		if ((trial_obj>v+LINE_SEARCH_C1*alpha_j*suff_decrease_value)||(trial_obj>=f_lo))
		{
			alpha_hi = alpha_j; 
			f_hi = trial_obj; 
			g_hi = trial_decrease_value; 
		} else {
			if (fabs(trial_decrease_value)<=-LINE_SEARCH_C2*suff_decrease_value)
			{
				return alpha_j;
			}
			if (trial_decrease_value*(alpha_hi-alpha_lo)>=0) 
			{
				alpha_hi = alpha_lo;
				f_hi = f_lo;
				g_hi = g_lo;
			}
			alpha_lo = alpha_j;
			f_lo = trial_obj;
			g_lo = trial_decrease_value;

		} // end else

	} // end while 

	// mar 12: allow condition 2 of Wolfe to fail
	if (trial_obj<v+LINE_SEARCH_C1*alpha_j*suff_decrease_value)
	{
		return alpha_j; 
	}

	//cout << "zoom() failed!" << endl; 
	cerr << "zoom() failed!" << endl;
	exit(1); 
}

double wolfe_line_search(const DenseVector &w, double v, const DenseVector &g, const DenseVector &r, const vector<Example> &sample, const Sparm &sparm)
{
	double suff_decrease_value;
	double a_l, f_l, g_l, a_c, f_c, g_c; 

	suff_decrease_value = sprod_nn(r,g); 
	size_t n = w.dim(); 

	int iter = 0; 
	a_l = 0; 
	f_l = v;
	g_l = suff_decrease_value; 

	a_c = 1;
	while (iter<100)
	{
		iter++; 
		DenseVector trial_pt(w); // copy from w
		multadd_nn(trial_pt, r, a_c); 

		DenseVector trial_gradient(n); 

		f_c = compute_gradient_obj_l2(trial_pt, sample, trial_gradient, sparm); 
		g_c = sprod_nn(r, trial_gradient); 

		cout << "in wolfe, v: " << v << ", alpha: " << a_c << ", trial_obj: " << f_c << endl; 

		if ((f_c>v+LINE_SEARCH_C1*a_c*suff_decrease_value)||((iter>1)&&(f_c>f_l)))
		{
			return zoom(a_l,a_c,f_l,f_c,g_l,g_c,w,r,sample,v,suff_decrease_value,sparm); 
		}

		if (fabs(g_c)<=-LINE_SEARCH_C2*suff_decrease_value)
		{
			return a_c; 
		}

		if (g_c>=0)
		{
			return zoom(a_c,a_l,f_c,f_l,g_c,g_l,w,r,sample,v,suff_decrease_value,sparm);
		}

		a_l = a_c;
		f_l = f_c; 
		g_l = g_c;

		a_c = min(2*a_c, (double) ALPHA_MAX);
	}
	//cout << "wolfe_line_search() failed! " << endl; 
	cerr << "wolfe_line_search() failed!" << endl; 
	exit(1); 

}






vector<vector<double> > compute_hidden_marginals(const DenseVector &w, const vector<Example> &sample, const Sparm &sparm)
{
	vector<vector<double> > p_h; 
	double b_constant = 100;
	if (sparm.GetRegularizer()==1) b_constant = 1; 

	int maxMonth = sparm.GetMaxMonth(); 
	for (size_t example_id=0; example_id<sample.size(); ++example_id)
	{
		int event_time = sample[example_id].m_y.event_time;  // 0 <= t <= maxMonth; t==maxMonth means no event during the observed period
		if (sample[example_id].m_y.censored)
		{
			vector<double> scores(maxMonth+1, 0);
			double sum=0; 
			for (int i=1; i<maxMonth+1; i++)
			{
				size_t bias_offset = maxMonth*sparm.GetOriginalSizePsi() + (i-1); 
				size_t offset = (i-1)*sparm.GetOriginalSizePsi(); 
				scores[i] = sprod_ns(w,sample[example_id].m_x.features, offset) + w[bias_offset]*b_constant; // need offset here
				// sum stores the scores of the sequence of all 1's
				sum += scores[i]; 
			}

			vector<double> running_scores(maxMonth+1,0); 
			running_scores[0] = sum; 
			double max_score = running_scores[0]; 
			for (int i=1; i<maxMonth+1; i++)
			{
				running_scores[i] = running_scores[i-1] - scores[i]; // flipping signs
				if (running_scores[i]>max_score) max_score = running_scores[i]; 
			}

			double z_h = 0;
			double log_z_h = 0;
			double max_score_h = running_scores[event_time]; 
			for (int i=event_time+1; i<maxMonth+1; i++)
			{
				if (running_scores[i]>max_score_h) max_score_h = running_scores[i]; 
			}
			for (int i=event_time; i<maxMonth+1; i++)
			{
				z_h += exp(running_scores[i]-max_score_h); 
			}
			log_z_h = log(z_h) + max_score_h; 

			vector<double> marginals_h;
			for (int i=event_time; i<maxMonth+1; i++)
			{
				marginals_h.push_back(exp(running_scores[i]-log_z_h)); 
			}
			p_h.push_back(marginals_h); 
		} else {
			p_h.push_back(vector<double>()); // empty
		}
	}
	return(p_h);
}





double compute_gradient_obj_l1(DenseVector &w, const vector<Example> &sample, DenseVector &g, const vector<vector<double > > &h_marginals, const Sparm &sparm)
{
	double v=0;
	int maxMonth = sparm.GetMaxMonth(); 
	double b_constant = 100;
	if (sparm.GetRegularizer()==1) b_constant=1; 

	for (size_t example_id=0; example_id<sample.size(); ++example_id)
	{
		int event_time = sample[example_id].m_y.event_time;  // 0 <= t <= maxMonth; t==maxMonth means no event during the observed period

		vector<double> scores(maxMonth+1, 0);
		double sum=0; 
		for (int i=1; i<maxMonth+1; i++)
		{
			size_t bias_offset = maxMonth*sparm.GetOriginalSizePsi() + (i-1); 
			size_t offset = (i-1)*sparm.GetOriginalSizePsi(); 
			scores[i] = sprod_ns(w,sample[example_id].m_x.features, offset) + w[bias_offset]*b_constant; // need offset here
			// sum stores the scores of the sequence of all 1's
			sum += scores[i]; 
		}

		vector<double> running_scores(maxMonth+1,0); 
		running_scores[0] = sum; 
		double max_score = running_scores[0]; 
		for (int i=1; i<maxMonth+1; i++)
		{
			running_scores[i] = running_scores[i-1] - scores[i]; // flipping signs
			if (running_scores[i]>max_score) max_score = running_scores[i]; 
		}

		double z=0; 
		vector<double> marginals(maxMonth+1, 0);
		for (int i=0; i<maxMonth+1; i++)
		{
			z += exp(running_scores[i]-max_score); 
			marginals[i] = z; // or cdf; marginals[i] stores p(t<i)
		}
		double log_z = log(z) + max_score; 
		v += log_z; 

		// gradient
		// for censored targets
		if ((sample[example_id].m_y.censored)&&(sparm.GetTrainUncensored()==0))
		{
			for (int i=event_time; i<maxMonth+1; i++)
			{
				v -= h_marginals[example_id][i-event_time]*running_scores[i]; 
			}

		} else {
			v -= running_scores[event_time]; 
		}

		double sum_h = 0;
		if ((sample[example_id].m_y.censored)&&(sparm.GetTrainUncensored()==0))
		{
			sum_h = h_marginals[example_id][0]; 
		}

		for (int i=1; i<maxMonth+1; i++) 
		{
			size_t offset = (i-1)*sparm.GetOriginalSizePsi(); 
			
			double p = exp(log(marginals[i-1])+max_score-log_z);
			double factor = p; 
			if ((sample[example_id].m_y.censored)&&(sparm.GetTrainUncensored()==0))
			{
				if (event_time<i) {
					sum_h += h_marginals[example_id][i-event_time];
					factor -= sum_h; 
				}
			} else {
				if (event_time<i) factor -= 1; 
			}

			multadd_ns(g, sample[example_id].m_x.features, factor, offset); // need offset here
			size_t bias_offset = maxMonth*sparm.GetOriginalSizePsi() + (i-1); 
			g[bias_offset] += factor*b_constant; 
		}
	}
	v = v*sparm.GetC1()/sample.size(); 
	smult_n(g, sparm.GetC1()/sample.size()); 

	v += 0.5*sparm.GetC2()*smoothing_obj(w, sparm);
	smoothing_grad(w, g, sparm); 

	return(v); 

}


double take_threshold(DenseVector &w, double t, const Sparm &sparm)
{
	assert(t>0);
	double ans = 0;
	//for (size_t i=0; i<w.dim(); i++)
	for (size_t i=0; i<sparm.GetMaxMonth()*sparm.GetOriginalSizePsi(); i++)
	{
		if (w[i]>t)
		{
			w[i] = w[i] - t;
			ans += w[i];
		} else {
			if (w[i]<-t) 
			{
				w[i] = w[i] + t;
				ans -= w[i];
			} else {
				w[i] = 0;
			}
		}
	}
	return(ans);
}


//double compute_l1_norm(const DenseVector &w)
double compute_l1_norm(const DenseVector &w, const Sparm &sparm)
{
	double ans=0;
	for (size_t i=0; i<sparm.GetMaxMonth()*sparm.GetOriginalSizePsi(); i++)
	{
		ans += fabs(w[i]);
	}
	return(ans);
}


double backtrack_ista(DenseVector&w, const DenseVector &r, const vector<Example> &sample, double v, double &L, const vector<vector<double > > &h_marginals, const Sparm &sparm)
{
// v: current obj
// L: Lipschitz constant (starting)
// during backtrack search you only need the objective, not the gradient; can change to a cheaper version later
	double Lk = L;
	double temp_v = 0;
	DenseVector trial_pt(w.dim()); 
	int iter=0; 

	DenseVector temp_g(w.dim());
	double fw = compute_gradient_obj_l1(w, sample, temp_g, h_marginals, sparm); 

	while (iter<200) 
	{
		trial_pt = w;
		DenseVector trial_g(w.dim());

		multadd_nn(trial_pt, r, -1.0/Lk);
		// threshold and return l1 norm
		double l1_norm = take_threshold(trial_pt, 1.0/Lk, sparm); 
		double trial_obj = compute_gradient_obj_l1(trial_pt, sample, trial_g, h_marginals, sparm);

		DenseVector diff(trial_pt);
		multadd_nn(diff, w, -1.0); 

		// original form
		//if (l1_norm + trial_obj < v) 
		//if (trial_obj < fw - 0.5/Lk*g_norm) 
		if (trial_obj < fw + sprod_nn(diff,r) + 0.5*Lk*sprod_nn(diff,diff) + 1E-10)
		{
			temp_v = l1_norm + trial_obj; 
			break;
		}
		Lk = Lk*ETA;
		iter++;
	}
	if (iter==200)
	{
		cerr << "Backtrack search failed in ISTA!" << endl; 
		exit(1);
	}

	w = trial_pt;
    L = Lk;
	return(temp_v);
}






