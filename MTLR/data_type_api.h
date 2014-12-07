
#pragma once

#include "SparseVector.h"
#include "Sparm.h"


struct Pattern
{
	SparseVector features;
	Pattern(const SparseVector& input_features): features(input_features) {}; 
};

struct Label
{
  int event_time; // from 0 up to MAX_MONTH; MAX_MONTH could be length of study, or some upper bound on the survival period that we're interested in
  bool censored; 
  double original_survival_time; // original survival time without quantization; only needed for testing
  Label(int n, bool c, double s): event_time(n), censored(c), original_survival_time(s) {}; 
};

struct Example
{
	Pattern m_x;
	Label m_y; 
	Example(const SparseVector &input_features, int n, bool c, double s): m_x(input_features), m_y(n,c,s) {}
}; 

