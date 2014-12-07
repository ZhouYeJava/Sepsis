#include "DenseVector.h"

#include <cassert>
#include <iostream>
#include <stdio.h>

DenseVector::DenseVector(size_t n): m_dim(n)
{
  m_dvector.reserve(n);
  for (size_t i=0; i<n; i++)
  {
	  m_dvector.push_back(0);
  }
}

DenseVector::DenseVector(void)
:m_dim(0)
{
}

DenseVector::~DenseVector(void)
{
}


void DenseVector::clear()
{
  for (size_t i=0; i<m_dim; i++)
  {
	  m_dvector[i] = 0;
  }
}


double sprod_nn(const DenseVector &a, const DenseVector &b) 
{
	double ans = 0; 
	assert(a.dim() == b.dim());
	size_t n = a.dim(); 

	for (size_t i=1; i<n; i++)
	{
		ans += a[i]*b[i]; 
	}
	return ans; 
}

void multadd_nn(DenseVector &w, const DenseVector &a, double factor)
{
	assert(w.dim()==a.dim());
	size_t n = w.dim();

	for (size_t i=1; i<n; i++)
	{
		w[i] += factor*a[i];
	}

}


void smult_n(DenseVector &w, double factor)
{
	size_t n = w.dim(); 
	for (size_t i=1; i<n; i++)
	{
		w[i] *= factor;
	}
}


