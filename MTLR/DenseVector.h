#pragma once

#include <vector>

class DenseVector
{
protected:
	std::vector<double> m_dvector;
	size_t m_dim; 
public:
	typedef std::vector<double>::const_iterator const_iterator; 
	DenseVector(void);
	~DenseVector(void);

	const_iterator begin() const {return m_dvector.begin(); }
	const_iterator end() const {return m_dvector.end(); }

	DenseVector(size_t n); 

	void push_back(double v)
	{
		m_dvector.push_back(v);
		m_dim++; 
	}

	double const& operator[](const size_t i) const
	{
		return(m_dvector[i]); 
	}

	double& operator[](const size_t i)
	{
		return(m_dvector[i]); 
	}

	void clear(); 

	size_t dim() const
	{
		return(m_dim);
	}

};

double sprod_nn(const DenseVector &a, const DenseVector &b);

void multadd_nn(DenseVector &w, const DenseVector &a, double factor); 

void smult_n(DenseVector &w, double factor); 
