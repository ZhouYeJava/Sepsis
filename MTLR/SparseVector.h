#pragma once

#include <vector>
#include <iterator>
#include <iostream>

class DenseVector;

using namespace std;

class SparseVector
{
	friend std::ostream& operator<<(std::ostream&, const SparseVector&); 
private:
	vector<pair<size_t, double> > m_svector; 
public:
	typedef vector<pair<size_t, double> >::const_iterator const_iterator; 
	const_iterator begin() const {return m_svector.begin(); }
	const_iterator end() const {return m_svector.end(); }

	void push_back(size_t key, double value)
	{
		m_svector.push_back(make_pair(key, value)); 
	}

	SparseVector(void); 
	SparseVector(vector<pair<size_t, double> > words);

	SparseVector(const SparseVector &copy); 

	~SparseVector(void);

	void append(const SparseVector &to_add);

};


SparseVector multadd_ss(const SparseVector &a, const SparseVector &b, double factor); 

double sprod_ns(const DenseVector &w, const SparseVector &b, size_t offset=0); 

void multadd_ns(DenseVector &w, const SparseVector &b, double factor, size_t offset=0); 
