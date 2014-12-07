#include "SparseVector.h"
#include "DenseVector.h"
#include <stdio.h>

SparseVector::SparseVector(vector<pair<size_t, double> > words)
:m_svector(words)
{
}

SparseVector::SparseVector(void)
{
}

SparseVector::SparseVector(const SparseVector &copy)
:m_svector(copy.m_svector)
{
}


SparseVector::~SparseVector(void)
{
}

void SparseVector::append(const SparseVector &to_add)
{
	m_svector.insert(m_svector.end(), to_add.m_svector.begin(), to_add.m_svector.end()); 
}	


SparseVector multadd_ss(const SparseVector &a, const SparseVector &b, double factor)
{
	vector<pair<size_t,double> > words; 

	SparseVector::const_iterator iter_a = a.begin();
	SparseVector::const_iterator iter_b = b.begin(); 

	while (iter_a!=a.end() && iter_b!=b.end()) 
	{
		if (iter_a->first > iter_b->first)
		{
			words.push_back(make_pair(iter_b->first, factor*iter_b->second));
			iter_b++;
		} else {
			if (iter_a->first < iter_b->first)
			{
				words.push_back(*iter_a);
				iter_a++; 
			} else {
				// indices equal
				double weight = iter_a->second + factor*iter_b->second; 
				if (weight!=0)
					words.push_back(make_pair(iter_a->first, weight));
				iter_a++;
				iter_b++; 
			}
		}
	}
	while (iter_b!=b.end())
	{
		words.push_back(make_pair(iter_b->first, factor*iter_b->second));
		iter_b++;
	}
	while (iter_a!=a.end())
	{
		words.push_back(*iter_a);
		iter_a++;
	}

	return(SparseVector(words));

}

double sprod_ns(const DenseVector &w, const SparseVector &b, size_t offset)
{
	double ans=0; 
	for (SparseVector::const_iterator iter = b.begin(); iter!=b.end(); ++iter)
	{
		ans += w[iter->first+offset]*iter->second;
	}
	return ans;
}

void multadd_ns(DenseVector &w, const SparseVector &b, double factor, size_t offset)
{
	for (SparseVector::const_iterator iter = b.begin(); iter!=b.end(); ++iter)
	{
		w[iter->first+offset] += factor*iter->second;
	}
}


// friend
ostream& operator<<(ostream &out, const SparseVector &sp)
{
	out << "["; 
	for (SparseVector::const_iterator iter = sp.begin(); iter!=sp.end(); ++iter)
	{
		out << "(" << iter->first << "," << iter->second << "), ";
	}
	out << "]";
	return out;
}
