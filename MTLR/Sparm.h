#pragma once

#include <string>
#include <vector>

class Sparm
{
private:
	std::string m_inputFilename; 
	std::string m_modelFilename; 
	int m_maxMonth;
	size_t m_sizePsi;          // feature dimension [by stacking different weight vectors for each month]
	size_t m_original_sizePsi; // original feature dimension
	double m_C1; 
	double m_C2; 
	std::string m_lossType; 
	bool m_printProb;          // print probabilities during prediction (for plotting survival cdf)
	double m_threshold;        // threshold for classification
	int m_bundleSize; 
	int m_trainUncensored; 
	std::string m_initialWeightFile; 
	std::vector<double> m_quantTime; 
	int m_regularizer; 
	std::vector<double> m_xqueries; // query point for finding out survival probability for a given time  (used in testing only)
	std::vector<double> m_yqueries; // query point for finding out time elapsed given survival probability (used in testing only)
public:
	Sparm(void);
	~Sparm(void);

	void ReadParam(int argc, char ** argv); 

	std::string GetInputFile() const
	{
		return m_inputFilename; 
	}
	std::string GetModelFile() const
	{
		return m_modelFilename; 
	}
	int GetMaxMonth() const
	{
		return m_maxMonth; 
	}
	size_t GetSizePsi() const
	{
		return m_sizePsi;
	}
	size_t GetOriginalSizePsi() const
	{
		return m_original_sizePsi; 
	}
	void SetSizePsi(size_t n)
	{
		m_original_sizePsi = n; 
		m_sizePsi = m_maxMonth*n + m_maxMonth; 
	}
	double GetC1() const
	{
		return m_C1;
	}
	double GetC2() const
	{
		return m_C2;
	}
	std::string GetLossType() const
	{
		return m_lossType; 
	}
	bool GetPrintProb() const 
	{
		return m_printProb; 
	}
	double GetThreshold() const
	{
		return m_threshold;
	}
	int GetBundleSize() const
	{
		return m_bundleSize;
	}
	int GetTrainUncensored() const
	{
		return m_trainUncensored;
	}
	std::string GetInitialWeightFile() const
	{
		return m_initialWeightFile; 
	}
	const std::vector<double> &GetQuantTime() const
	{
		return m_quantTime; 
	}
	// may 9: new default for q-file
	void SetQuantTime(std::vector<double> &q)
	{
		m_quantTime = q; 
		m_maxMonth = q.size(); 
	}
	// end may 9
	void ReadQuantTime(std::string filename); 
	int GetRegularizer() const
	{
		return m_regularizer; 
	}
	// aug 2
	void SetRegularizer(int r)
	{
		m_regularizer = r; 
	}
	// end aug 2
	// june 25
	const std::vector<double> &GetXQueries() const
	{
		return m_xqueries;
	}
	const std::vector<double> &GetYQueries() const
	{
		return m_yqueries;
	}
};

