/*
 * MTB_ProbabilityCalculation.h
 *
 *  Created on: 16.05.2017
 *      Author: schnell-42
 */

#ifndef MTB_PROBABILITYCALCULATION_H_
#define MTB_PROBABILITYCALCULATION_H_
#include <math.h>

#include "MTB_EMAlgorithm.h"
#include "MergingConfiguration.h"

#include "Util.h"

using namespace std;

class MTB_ProbabilityCalculation
{

	/**
	 * @param mtbs
	 * @param name
	 * @param description
	 */
public:
	MTB_ProbabilityCalculation();
	~MTB_ProbabilityCalculation();
	/**
	 * counts frequencies
	 *
	 * @param mergingConfigurationList
	 * @param blocks
	 * @return returns a vector frequencyCount with two values, every value of the data1 is compared to every value in data2,
	 *  if they are equal, frequencyCount[1] is increased, else  frequencyCount[0] is increased
	 */
//	vector<int> countFrequencies(
//			vector<MergingConfiguration> mergingConfigurationList,
//			vector<MTB_ObservationSetPair> blocks);

	vector<int> countFrequencies(
			vector<string> data1, vector<string> data2);

	string getDescription() {
		return "EM algorithm...";
	}

	string getName() {
		return this->NAME;
	}

	vector<MergingConfiguration> run(vector<MergingConfiguration> mc,
			vector<string> data1, vector<string> id1, vector<string> data2,
			vector<string> id2, MTB_EMAlgorithm em);

	vector<double> getUArray(vector<int> frequencies, int compFreq);

	vector<MergingConfiguration> setMArray(vector<double> mArray,
			vector<MergingConfiguration> mergingConfigurationList);

	vector<MergingConfiguration> setUArray(vector<double> mArray,
				vector<MergingConfiguration> mergingConfigurationList);


//
//	@Override
//	protected void onFinish()
//	{
//		final ArrayList<MergingConfiguration> vatList = this.getMtbs().getMergingConfigurationList().getMergingConfigurationList();
//		if (this.em != null)
//		{
//			this.reportModule.add(new KeyValuePair(KEY_EM_ITERATIONS, this.em.getIterations()));
//		}
//		for (int i = 0; i < vatList.size(); i++)
//		{
//			final double m = vatList.get(i).getM();
//			final double u = vatList.get(i).getU();
//			final double log_m_div_u = Math.log(m / u) / Math.log(2);
//			final double log_1_minus_m_div_1_minus_u = Math.log((1 - m) / (1 - u)) / Math.log(2);
//			final string m_u_values = " m: " + m + ", u: " + u + ", log(m/u): " + log_m_div_u + ", log((1-m)/(1-u)): " + log_1_minus_m_div_1_minus_u;
//			if (this.mtbs.getDeduplicate().getValue())
//			{
//				this.reportModule.add(new KeyValuePair(KEY_M_U_PROBS, vatList.get(i).getVar1().getName() + m_u_values));
//			}
//			else
//			{
//				this.reportModule.add(new KeyValuePair(KEY_M_U_PROBS, vatList.get(i).getVar1().getName() + ", " + vatList.get(i).getVar2().getName()
//						+ m_u_values));
//			}
//		}
//	}

private:

	string NAME = "Probabilities";
	string KEY_EM_ITERATIONS = "EM Iterations";
	string KEY_EM_FAILED = "EM estimation failure";
	string KEY_M_U_PROBS = "m/u";
	string VALUE_TOO_FEW_OBS =
			"Too few observation pairs for chosen number of matching variables";
	string VALUE_U_NAN = "Some u value is not a number";
	string VALUE_U_1 = "Some u value is 1";
	string error = "";
	int value = 0;

//
//	private double[] getInitialMArray(ArrayList<MergingConfiguration> mergingConfigurationList)
//	{
//		double[] mArray = null;
//		if (mergingConfigurationList != null)
//		{
//			final int vstsSize = mergingConfigurationList.size();
//			mArray = new double[vstsSize];
//			for (int i = 0; i < vstsSize; i++)
//			{
//				mArray[i] = mergingConfigurationList.get(i).getInitialM();
//			}
//		}
//		return mArray;
//	}
//
//	private double[] getInitialUArray(ArrayList<MergingConfiguration> mergingConfigurationsList)
//	{
//		double[] uArray = null;
//		if (mergingConfigurationsList != null)
//		{
//			final int vstsSize = mergingConfigurationsList.size();
//			uArray = new double[vstsSize];
//			for (int i = 0; i < vstsSize; i++)
//			{
//				uArray[i] = mergingConfigurationsList.get(i).getInitialU();
//			}
//		}
//		return uArray;
//	}
//


};

#endif /* MTB_PROBABILITYCALCULATION_H_ */
