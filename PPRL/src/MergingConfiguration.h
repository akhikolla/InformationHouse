/*
 * MergingConfiguration.h
 *
 *  Created on: 17.05.2017
 *      Author: schnell-42
 */

#ifndef MERGINGCONFIGURATION_H_
#define MERGINGCONFIGURATION_H_

//using namespace std;

#include "MTB_Similarity.h"
#include "MTB_Levenshtein.h"
#include "MTB_DamerauLevenshtein.h"
#include "MTB_Exact.h"
#include "MTB_Jaro.h"
#include "MTB_Ngram.h"
#include "MTB_JaroWinklerMcLaughlinWinklerLynch.h"
#include "MTB_LCS.h"
#include "MTB_Baystat.h"
#include "MTB_Reth.h"
#include "MTB_Soundex.h"
#include "MTB_DoubleMetaphone.h"
#include "MTB_Metaphone.h"
#include "MTB_Tanimoto.h"
#include <string>
#include <Rcpp.h>

/**
* This class is used to set the configuration for one specific Merge set
*/
class MergingConfiguration // extends AbstractVariablesSimilarityTriple
{
private:

  double initialM = 0.9;
  double initialU = 0.1;
  double maxM = 0.9999;
  double minU = 0.0001;
  double epsilon = 0.0004;
  double p = 0.5;
  //bool useValueSpecificM = false;
  //bool useValueSpecificU = false;
  //bool useAutomatchArrayMatching = false;
  //MergeVariable	masterMergeVariable1;
  //MergeVariable	masterMergeVariable2;
  double m = 0.0;
  double u = 0.0;
  double lower = 0.0;
  double upper = 0.0;
  double missingWeight = 0.0;
  int frequenciesSum = 0;
  double jaroWeightFactor = 1.0;
  bool useGlobalM = false;
  bool useGlobalU = false;
  bool ind_c0 = false;
  bool ind_c1 = false;
  float threshold = 1.0; // in case of mbt = minSimilarity
  MTB_Similarity *similarity = NULL;
  string algorithm;
  string blocking = "0";
  float looseThreshold = 0.7;
  float tightThreshold = 0.8;
  int windowSize = 5;
  int lenNgram = 2;

public:
  MergingConfiguration(double m, double u, string algorithm);
  MergingConfiguration();
  ~MergingConfiguration();

  int getFrequenciesSum() {
    return frequenciesSum;
  }
  void setAlgorithm(string algorithm_);

  MTB_Similarity* getAlgorithm();

  string getAlgorithmName() {
    return this->algorithm;
  }

  double getInitialM() {
    return this->initialM;
  }

  double getInitialU() {
    return this->initialU;
  }

  double getM() {
    if (this->m != 0.0) {
      return this->m;
    }
    return getInitialM();
  }

  double getU() {
    if (this->u != 0.0) {
      return this->u;
    }
    return getInitialU();
  }

  double getMaxM() {
    return this->maxM;
  }

  double getMinU() {
    return this->minU;
  }

  double getP() {
    return this->p;
  }

  double getEpsilon() {
    return this->epsilon;
  }

  double getLower() {
    return this->lower;
  }

  double getUpper() {
    return this->upper;
  }

  double getMissingWeight() {
    return this->missingWeight;
  }
  void setM(double m_) {
    this->m = m_;
  }

  void setU(double u_) {
    this->u = u_;
  }

  void setEpsilon(double epsilon_) {
    this->epsilon = epsilon_;
  }

  void setP(double p_) {
    this->p = p_;
  }

  void setLower(double lower_) {
    this->lower = lower_;
  }

  void setUpper(double upper_) {
    this->upper = upper_;
  }

  double getJaroWeightFactor() {
    return this->jaroWeightFactor;
  }

  bool getUseGlobalM() {
    return this->useGlobalM;
  }

  void setUseGlobalM(bool useGlobalM_) {
    this->useGlobalM = useGlobalM_;
  }

  bool getUseGlobalU() {
    return this->useGlobalU;
  }

  void setUseGlobalU(bool useGlobalU_) {
    this->useGlobalU = useGlobalU_;
  }

  string getBlocking() {
    return this->blocking;
  }

  void setBlocking(string blocking_) {
    this->blocking = blocking_;
  }

  void setJaroWeightFactor(double jaroWeightFactor_) {
    this->jaroWeightFactor = jaroWeightFactor_;
  }

  void setInd_c0(bool ind_c0_) {
    this->ind_c0 = ind_c0_;
  }

  void setInd_c1(bool ind_c1_) {
    this->ind_c1 = ind_c1_;
  }

  void setThreshold(float threshold_) {
    this->threshold = threshold_;
  }

  float getThreshold() {
    return this->threshold;
  }

  void setWindowSize(int windowSize_){
    this->windowSize = windowSize_;
  };

  int getWindowSize(){
    return this->windowSize;
  };

  void setLooseThreshold(float looseThreshold_){
    this->looseThreshold = looseThreshold_;
  }

  void setTightThreshold(float tightThreshold_){
    this->tightThreshold = tightThreshold_;
  }

  float getLooseThreshold(){
    return this->looseThreshold;
  }

  float getTightThreshold(){
    return this->tightThreshold;
  }

  void setLenNgram(int lenNgram){
    this->lenNgram = lenNgram;
  }
  //int l1,l2;

  //	public MergeVariable getMasterMergeVariable1()
  //	{
  //		return this.masterMergeVariable1;
  //	}
  //
  //	public MergeVariable getMasterMergeVariable2()
  //	{
  //		return this.masterMergeVariable2;
  //	}

  //	public double getMissingWeight()
  //	{
  //		return this.missingWeight;
  //	}
  //

  //
  //	public boolean isArrayTriple()
  //	{
  //		return (this.mergeVariables1.size() > 1) || (this.mergeVariables2.size() > 1);
  //	}
  //
  //	public boolean isUseAutomatchArrayMatching()
  //	{
  //		return this.useAutomatchArrayMatching;
  //	}
  //
  //	public boolean isUseValueSpecificM()
  //	{
  //		return this.useValueSpecificM;
  //	}
  //
  //	public boolean isUseValueSpecificU()
  //	{
  //		return this.useValueSpecificU;
  //	}
  //
  //	public void setFrequencies(HashMap<string, Integer> frequencies)
  //	{
  //		this.frequenciesSum = 0;
  //		final Iterator<Integer> i = frequencies.values().iterator();
  //		while (i.hasNext())
  //		{
  //			this.frequenciesSum += i.next();
  //		}
  //	}
  //
  //	public void setInitialM(double initialM)
  //	{
  //		if (this.initialM != initialM)
  //		{
  //			final double oldInitialM = this.initialM;
  //			this.initialM = initialM;
  //			this.changes.firePropertyChange("initialM", oldInitialM, initialM);
  //		}
  //	}
  //
  //	public void setInitialU(double initialU)
  //	{
  //		if (this.initialU != initialU)
  //		{
  //			final double oldInitialU = this.initialU;
  //			this.initialU = initialU;
  //			this.changes.firePropertyChange("initialU", oldInitialU, initialU);
  //		}
  //	}
  //
  //	public void setMasterMergeVariable1(MergeVariable masterMergeVariable1)
  //	{
  //		if (this.mergeVariables1.contains(masterMergeVariable1) && (this.masterMergeVariable1 != masterMergeVariable1))
  //		{
  //			final MergeVariable oldMasterMergeVariable1 = this.masterMergeVariable1;
  //			this.masterMergeVariable1 = masterMergeVariable1;
  //			this.changes.firePropertyChange("masterMergeVariable1", oldMasterMergeVariable1, masterMergeVariable1);
  //		}
  //	}
  //
  //	public void setMasterMergeVariable2(MergeVariable masterMergeVariable2)
  //	{
  //		if (this.mergeVariables2.contains(masterMergeVariable2) && (this.masterMergeVariable2 != masterMergeVariable2))
  //		{
  //			final MergeVariable oldMasterMergeVariable2 = this.masterMergeVariable2;
  //			this.masterMergeVariable2 = masterMergeVariable2;
  //			this.changes.firePropertyChange("masterMergeVariable2", oldMasterMergeVariable2, masterMergeVariable2);
  //		}
  //	}
  //
  //	public void setMaxM(double maxM)
  //	{
  //		if (this.maxM != maxM)
  //		{
  //			final double oldMaxM = this.maxM;
  //			this.maxM = maxM;
  //			this.changes.firePropertyChange("maxM", oldMaxM, maxM);
  //		}
  //	}
  //
  //	@Override
  //	public void setMergeVariables1(ArrayList<MergeVariable> mergeVariables1)
  //	{
  //		super.setMergeVariables1(mergeVariables1);
  //		if (mergeVariables1 == null)
  //		{
  //			setMasterMergeVariable1(null);
  //		}
  //		else
  //			if (!this.mergeVariables1.contains(this.masterMergeVariable1))
  //			{
  //				setMasterMergeVariable1(getMergeVar1());
  //			}
  //	}
  //
  //	@Override
  //	public void setMergeVariables2(ArrayList<MergeVariable> mergeVariables2)
  //	{
  //		super.setMergeVariables2(mergeVariables2);
  //		if (mergeVariables2 == null)
  //		{
  //			setMasterMergeVariable2(null);
  //		}
  //		else
  //			if (!this.mergeVariables2.contains(this.masterMergeVariable2))
  //			{
  //				setMasterMergeVariable2(getMergeVar2());
  //			}
  //	}
  //
  //	public void setMinU(double minU)
  //	{
  //		if (this.minU != minU)
  //		{
  //			final double oldMinU = this.minU;
  //			this.minU = minU;
  //			this.changes.firePropertyChange("minU", oldMinU, minU);
  //		}
  //	}
  //
  //	public void setMissingWeight(double missingWeight)
  //	{
  //		if (this.missingWeight != missingWeight)
  //		{
  //			final double oldMissingWeight = this.missingWeight;
  //			this.missingWeight = missingWeight;
  //			this.changes.firePropertyChange("missingWeight", oldMissingWeight, missingWeight);
  //		}
  //	}
  //

  //
  //	public void setUseAutomatchArrayMatching(boolean useAutomatchArrayMatching)
  //	{
  //		this.useAutomatchArrayMatching = useAutomatchArrayMatching;
  //	}
  //
  //	public void setUseValueSpecificM(boolean useValueSpecificM)
  //	{
  //		if (this.useValueSpecificM != useValueSpecificM)
  //		{
  //			final boolean oldUseValueSpecificM = this.useValueSpecificM;
  //			this.useValueSpecificM = useValueSpecificM;
  //			this.changes.firePropertyChange("useValueSpecificM", oldUseValueSpecificM, useValueSpecificM);
  //		}
  //	}
  //
  //	public void setUseValueSpecificU(boolean useValueSpecificU)
  //	{
  //		if (this.useValueSpecificU != useValueSpecificU)
  //		{
  //			final boolean oldUseValueSpecificU = this.useValueSpecificU;
  //			this.useValueSpecificU = useValueSpecificU;
  //			this.changes.firePropertyChange("useValueSpecificU", oldUseValueSpecificU, useValueSpecificU);
  //		}
  //	}
};

#endif /* MERGINGCONFIGURATION_H_ */
