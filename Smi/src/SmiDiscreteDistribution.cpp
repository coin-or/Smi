#include "SmiDiscreteDistribution.hpp"

SmiDiscreteDistribution::~SmiDiscreteDistribution() {
		for (size_t i=0; i<smiDiscrete_.size(); ++i)
			delete smiDiscrete_[i];
	}

void SmiDiscreteRV::addEvent(CoinPackedMatrix &matrix,
				CoinPackedVector &dclo, CoinPackedVector &dcup,
				CoinPackedVector &dobj,
				CoinPackedVector &drlo, CoinPackedVector &drup, double prob)
	{
		SmiLinearData d(matrix,dclo,dcup,dobj,drlo,drup);
		SmiDiscreteEvent *e = new SmiDiscreteEvent(d,prob);
		this->events_.push_back(e); 
		prob_+=prob;
	}