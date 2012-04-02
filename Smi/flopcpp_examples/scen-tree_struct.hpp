/**
	\brief Scenario-tree classes for the FlopCpp-Smi modelling framework
	\author Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file includes declarations for classes describing different
	scenario-tree structures.
**/

#ifndef SCENTREE_STRUCT_HPP
#define SCENTREE_STRUCT_HPP

#include <iostream>
#include <vector>
#include <cassert>

namespace FlopSmiEx {
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

/// Base class for scenario-trees
class ScenTreeStruct {
protected:
	int nmbScens;           ///< number of scenarios
	vector<int> leaves;     ///< list of leaf nodes - they define scenarios

public:
	int nmbStages;          ///< number of stages, again counted from zero

	/// Constructor
	ScenTreeStruct(int const nScens, int const nStages = 0)
	: nmbScens(nScens), leaves(nScens, -1), nmbStages(nStages)
	{}
	/// Destructor
	virtual ~ScenTreeStruct() {}

	/// Get the parent of a given node
	/** In a general case, this would be given by a table, for balanced trees
	    one can use a simple formula. The question is what to do with the root:
	    should the function return 0, -1, or throw an exception? **/
	virtual int get_parent_node(int n) const = 0;

	/// Get the number of stages
	int get_nmb_stages() const { return nmbStages; }

	/// Get the number of scenarios = the number of leaves
	int get_nmb_scens() const { return static_cast<int>(leaves.size()); }

	/// Get probability of a given scenarios
	virtual double get_scen_prob(int const sc) const = 0;

	/// Get the vector of nodes of a given scenarios
	virtual int const * get_scen_nodes(int const sc) = 0;

	/// Get the vector of nodes of a core scenario - just a wrapper
	inline int const * get_core_scen() { return get_scen_nodes(0); }

	/// Get vector of nodes of a next scenario in the list
	/** This assumes that the class itself will keep track of the calls
	    and hence knows what is the next scenario.
	    This is used for Smi, when each scenario is described as a difference
	    from a given parent scenario.
	    \param[out] scen number of the new scenario
	    \param[out] parentScen parent scenario
	    \param[out] branchStage stage where \a scen diverges from \a parentScen
	    \param[out] prob probability of the new scenario
	    \return vector of node indices of the new scen, NULL at the end **/
	virtual int const * get_next_scen(int &scen, int &parentScen,
	                                  int &branchStage, double &prob) = 0;
};


/// Class for balanced binary trees
class BinTreeStruct : public ScenTreeStruct {
protected:
	int * scenNodeNmb; ///< vector of nodes of a scenario - for internal use
	int nextScen;      ///< next scenario to be processed by gen_next_scen
	double scenProb;   ///< probability of each scenario (equiprobable)

	/// Fill \c scenNodeNmb with nodes of a given scenarios
	/// \return first stage when the scen. differs from current \c scenNodeNmb
	int set_scen_nodes(int const sc);

public:
	/// Constructor
	/** A binary tree with T periods has 2^(T-1) scenarios.
	    We number the nodes in a breadth-first manner, so the leaves are at
	    the end, with first leaf at node 2^(T-1)-1 **/
	BinTreeStruct(int const T);

	~BinTreeStruct(){
		delete[] scenNodeNmb;
	}

	int get_parent_node(int n) const {
		return (n-1) / 2;    // This gives: get_parent_node(0) = 0
	}

	double get_scen_prob(int const sc) const { return scenProb; }

	int const * get_scen_nodes(int const sc) {
		set_scen_nodes(sc);
		return scenNodeNmb;
	}

	int const * get_next_scen(int &scen, int &parentScen,
	                          int &branchStage, double &prob);
};


/// Class for a 2-stage tree (bush)
class TwoStageTree : public ScenTreeStruct {
protected:
	int nmbScens;         ///< number of scenarios
	vector<int> leaves;   ///< list of leaf nodes - they define scenarios
	vector<double> probs; ///< scenario probabilities
	int scenNodeNmb[2];   ///< vector of nodes of a scenario - for internal use
	int nextScen;         ///< next scenario to be processed by gen_next_scen

	/// Fill \c scenNodeNmb with nodes of a given scenarios
	/// \return first stage when the scen. differs from current \c scenNodeNmb
	int set_scen_nodes(int const sc) {
		scenNodeNmb[1] = 1 + sc; // scenNodeNmb[0] == 0, the root
		return 1;
	}

public:
	int nmbStages;      ///< number of stages, again counted from zero

	/// Constructor
	TwoStageTree(int const nScens);

	/// Destructor
	 ~TwoStageTree() {}

	/// Get the parent of a given node
	int get_parent_node(int n) const { return 0; }

	/// Set scenario probabilities (if non-equiprobable)
	void set_scen_prob(double * const pr);

	/// Get probability of a given scenarios
	double get_scen_prob(int const sc) const { return probs[sc]; }

	/// Get the vector of nodes of a given scenarios
	int const * get_scen_nodes(int const sc) {
		set_scen_nodes(sc);
		return scenNodeNmb;
	}

	/// Get vector of nodes of a next scenario in the list
	/** This assumes that the class itself will keep track of the calls
	    and hence knows what is the next scenario.
	    This is used for Smi, when each scenario is described as a difference
	    from a given parent scenario.
	    \param[out] scen number of the new scenario
	    \param[out] parentScen parent scenario
	    \param[out] branchStage stage where \a scen diverges from \a parentScen
	    \param[out] prob probability of the new scenario
	    \return vector of node indices of the new scen, NULL at the end **/
	int const * get_next_scen(int &scen, int &parentScen,
	                          int &branchStage, double &prob);
};


} // namespace FlopSmiEx
#endif
