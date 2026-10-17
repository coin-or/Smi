/*
	Scenario-tree classes for the FlopCpp-Smi modelling framework
	authors: Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file includes implementations of methods from classes describing
	different scenario-tree structures.
*/

#include <cstdlib>
#include <cmath>
#include "scen-tree_struct.hpp"

using namespace std;
using namespace FlopSmiEx;

// ----------------------------------------------------------------------------
// class BinTreeStruct

int BinTreeStruct::set_scen_nodes(int const sc)
{
	int n = leaves[sc];
	int t = nmbStages-1;
	// start from the leaf, stop when the scenario node is the same
	// as the node in scenNodeNmb[t] - no need to continue then..
	while (n != scenNodeNmb[t]) {
		assert (n > 0 && t > 0 && "We should have scenNodeNmb[0] = 0!");
		scenNodeNmb[t] = n;
		n = get_parent_node(n);
		t--;
	}
	return t + 1; // t points to the last common stage
}


BinTreeStruct::BinTreeStruct(int const T)
: ScenTreeStruct((int) floor(pow(2.0, T-1) + 0.5), T),
  nextScen(0), scenProb(1.0 / nmbScens)
{
	int firstLeaf = (int) floor(pow(2.0, T-1) + 0.5) - 1;
	assert ((int) leaves.size() == nmbScens && "done in the parent class");
	for (int s = 0; s < nmbScens; s++) {
		leaves[s] = firstLeaf + s;
	}
	scenNodeNmb = new (nothrow) int[T];
	if (scenNodeNmb == NULL) {
		cerr << "Memory allocation failed!" << endl; exit(1);
	}
	// initialize the array with the first scenario
	scenNodeNmb[T - 1] = firstLeaf;
	for (int t = T - 1; t > 0; --t)
		scenNodeNmb[t-1] = get_parent_node(scenNodeNmb[t]);
	assert (scenNodeNmb[0] == 0 && "Tree must start in the root");
}


int const * BinTreeStruct::get_next_scen(int &scen, int &parentScen,
                                         int &branchStage, double &prob)
{
	if (nextScen == nmbScens) {
		return NULL; // no more scenarios
	}
	scen = nextScen;
	parentScen = (scen == 0 ? 0 : scen - 1);
	if (scenNodeNmb[nmbStages - 1] != leaves[parentScen]) {
		cerr << "Warning: scenNodeNmb does not hold the prev. scenario!" << endl;
		set_scen_nodes(parentScen);
	}
	branchStage = (scen == 0 ? 1 : set_scen_nodes(scen));
	prob = scenProb;
	nextScen++;
	return scenNodeNmb;
}


// ----------------------------------------------------------------------------
// class TwoStageTree

TwoStageTree::TwoStageTree(int const nScens)
: ScenTreeStruct(nScens, 2),
  nmbScens(nScens), leaves(nScens, -1.0),
  probs(nScens, 1.0 / static_cast<double>(nScens)),
  nextScen(0), nmbStages(2)
{
	scenNodeNmb[0] = 0;
	for (int s = 0; s < nmbScens; s++) {
		leaves[s] = (1 + s); // root = 0, first leaf = 1
	}
}

void TwoStageTree::set_scen_prob(double * const pr)
{
	//assert (pr.length() == nmbScens && "dimenstion check");
	for (int s = 0; s < nmbScens; s++) {
		probs[s] = pr[s];
	}
}

int const * TwoStageTree::get_next_scen(int &scen, int &parentScen,
	                                     int &branchStage, double &prob)
{
	if (nextScen == nmbScens) {
		return NULL; // no more scenarios
	}
	scen = nextScen;
	parentScen = (scen == 0 ? 0 : scen - 1);
	branchStage = (scen == 0 ? 1 : set_scen_nodes(scen));
	prob = probs[scen];
	nextScen++;
	return scenNodeNmb;
}
