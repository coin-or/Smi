/*
	Core-node base class for the FlopCpp-Smi modelling framework
	authors: Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file includes implementations of methods from the general,
	problem independent, core-node base class.
*/

#include "corenode_base.hpp"

#include <vector>

using namespace FlopSmiEx;

// Initialize the static members - this has to be done outside the class
CoreNodeBase * CoreNodeBase::p2activeNode = NULL;
double const * CoreNodeBase::p2varValues = NULL;

void CoreNodeBase::make_obj_func_rec() {
	// recursive sum of the obj. func from the node and its child
	objFuncRec = objFuncNode;
	if (p2child) {
		p2child->make_obj_func_rec();
		objFuncRec = objFuncRec + p2child->objFuncRec;
	}
}
