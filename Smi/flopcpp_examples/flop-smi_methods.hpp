/**
	\brief Scenario-tree classes for the FlopCpp-Smi modelling framework
	\author Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file includes declarations of misc. methods used for the examples,
	such as creating the SMI core models or solving the deterministic equivalent.
**/

#ifndef FLOP_SMI_FUNCTIONS_H
#define FLOP_SMI_FUNCTIONS_H

// Change these two lines to use a different solver
#include "OsiClpSolverInterface.hpp"
#define OSI_SOLVER_INTERFACE OsiClpSolverInterface

#include "SmiScnModel.hpp"
#include "corenode_base.hpp"    // problem-independent core-node base class

namespace FlopSmiEx {

/// creates an \c SmiCoreData object from a vector of core-node models
/** Note that the input must be a vector of the base-class pointers
    \c CoreNodeBase*. If we have a vector of pointers to derived classes, we
    have to make a copy with base-class pointers first. **/
SmiCoreData * create_smi_core (vector<CoreNodeBase *> & coreNodes,
                               std::string problemName = "flop-smi_ex");


/// solves the deterministic equivalent of an \c SmiScnModel model
/**
	\param stochModel the stochastic model to solve
	\param[in] minOrMax optimization direction as \c MP_model::MP_direction
	\param[out] p2solVector will point to a solution vector
	\param[in] problemName used for the output MPS files (in debug mode)
**/
double solve_det_equiv(SmiScnModel & stochModel,
                       MP_model::MP_direction const minOrMax,
                       double const * & p2solVector,
                       std::string const problemName);


/// wrapper for \c solve_det_equiv with \a minOrMax = \c MP_model::MINIMIZE
double minimize_det_equiv(SmiScnModel & stochModel,
                          double const * & p2solVector,
                          std::string const problemName = "flop-smi_ex");


/// wrapper for \c solve_det_equiv with \a minOrMax = \c MP_model::MAXIMIZE
double maximize_det_equiv(SmiScnModel & stochModel,
                          double const * & p2solVector,
                          std::string const problemName = "flop-smi_ex");

} // namespace FlopSmiEx
#endif
