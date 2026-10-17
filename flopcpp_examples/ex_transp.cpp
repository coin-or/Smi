/*
	Scenario-tree classes for the FlopCpp-Smi modelling framework
	authors: Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file implements a production and transportation example from
	"Benders decomposition for stochastic programming with GAMS" by
	Erwin Kalvelagen; http://www.amsterdamoptimization.com/pdf/stochbenders.pdf
	It is the same example as in stochbenders.cpp from the FlopC++ package.
**/

#include <iomanip>              // needed for setw() in cout

#include "corenode_base.hpp"    // problem-independent core-node base class
#include "scen-tree_struct.hpp" // classes for scenario-tree structures
#include "flop-smi_methods.hpp" // misc. methods

using namespace FlopSmiEx;


// -------------------------------------------------------------------------- //
//           Classes and Definitions for the Transportation Model             //
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Core-nodes base class

/// This is the base class for all stage models.
class CoreNodeTr : public CoreNodeBase {
public:
	/// constructor
	CoreNodeTr(CoreNodeTr * p2pred, int const nmbFact, int const nmbDistC,
	           int t = -1)
	: CoreNodeBase(p2pred, t), FACT(nmbFact), DEST(nmbDistC)
	{}
	virtual ~CoreNodeTr() {}     ///< destructor

	MP_set
		FACT,       ///< set of factories
		DEST;       ///< set of destinations (distribution centers)
	MP_index i, j; ///< indices used in formulas (i for factories)
};


// -------------------------------------------------------------------------- //
// Core-nodes classes used in the code

/// Class for the root node of the tree
class RootNodeTr : public CoreNodeTr {
protected:
	MP_data prodCap; ///< production capacities
	MP_data transpC; ///< transportation costs

	double prodCost; ///< production costs

public:
	RootNodeTr(int const nmbFact, int const nmbDest, double const prCost,
	  double * const prCap, double const ** const trC)
	: CoreNodeTr(NULL, nmbFact, nmbDest, 0),
	  prodCap(prCap, FACT), transpC(FACT, DEST),
	  prodCost(prCost),
	  prod(FACT), ship(FACT, DEST),
	  volumeBalance(FACT)
	{
		volumeBalance(i) = sum(DEST(j), ship(i,j)) == prod(i);

		// objective - only costs & maximizing -> negative
		objFuncNode = -1 * sum(FACT(i), prodCost * prod(i)
		                       + sum(DEST(j), transpC(i,j) * ship(i,j)));

		for (int f = 0; f < nmbFact; ++f) {
			prod.upperLimit(f) = prodCap(f);
			for (int d = 0; d < nmbDest; ++d) {
				transpC(f,d) = trC[f][d];
			}
		}
	}

	SP_variable prod; ///< production per factory
	SP_variable ship; ///< shipping between factories and dist. centers

	/// A product-volume balancing constraint
	SP_constraint volumeBalance;

	/// reporting of the results
	void report_results() {
		printf("Optimal production:");
		for (int f = 0; f < FACT.size(); ++f) {
			printf("  %3g", prod.value(f));
		}
		printf("\nOptimal shipping:\n");
		for (int f = 0; f < FACT.size(); ++f) {
			for (int d = 0; d < DEST.size(); ++d) {
				if (ship.value(f, d) > 1e-3) {
					printf (" %1d -> %1d .. %3g\n", f, d, ship.value(f, d));
				}
			}
		}
	}
};


/// This is the class for the leaves, i.e. the last-stage nodes
class LeafNodeTr : public CoreNodeTr {
protected:
	SP_variable sales;
	SP_variable waste;

	double price;
	double wasteC;

public:
	/// Constructs a \a LeafNode object
	LeafNodeTr(RootNodeTr *p2pred, double const pr, double const wC)
	: CoreNodeTr(p2pred, p2pred->FACT.size(), p2pred->DEST.size(), 1),
	  sales(DEST), waste(DEST), price(pr), wasteC(wC),
	  volumeBalance(DEST)
	{
		volumeBalance(j) = sum(FACT(i), p2pred->ship(i,j)) == sales(j) + waste(j);

		objFuncNode = sum(DEST(j), price * sales(j) - wasteC * waste(j));
	}

	/// A product-volume balancing constraint
	SP_constraint volumeBalance;

	/// Create a vector with elements that differ between scenarios
	/** In this case, the only difference between the scenarios is the
	    demand, which forms the upper bounds on the sales variables **/
	void make_diff_vect(CoinPackedVector &bDiff, double const scDemand[])
	{
		bDiff.clear();
		for (int j = 0; j < DEST.size(); ++j) {
			bDiff.insert(sales(j).getColumn(), scDemand[j]);
		}
	}
};


// -------------------------------------------------------------------------- //
//                               The Code                                     //
// -------------------------------------------------------------------------- //

int main()
{
	std::string const problemName = "transp-ex";

	// ----------------------------------------------------------------------- //
	// DATA - would be normally read from a file

	// set sizes etc
	int const nmbStages = 2;
	int const nmbScens  = 3;
	int const nmbFact   = 3;
	int const nmbDest   = 5;

	double shipC[nmbFact][nmbDest] = {
		{2.49, 5.21, 3.76, 4.85, 2.07},
		{1.46, 2.54, 1.84, 1.86, 4.76},
		{3.26, 3.08, 2.60, 3.76, 4.45}
	};
	double const * p2shipC[nmbFact];
	for (int i = 0; i < nmbFact; ++i)
		p2shipC[i] = shipC[i];

	double const demand[nmbScens][nmbDest] = {
		{150, 100, 250, 300, 600},
		{160, 120, 270, 325, 700},
		{170, 135, 300, 350, 800}
	};

	double prodCap[nmbFact] = {500, 450, 650};

	double scProb[nmbScens] = {0.25, 0.50, 0.25};
	double const prodCost  = 14;
	double const price     = 24;
	double const wasteCost = 4;

	// init the scenario tree
	TwoStageTree scTree(nmbScens);
	scTree.set_scen_prob(scProb);

	assert (scTree.get_nmb_stages() == nmbStages
	        && "Checking that get_nmb_stages() returns what it should.");


	// ----------------------------------------------------------------------- //
	// CREATE THE CORE OBJECT

	// Create scenario tree for the core model, using data for the 1st scenario
	vector<CoreNodeBase *> coreNodes(nmbStages);
	coreNodes[0] = new RootNodeTr(nmbFact, nmbDest, prodCost, prodCap, p2shipC);
	RootNodeTr * p2root = dynamic_cast<RootNodeTr *>(coreNodes[0]);
	coreNodes[1] = new LeafNodeTr(p2root, price, wasteCost);
	LeafNodeTr * p2leaf = dynamic_cast<LeafNodeTr *>(coreNodes[1]);

	SmiCoreData * p2stochCore = create_smi_core(coreNodes, problemName);


	// ----------------------------------------------------------------------- //
	// START BUILDING THE STOCHASTIC MODEL

	// This is done in an SMPS-like fashion, ie. each scenario has a parent
	// scenario it branches from. We then have to specify the branching stage
	// and all the data that are different than the parent's.
	// In our case, the only difference is in the matrix A. We only need to
	// specify the elements that differ from the parent, that is the returns.
	SmiScnModel stochModel;

	// the vector of differences wrt. the previous (parent) scenario
	CoinPackedVector bDiff;

	// Get the node numbers for the Core scenario
	// No allocation here, as this points to the vector in BinTreeStruct!
	int const * scenNodeNmb = scTree.get_core_scen();

	// Add scenarios, one by one.
	int scen, parentScen, branchStage;
	double scenProb;
	// get_next_scen() returns the nodes of the next scen, or NULL at the end
	while ((scenNodeNmb = scTree.get_next_scen(scen, parentScen, branchStage,
	                                           scenProb)) != NULL) {
		// load modified data (upper bounds) into bDiff
		p2leaf->make_diff_vect(bDiff, demand[scen]);

		stochModel.generateScenario(p2stochCore, NULL, NULL, &bDiff,
		                            NULL, NULL, NULL,
		                            branchStage, parentScen, scenProb);
	}


	// ----------------------------------------------------------------------- //
	// SOLVING THE PROBLEM
	/* Since there is no stochastic solver in COIN-OR yet, so we have to solve
	   the model using a deterministic solver on the determ. equivalent. */

	double const * p2OsiSolution = NULL; // solution vector
	double maxDetEqObj = maximize_det_equiv(stochModel, p2OsiSolution,
	                                        problemName);


	// ----------------------------------------------------------------------- //
	// REPORTING
	cout << endl << "Optimal expected profit: " << maxDetEqObj << endl;

	// Have to point the root to the solution vector
	p2root->set_var_values(p2OsiSolution);
	p2root->report_results(); // report using the method from the root node

	// cleaning
	for (int t = 0; t < (int) coreNodes.size(); t++) {
		delete coreNodes[t];
	}
	delete p2stochCore;


	return 0;
}
