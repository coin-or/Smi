/*
	Scenario-tree classes for the FlopCpp-Smi modelling framework
	authors: Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file implements the investment example from "Introduction to
	Stochastic Programming", Birge and Louveaux 1997, pp. 20-28.
*/

#include <iomanip>              // needed for setw() in cout

#include "corenode_base.hpp"    // problem-independent core-node base class
#include "scen-tree_struct.hpp" // classes for scenario-tree structures
#include "flop-smi_methods.hpp" // misc. methods

using namespace FlopSmiEx;


// -------------------------------------------------------------------------- //
//             Classes and Definitions for the Investment Model               //
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Core-nodes base classes

/// This is the base class for all stage models.
class CoreNodeInv : public CoreNodeBase
{
public:
	/// constructor
	CoreNodeInv(CoreNodeInv * p2pred, const int nmbAssets, const int t = -1)
	: CoreNodeBase(p2pred, t), ASSETS(nmbAssets)
	{}
	virtual ~CoreNodeInv() {}     ///< destructor

	MP_set ASSETS;             ///< set of assets
	MP_index a;                ///< index used in formulas

	/// A cash-balancing constraint
	/** There is a cash-balancing constraint at each stage, even if they have
	    different meaning in root (initial wealth), leaves (counting deviation
	    from the target) and the rest of the nodes (cash-flow balance).
	    The advantage of declaring it here is that we get easily access all of
	    them using the same name. **/
	SP_constraint cashBalance;

	/// Function that reports the wealth in the node - just for illustration
	virtual double get_wealth() = 0;

	/// This creates the matrix with changes of the LP matrix for a given scen.
	virtual void make_diff_mat(CoinPackedMatrix &ADiff, double *retData) {}
};


/// Base class for all non-leaf nodes - they have the "buy" variables
/**
	We define this class as all non-leaf classes have variable \c x.
	In addition, the parent pointer is re-cast to NonLeafNodeInv to get access
	to \c x().
**/
class NonLeafNodeInv : public CoreNodeInv {
public:
	SP_variable x; ///< the "buy" variable, defined on ASSETS

	/// constructor
	NonLeafNodeInv(CoreNodeInv * p2pred, const int nmbAssets, const int t = -1)
	: CoreNodeInv(p2pred, nmbAssets, t), x(ASSETS)
	{
		// non-leaf nodes do no contribute to the objective
		objFuncNode = 0.0 * x(0); // Is there a better way to represent zero??
	}
	virtual ~NonLeafNodeInv() {} ///< destructor

	/// get the parent pt. re-casted to \c NonLeafNodeInv* (for balance constr.)
	NonLeafNodeInv * get_parent() {
		return static_cast<NonLeafNodeInv *>(p2parent);
	}

	/// wealth at these node is the sum of the positions
	virtual double get_wealth() {
		double w = 0.0;
		for (int a = 0; a < ASSETS.size(); ++a) {
			w += x.value(a); // remember to call set_var_values() first!
		}
		return w;
	}

	/// return the value of the \c x variables
	void get_x_values(double * solX) {
		for (int a = 0; a < ASSETS.size(); ++a) {
			solX[a] = x.value(a); // remember to call set_var_values() first!
		}
	}
};


// -------------------------------------------------------------------------- //
// Core-nodes classes used in the code

/// Class for the root node of the tree
class RootNodeInv : public NonLeafNodeInv {
protected:
	double initWealth;
public:
	RootNodeInv(const int nmbAssets, const double initCap)
	: NonLeafNodeInv(NULL, nmbAssets, 0), initWealth(initCap)
	{
		cashBalance() = sum(ASSETS, x(ASSETS)) == initWealth;
	}

	/// wealth at the root is simply the initial wealth
	/** This is not necessary, since \c NonLeafNodeInv.get_wealth() works as
	    well - but this is faster. **/
	double get_wealth() { return initWealth; }
};


/// Class for the middle-nodes, i.e. all the nodes between the root and leaves
/** This class share some members with the \c LeafNodeInv class, which could be
    eliminated by creating an addition \c NonRootNode class and derive both
    \c MidStageNode and \c LeafNodeInv from that one.
    This would, however, mean that we would have to deal with multiple
    inheritance (of the 'diamond type') for \c MidStageNode, something we
    wanted to avoid in this example. **/
class MidStageNode : public NonLeafNodeInv {
public:
	MP_data Return;                ///< returns of the assets at this node

	/// Constructs a \c MidStageNode object
	/** \warning We use a <em>shallow copy</em> for the \c MP_data \a Return,
	    i.e. the return values in the constraints will be linked to the array
	    \a retVect. This means that if the external array changes before we
	    build the OSI object (using the \c attach method), the constraints will
	    be changed as well - and if the external array is deallocated, the
	    program will crash on calling the \c attach method! **/
	MidStageNode(CoreNodeInv *p2pred, double *p2retVect, const int t = -1)
	: NonLeafNodeInv(p2pred, p2pred->ASSETS.size(), t),
	  Return(p2retVect, ASSETS)
	{
		// This shows the use of MP_index in a formula
		cashBalance() = sum(ASSETS(a), get_parent()->x(a) * Return(a))
		                == sum(ASSETS(a), x(a));
	}

	/// Create a matrix with elements that differ between scenarios
	/** In this case, the only difference between the scenarios are the
	    asset returns, which appear in the \c cashBalance constraints
	    as coefficients of the \c x variables from the \em parent node.
	    To construct the matrix, we thus have to find the indices (row and
	    column numbers) of the relevant constraints and variables. **/
	void make_diff_mat(CoinPackedMatrix &ADiff, double *retData)
	{
		int i = cashBalance.row_number();
		for (int a = 0; a<ASSETS.size(); a++) {
			int j = get_parent()->x(a).getColumn();
			ADiff.modifyCoefficient(i, j, retData[a]);
		}
	}
};


/// This is the class for the leaves, i.e. the last-stage nodes
class LeafNodeInv : public CoreNodeInv {
public:
	SP_variable w;              ///< deficit/shortage variable
	SP_variable y;              ///< excess/surplus variable
	MP_data Return;             ///< returns of the assets at this node
	double capTarget;           ///< the capital target parameter

	/// Constructs a \a LeafNodeInv object
	/** Unlike \c MidStageNode, this class uses a <em>deep copy</em> for the
	    \c MP_data \a Return, i.e. the return values in the constraints are
	    copied from the \a retVect array to the constraints. This means that
	    the \a retVect array can be safely changed or deleted afterwards. **/
	LeafNodeInv(CoreNodeInv *p2pred, const double *p2retVect, const double capTg,
	         const int t = -1)
	: CoreNodeInv(p2pred, p2pred->ASSETS.size(), t),
	  Return(ASSETS), capTarget(capTg)
	{
		Return.value(p2retVect); // this copis values from retVect to Return
		// This shows a formula without using an additional MP_index
		cashBalance() = sum(ASSETS, get_parent()->x(ASSETS) * Return(ASSETS))
		                + w() - y() == capTarget;

		// we minimize -> objective = 4.0 * shortage - 1.0 * surplus
		objFuncNode = 4.0 * w() - 1.0 * y();
	}

	/// wealth at the leaf is the target +/- the deviations
	/** This is not necessary, since \c NonLeafNodeInv.get_wealth() works as
	    well - but this should be faster. **/
	double get_wealth() {
		return capTarget -w.value() + y.value();
	}

protected:
	/// pointer to the parent, re-casted into \c NonLeafNodeInv* to get \c x()
	NonLeafNodeInv * get_parent() {
		return dynamic_cast<NonLeafNodeInv *>(p2parent);
	}

	/// Create a matrix with elements that differ between scenarios
	/** For more info, see \c MidStageNode.make_diff_mat(). **/
	void make_diff_mat(CoinPackedMatrix &ADiff, double *retData)
	{
		int i = cashBalance.row_number();
		for (int a=0; a<ASSETS.size(); a++) {
			int j = this->get_parent()->x(a).getColumn();
			ADiff.modifyCoefficient(i, j, retData[a]);
		}
	}
};


// -------------------------------------------------------------------------- //
//                               The Code                                     //
// -------------------------------------------------------------------------- //

int main()
{
	std::string const problemName = "invest-ex";
	int t;

	// ----------------------------------------------------------------------- //
	// DATA - would be normally read from a file

	// binary scenario tree with 4 stages: 15 nodes, firstLeaf = 7
	const int nmbStages = 4;
	BinTreeStruct scTree(nmbStages);

	// model parameters
	enum {stocks, bonds, nmbAssets}; // assets to invest to; sets nmbAssets = 2
	double initBudget = 55;          // initial budget
	double capTarget = 80;           // capital target

	// vector of returns at the 14 non-root nodes
	double retData[][nmbAssets] = {
		{1.25, 1.14},
		{1.06, 1.12},
		{1.25, 1.14},
		{1.06, 1.12},
		{1.25, 1.14},
		{1.06, 1.12},
		{1.25, 1.14},
		{1.06, 1.12},
		{1.25, 1.14},
		{1.06, 1.12},
		{1.25, 1.14},
		{1.06, 1.12},
		{1.25, 1.14},
		{1.06, 1.12},
	};

	assert (nmbStages == scTree.get_nmb_stages()
	        && "Checking that get_nmb_stages() returns what it should.");


	// ----------------------------------------------------------------------- //
	// CREATE THE CORE OBJECT

	// Get the node numbers for the Core scenario
	// No allocation here, as this points to the vector in BinTreeStruct!
	int const * scenNodeNmb = scTree.get_core_scen();

	// Create scenario tree for the core model, using data for the 1st scenario
	vector<CoreNodeInv *> coreNodes(nmbStages);
	coreNodes[0] = new RootNodeInv(nmbAssets, initBudget);
	for (t = 1; t < nmbStages-1; ++t) {
		coreNodes[t] = new MidStageNode(coreNodes[t-1],
		                                retData[scenNodeNmb[t]-1], t);
	}
	assert (t == nmbStages-1 && "t should be nmbStages-1 after the loop");
	coreNodes[t] = new LeafNodeInv(coreNodes[t-1], retData[scenNodeNmb[t]-1],
	                            capTarget, t);

	// create a "shortcut object" for the root
	RootNodeInv * p2root = dynamic_cast<RootNodeInv *>(coreNodes[0]);

	/* create_smi_core is a general method, so it needs vector<CoreNodeBase *>;
	   Since it is not possible to directly cast from vector<CoreNodeInv *>,
	   we need to create a temp. vector of base-class pointers */
	vector<CoreNodeBase *> baseCoreNodes(coreNodes.begin(), coreNodes.end());
	SmiCoreData * p2stochCore = create_smi_core(baseCoreNodes, problemName);


	// ----------------------------------------------------------------------- //
	// START BUILDING THE STOCHASTIC MODEL

	/* This is done in an SMPS-like fashion, ie. each scenario has a parent
	   scenario it branches from. We then have to specify the branching stage
	   and all the data that are different than the parent's.
		This means that the structure of this part is always the same, but the
		differences (what differs and what data are needed) are case dependent.

	   In our case, the only difference is in the matrix A, caused by the
	   different values of returns. */
	SmiScnModel stochModel;

	// the sparse matrix of differences wrt. the previous (parent) scenario
	CoinPackedMatrix ADiff;
	// The default constructor creates a column-ordered matrix, while Smi uses
	// row-ordering; it would be done automatically later, but this is faster...
	ADiff.reverseOrdering();

	// Add scenarios, one by one.
	int scen, parentScen, branchStage;
	double scenProb;
	// get_next_scen() returns the nodes of the next scen, or NULL at the end
	while ((scenNodeNmb = scTree.get_next_scen(scen, parentScen, branchStage,
	                                           scenProb)) != NULL) {
		#ifndef NDEBUG
			cout << "Nodes in scen " << scen + 1
			     << " that differ from the prev. scen: ";
			for (t=branchStage; t<nmbStages; t++) {
				cout << setw(2) << scenNodeNmb[t] << " ";
			}
			cout << endl;
		#endif

		// load modified data into matrix ADiff
		ADiff.clear(); // clean the matrix of differences - must reset dimensions!
		ADiff.setDimensions(p2stochCore->getNumRows(), p2stochCore->getNumCols());
		for (t=branchStage; t<nmbStages; t++){
			coreNodes[t]->make_diff_mat(ADiff,retData[scenNodeNmb[t]-1]);
		}

		stochModel.generateScenario(p2stochCore, &ADiff,
		                            NULL, NULL, NULL, NULL, NULL,
		                            branchStage, parentScen, scenProb);
	}

	cout << endl << "The stochastic model has " << stochModel.getNumScenarios()
	     << " scenarios." << endl;
	assert (stochModel.getNumScenarios() == scTree.get_nmb_scens()
	        && "Check number of scens.");


	// ----------------------------------------------------------------------- //
	// SOLVING THE PROBLEM
	/* Since there is no stochastic solver in COIN-OR yet, so we have to solve
	   the model using a deterministic solver on the determ. equivalent. */

	double const * p2OsiSolution = NULL; // solution vector
	#ifndef NDEBUG
		double minDetEqObj = minimize_det_equiv(stochModel, p2OsiSolution,
		                                        problemName);
	#else
		// the returned obj. value is currently used only for debugging
		minimize_det_equiv(stochModel, p2OsiSolution, problemName);
	#endif


	// ----------------------------------------------------------------------- //
	// REPORTING
	/* Even if we use a deterministic solver, we can still get information about
	   the solution on the scenario tree, using the SMI (stochastic) model.
	   This part is completely case dependent. */

	// Report the wealth at each node of the tree, plus the obj. value
	vector<double> nodeWealth(nmbStages, 0);
	double objValue = 0.0;

	// Compute the wealth at each node, by traversing the tree from leafs up
	for (SmiScenarioIndex sc = 0; sc < scTree.get_nmb_scens(); sc ++) {

		// Get the leaf node of scenario sc:
		SmiScnNode * p2node = stochModel.getLeafNode(sc);
		int nodeStage = nmbStages;

		double scProb = p2node->getModelProb(); // probability of the leaf
		double scenObjVal = stochModel.getObjectiveValue(sc); // objective value
		objValue += scProb * scenObjVal;
		printf ("scen %d: prob = %.3f  obj =%7.2f", sc+1, scProb, scenObjVal);

		// This loop traverses the tree, from the leaf to the root
		while (p2node != NULL) {
			assert ((unsigned) p2node->getNumCols()
			        == (coreNodes[nodeStage-1]->all_variables.size())
			        && "The OSI and CoreNodeInv models have the same #variables.");

			// point the CoreNodeInv model to the correct solution vector
			int firstColInOsi = p2node->getColStart();
			coreNodes[nodeStage-1]->set_var_values(&(p2OsiSolution[firstColInOsi]));

			// get the wealth, using nodes of the core model
			nodeWealth[nodeStage-1] = coreNodes[nodeStage-1]->get_wealth();

			// get the parent node (Root will return NULL, stopping the loop)
			p2node = p2node->getParent();
			nodeStage--;
		}

		printf (";  wealth:");
		for (int t = 0; t<nmbStages-1; t++)
			printf ("%6.2f ->", nodeWealth[t]);
		printf ("%6.2f\n", nodeWealth[nmbStages-1]);
	}
	printf ("%15s Total obj = %7.3f\n\n", "", objValue);
	assert (fabs(objValue - minDetEqObj) < 1e-6 && "Check obj.");

	// report the portfolio at the root node
	double * rootSol = new double[nmbAssets];
	p2root->get_x_values(rootSol);
	printf ("Portfolio at the root node:\n");
	for (int a = 0; a < nmbAssets; ++a) {
		printf ("x(%d) = %6.2f\n", a, rootSol[a]);
	}
	printf ("\n");

	// cleaning
	delete[] rootSol;
	for (t = 0; t < (int) coreNodes.size(); t++) {
		delete coreNodes[t];
	}
	delete p2stochCore;

	return 0;
}
