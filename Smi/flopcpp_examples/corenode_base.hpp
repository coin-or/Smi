/**
	\brief Core-node base class for the FlopCpp-Smi modelling framework
	\author Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file includes definitions of the general, problem independent,
	core-node base class.

	Note that the code is meant as an illustrative example that mixes different
	styles to show more ways of doing things, something you most likely do
	\e not want to do in a real code.
	In addition, in a real code one would probably made many of the members
	private and write get/set methods where needed.
**/

#ifndef CORENODE_BASE_HPP
#define CORENODE_BASE_HPP


/** Note that the COIN-OR classes should be included as "coin/class-name",
    but some of the FlopC++ classes do not do that, so we have to call
    the class without the prefix and add the coin dir. to the search path. **/
#include "flopc.hpp"

namespace FlopSmiEx {
using namespace flopc;
using namespace std;



/// Problem-independent base class for an LP/MIP model in one node of a tree.
class CoreNodeBase {
protected:
	/// static pointer to the current CoreNodeBase object
	/** Each time a variable or constraint is constructed, it should register
	    itself with the node it belongs to. However, this means that they need
	    a pointer to the node object. The easiest option is to pass a pointer
	    to the constructor, but this would clutter the code. Instead, we declare
	    this \em static variable and let the \c CoreNodeBase's constructor set
	    it to itself. Since every node calls the constructor of the base class
	    first, \c p2activeNode will point to the correct node at the time the
	    node's variables and constraints get constructed. **/
	static CoreNodeBase * p2activeNode;

	/// vector of scenario solution values, for use in \c SP_variable
	/** The \c value() member of \c SP_variable needs to have access to scenario
	    solution values. One option would be to pass a pointer, but this would
	    make the function a bit more difficult to read. Instead, we point this
	    \em static pointer to the solution and then use this one in \c value()
	    instead. **/
	static double const * p2varValues;

	MP_expression objFuncNode; ///< objective function in this node
	MP_expression objFuncRec;  ///< objective function in this node and below

public:
	/// constructor
	CoreNodeBase(CoreNodeBase *p2pred, int nodeStage = -1)
	:  p2parent(p2pred), p2child(NULL), stage(nodeStage)
	{
		if (p2parent != NULL) {
			assert (p2parent->p2child == NULL && "only one child per node in Core");
			p2parent->p2child = this; // register with the parent
		}
		CoreNodeBase::p2activeNode = this; // point the static pointer to itself
	}

	virtual ~CoreNodeBase() {} ///< destructor

	CoreNodeBase * p2parent; ///< pointer to the node's parent
	CoreNodeBase * p2child;  ///< pointer to the node's child
	int stage;               ///< stage of the node (used only for reporting)

	/// \name references to variables and constraints
	/** The idea is to have "meta objects" with all variables and constraints
	    in a node. This is important in creation of the Smi object, where we
	    have to associate the nodes variables and constraints to stages.
	    Without these new objects, we would have to access all the derived
	    classes independently, cluttering the code. **/
	///@{
	vector<VariableRef const *> all_variables; ///< list of ref. to variables
	vector<MP_constraint *> all_constraints;   ///< list of ref. to constraints
	vector<int> constr_row_offsets;            ///< list of row offsets
	///@}

	/// point a given pointer to the vector of scenario solution values
	void set_var_values(double const *p2values) {
		p2varValues = p2values;
	}

	/// make and return the recursive objective function
	/** This probably makes sense only in root nodes. **/
	MP_expression & get_obj_func() {
		if (p2parent != NULL) {
			cerr << "Warning: get_obj_func called from a non-root node!" << endl;
		}
		make_obj_func_rec();
		return objFuncRec;
	}

protected:
	/// class for a stochastic (node-based) variable
	/** Basically, this is identical to the MP_variable, but it registers itself
	    with the scenario-tree node and adds a new \c value() method. <p>
	    In this example, we have implemented only variables indexed on simple
	    one- or two-dimensional sets. **/
	class SP_variable : public MP_variable {
	private:
		CoreNodeBase * p2node; ///< pointer to the node the var belongs to
		int varIndx; ///< index of the 1st el. in the vector of node's variables
	public:
		SP_variable(MP_set const &s1 = MP_set::getEmpty(),
		            MP_set const &s2 = MP_set::getEmpty())
		: MP_variable(s1, s2), p2node(CoreNodeBase::p2activeNode),
			varIndx(static_cast<int>(p2node->all_variables.size()))
		{
			#ifndef NDEBUG
			cout << "Stage " << p2node->stage << ": adding SP_variable with "
			     << s1.size() * s2.size() << " element(s), starting at col "
			     << varIndx << endl;
			#endif
			for (int i1 = 0; i1 < s1.size(); i1++) {
				for (int i2 = 0; i2 < s2.size(); i2++) {
					p2node->all_variables.push_back(& this->operator()(i1,i2));
				}
			}
		}

		/// get the value of a given variable
		/**
			For this to work, the CoreNodeBase::p2varValues must be pointed to the
			vector of scenario solution, using \c CoreNodeBase::set_var_values().

			Note that we use the \c f method inherited from \c MP_variable
			to compute the offset for a given combination of indices.
		**/
		double value(int const i1 = 0, int const i2 = 0) const {
			assert (CoreNodeBase::p2varValues != NULL
			        && "variable values must be set before they are accessed!");
			return CoreNodeBase::p2varValues[varIndx + f(i1, i2)];
		}
	};

	/// class for a stochastic (node-based) constraint
	/**
		Basically, this is identical to the MP_constraint, but it registers itself
		with the scenario-tree node and adds a new \c value() method. <p>
		So far, we implement only variables indexed on simple one-dimensional sets.

		Note that constraints are different from variables, since we cannot get
		a pointer to one constraint in the set (one row). This is because if
		we have an \c MP_constraint \c C, then
		- \c C(i) sets the internal indices to point to the right row
		  and returns a pointer to \c C (as pointer to \c C(i) does not exist);
		  by default, the internal indices point to the first row.
		- \c C() returns the row of the currently selected row, based on the
		  values of the internal indices; this is why \c C(i)() returns the
		  row corresponding to the i-th constraint!
		- \c C.row_number(i) is the same as \c(), but with bound checking.
		This implies that even if we put pointers/references to \c C(i) into
		\c all_constraints, we get repeated pointers to the same object. One
		way around this is to detect this when we process the list, but that
		would mean an extra burden on the user. Instead, we add a vector
		\c constr_row_offsets of offsets that are added to the rows.
	*/
	class SP_constraint : public MP_constraint {
	private:
		CoreNodeBase * p2node; ///< pointer to the node the constr. belongs to
		int constrIndx;
	public:
		SP_constraint (MP_set const &s1 = MP_set::getEmpty())
		: MP_constraint(s1), p2node(CoreNodeBase::p2activeNode),
		  constrIndx(static_cast<int>(p2node->all_constraints.size()))
		{
			#ifndef NDEBUG
			cout << "Stage " << p2node->stage << ": adding SP_constr.  with "
			     << s1.size() << " element(s), starting at row " << constrIndx
			     << endl;
			#endif
			for (int i1 = 0; i1 < s1.size(); i1++) {
				p2node->all_constraints.push_back(this);
				p2node->constr_row_offsets.push_back(i1);
			}
		}
	};

	/// create the objective function expression, recursively for all children
	/**
		This function is protected, as it only makes sense to call it in
		the root, to create the complete objective function.
	**/
	void make_obj_func_rec();
};


} // namespace FlopSmiEx
#endif
