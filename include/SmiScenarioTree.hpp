// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef SmiScenarioTree_H
#define SmiScenarioTree_H

/** Scenario Tree

  This class is used for storing and accessing scenario trees.

*/
#include <list>
#include <vector>
#include <map>
#include <iostream>
#include <assert.h>


using namespace std;

/** SmiTreeNode template class.

    Manages pointers to parent, child and sibling for tree navigation.
	Template class instance is a pointer to an object that
	must be created with "new" operator.
*/


template <class T> 
class SmiTreeNode  {
   friend void SmiTreeNodeUnitTest();
public:

	typedef map<int,SmiTreeNode<T>*> child_label_map;

	bool hasParent() { return (parent_!=NULL); }
	bool hasChild()  { return (child_!=NULL);  }
	bool hasSibling(){ return (sibling_!=NULL);}

	SmiTreeNode<T>  *getParent() { return parent_; }
	SmiTreeNode<T>  *getChild()  { return child_;  }
	SmiTreeNode<T>  *getSibling(){ return sibling_;}

	SmiTreeNode<T>  *getChildByLabel(int n){
		child_label_map::iterator begpos = child_labels_.begin();
	    child_label_map::iterator endpos = child_labels_.end();
		while(begpos!=endpos)
		{
			printf(" found label %d \n",begpos->first);
			++begpos;
		}
		child_label_map::iterator pos = child_labels_.find(n);
		if (pos!=child_labels_.end())
			return pos->second;
		else
			return NULL;
	}


	int depth() { return depth_; }
	int numChildren() { return nchild_; }
	

	SmiTreeNode<T> * addChild(T cd, int label=-1)
	{
		SmiTreeNode<T> *c = new SmiTreeNode(cd);
		if (label==-1) label=nchild_;
		child_labels_.insert(make_pair(label,c));
		//debug code
		child_label_map::iterator pos = child_labels_.find(label);
		assert (pos!=child_labels_.end());
		//
		c->parent_     = this;
		c->depth_      = depth_ + 1;
		c->sibling_    = child_;
		nchild_++;
		child_ = c;
		return c;
	}

	vector<SmiTreeNode<T> *> &getChildren()
	{
		SmiTreeNode<T> *pnode = this;
		int i = this->numChildren();
		vector<SmiTreeNode<T> *> *vec = new vector<SmiTreeNode<T> *>(i);
		(*vec)[--i] = pnode = pnode->getChild();
		while (i>0)
			(*vec)[--i] = pnode = pnode->getSibling();
		return *vec;
	}

	T getDataPtr() { return ptr_; }

//--------------------------------------------------------------------------
   /**@name Constructors, destructors and major modifying methods*/
   //@{
   /// Default Constructor creates an empty node
   SmiTreeNode<T>(){
	parent_=NULL;
	child_=NULL;
	sibling_= NULL;
	nchild_ = 0;
	depth_ = 0;
   }

   /// Constructor from P
   SmiTreeNode<T>(T p){
	parent_=NULL;
	child_=NULL;
	sibling_= NULL;
	nchild_ = 0;
	depth_ = 0;
	ptr_ = p;
   }

   /// Destructor 
   ~SmiTreeNode()
   {
	   delete sibling_;
	   delete child_;
	   delete ptr_ ;
   }
   //@}

protected:
	
	void setChild  (SmiTreeNode<T>  *c) {child_ = c;}
	void setSibling(SmiTreeNode<T>  *s) {sibling_ = s;}
	SmiTreeNode<T>  *getParentP() { return parent_; }
	SmiTreeNode<T>  *getChildP()  { return child_;  }
	SmiTreeNode<T>  *getSiblingP(){ return sibling_;}
	

private:
	SmiTreeNode<T>  *parent_;
	SmiTreeNode<T>  *child_;
	SmiTreeNode<T>  *sibling_;
	int nchild_;
	int depth_;
	T ptr_;
	child_label_map child_labels_;

};

//#############################################################################
/** A function that tests the methods in the SmiTreeNode class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
SmiTreeNodeUnitTest();


template<class T>
class SmiScenarioTree  {
   friend void SmiScenarioTreeUnitTest();
public:

//---------------------------------------------------------------------------
  /**@name Iterators */
  //@{


	/** begin */
	typename vector<T>::iterator treeBegin(){ return node_data.begin(); }
	/** end */
	typename vector<T>::iterator treeEnd(){ return node_data.end();}
	/** whole tree */
	vector<T> &wholeTree(){ return node_data;}

	/** scenario iterators
	TODO: native code for these iterators that does not
	depend on copying. */
	typename vector<T>::iterator scenBegin(int s){
		getScenario(s);
		return scen_data.begin();
	}

	typename vector<T>::iterator scenEnd(int s){
		getScenario(s);
		return scen_data.begin()+leaf_[s]->depth()+1;
	}




//---------------------------------------------------------------------------
  /**@name Query members */
  //@{

	/** Get root node.	*/
	SmiTreeNode<T> *getRoot() { return root_; }

	/** Get leaf node. */
	SmiTreeNode<T> *getLeaf(int scn) { return leaf_[scn];}

	/** Get node identified by scenario/stage.	*/
	SmiTreeNode<T> &find(unsigned int scenario, int stage)
	{
		assert (scenario < leaf_.size());		
		SmiTreeNode<T> * n = leaf_[scenario];
		assert (stage < n->depth() + 1);
		while (stage < n->depth())
			n = n->getParent();
		return *n;
	}

	/** Get node identified by longest match to array of labels */
	SmiTreeNode<T> &find(vector<int> &label)
	{
		assert(label.size()>0);
		SmiTreeNode<T> *n = root_,*next;
		int i=1;
		while ((i<label.size()) && (next=n->getChildByLabel(label[i])))
		{
			++i;
			n=next;
		}
		return *n;
	}


	/** Get vector of node data for given scenario */
	vector<T> &getScenario(int scenario)
	{
		assert (scenario < leaf_.size());
		SmiTreeNode<T> * n = leaf_[scenario];

//		if ( n->getDataPtr()==scen_data[n->depth()] ) return scen_data;

		int ns=n->depth()+1-scen_data.size();
		for (int j=0; j<ns;j++)
			scen_data.push_back(n->getDataPtr());

		int i=n->depth()+1;
		while(i>0)
		{
			scen_data[--i] = n->getDataPtr();
			n = n->getParent();
		}
		return scen_data;
	}


  //@}

//---------------------------------------------------------------------------
  /**@name Tree modification members */
  //@{
    /** Add path from node id'd by scenario and stage.
	    Responsibility for memory management of SmiTreeNodeData elements
		is assigned to SmiScenarioTree.
		SmiTreeNodeData elements must be created with "new" operator.
	*/
	int addPathtoLeaf(int brscenario, int stage, vector<T> &pathdata)
	{
		SmiTreeNode<T> *parent = NULL;		
		if (leaf_.size())
			parent = &find(brscenario,stage);
		for ( unsigned int i = 0; i < pathdata.size(); i++)
		{	
			if (parent)
			{
				// no label
				parent = parent->addChild(pathdata[i]);
			}
			else
			{
				parent = root_ = new SmiTreeNode<T>(pathdata[0]);
			}
			// add data to full node_data array
			node_data.push_back(pathdata[i]);
		}
		if (pathdata.size())
		{
			leaf_.push_back(parent);			
		}
		return leaf_.size()-1;
		
	}
	/** Add path from node id'd by matching node labels.
	    Responsibility for memory management of SmiTreeNodeData elements
		is assigned to SmiScenarioTree.
		SmiTreeNodeData elements must be created with "new" operator.
	*/
	int addPathtoLeaf(vector<int> &label, vector<T> &pathdata)
	{
		SmiTreeNode<T> *parent = NULL;		
		if (leaf_.size())
			parent = &find(label);
		int stage=0;
		if (parent)
			stage=parent->depth()+1;

		for (int i=0 ; i < pathdata.size(); i++)
		{	
			if (parent)
			{
				//has label
				//printf(" put label %d at depth %d\n ",label[i+stage],parent->depth());
				parent = parent->addChild(pathdata[i],label[i+stage]);
			}
			else
			{
				parent = root_ = new SmiTreeNode<T>(pathdata[0]);
			}
			// add data to full node_data array
			node_data.push_back(pathdata[i]);
		}
		if (pathdata.size())
		{
			leaf_.push_back(parent);			
		}
		return leaf_.size()-1;
	}

  //@}

//--------------------------------------------------------------------------
   /**@name Constructors, destructors and major modifying methods*/
   //@{
   /// Default Constructor creates an empty scenario tree
	SmiScenarioTree<T>(): leaf_(0),root_(NULL) {};

   /// Destructor 
   virtual ~SmiScenarioTree<T>() {delete root_;}
   //@}

private:
	vector<T> node_data;
	vector<T> scen_data;
	vector<SmiTreeNode<T> *> leaf_;
	SmiTreeNode<T> *root_;
};

//#############################################################################
/** A function that tests the methods in the SmiScenarioTree class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
SmiScenarioTreeUnitTest();

#endif //SmiScenarioTree_H
