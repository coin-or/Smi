// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef SmiTreeNode_H
#define SmiTreeNode_H

#include <iostream>

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

	bool hasParent() { return (parent_!=NULL); }
	bool hasChild()  { return (child_!=NULL);  }
	bool hasSibling(){ return (sibling_!=NULL);}

	SmiTreeNode<T>  *getParent() { return parent_; }
	SmiTreeNode<T>  *getChild()  { return child_;  }
	SmiTreeNode<T>  *getSibling(){ return sibling_;}

	int depth() { return depth_; }
	int numChildren() { return nchild_; }
	
	SmiTreeNode<T> * addChild(T cd)
	{
		SmiTreeNode<T> *c = new SmiTreeNode(cd);
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

};

//#############################################################################
/** A function that tests the methods in the SmiTreeNode class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
SmiTreeNodeUnitTest();

#endif //SmiTreeNode_H
