// =============================================================================
// == binNode.h
// == --------------------------------------------------------------------------
// == A Binary Tree Node class
// == --------------------------------------------------------------------------
// == Written by Jason Chang 02-16-2008
// =============================================================================

#ifndef BINNODE_H
#define BINNODE_H

#include "assert.h"
#include "string.h"

template <typename T>
class binNode
{
private:
   T data;
   binNode* parent;
   binNode* child1;
   binNode* child2;
   template <typename TT> friend class binTree;

public:
   // --------------------------------------------------------------------------
   // -- binNode
   // --   constructor; initializes the binary tree node to nothing
   // --------------------------------------------------------------------------
   binNode();

   // --------------------------------------------------------------------------
   // -- binNode
   // --   constructor; initializes the binary tree node to contain the data
   // --
   // --   parameters:
   // --     - new_data : the data to put in the new node
   // --------------------------------------------------------------------------
   binNode(T new_data);
   binNode(const binNode &p);
   binNode& operator=(const binNode &p);

   // --------------------------------------------------------------------------
   // -- clearNode
   // --   deletes dynamic data at node
   // --------------------------------------------------------------------------
   void clearNode();

   // --------------------------------------------------------------------------
   // -- getData
   // --   retrieves the data at the node
   // --
   // --   return_value: the data of the node
   // --------------------------------------------------------------------------
   T getData();
};
#endif /* BINNODE_H */


// =============================================================================
// == binTree.h
// == --------------------------------------------------------------------------
// == A Binary Tree class
// == --------------------------------------------------------------------------
// == Written by Jason Chang 12-14-2007
// =============================================================================

#ifndef BINTREE_H
#define BINTREE_H


#include <cstdlib>
#include <stdlib.h>
#include <iostream>

template <typename T>
class binTree
{
private:
   binNode<T>* root;

   // --------------------------------------------------------------------------
   // -- deleteSubTree
   // --   a helper function used to recursively delete the tree
   // --
   // --   parameters:
   // --     - rootSubTree : the root of the subtree to delete
   // --------------------------------------------------------------------------
   void deleteSubTree(binNode<T>* rootSubTree);

   // --------------------------------------------------------------------------
   // -- findLeaf
   // --   finds a random leaf
   // --
   // --   parameters:
   // --     - rootSubTree : the roof of the subtree to start looking from
   // --
   // --   return_value: a pointer to a random binary tree leaf
   // --------------------------------------------------------------------------
   binNode<T>* findLeafSubTree(binNode<T>* rootSubTree);

   // --------------------------------------------------------------------------
   // -- swapData
   // --   a helper function that swaps the the parent and child
   // --
   // --   parameters:
   // --     - orig_parent : the original parent
   // --     - orig_child : the original child
   // --------------------------------------------------------------------------
   void swapData(binNode<T>* orig_parent, binNode<T>* orig_child);

   // --------------------------------------------------------------------------
   // -- sortLeaf
   // --   sorts the tree from the new leaf given as a parameter
   // --
   // --   parameters:
   // --     - new_leaf : the newly added leaf
   // --
   // --   return_value : returns a pointer to the new root
   // --------------------------------------------------------------------------
   void sortLeaf(binNode<T>* new_leaf);

   // --------------------------------------------------------------------------
   // -- sortRoot
   // --   sorts the tree from the root (ignores the data in the actual root)
   // --
   // --   parameters:
   // --     - new_root : the root to look from
   // --------------------------------------------------------------------------
   binNode<T>* sortRoot(binNode<T>* new_root);

   // --------------------------------------------------------------------------
   // -- printSubTree
   // --   prints the tree.  used for debugging
   // --------------------------------------------------------------------------
   void printSubTree(binNode<T>* rootSubTree);

public:
   // --------------------------------------------------------------------------
   // -- binTree
   // --   constructor; initializes the binary tree to nothing
   // --------------------------------------------------------------------------
   binTree();
   
   // --------------------------------------------------------------------------
   // -- binTree
   // --   destructor
   // --------------------------------------------------------------------------
   ~binTree();

   // --------------------------------------------------------------------------
   // -- addLeaf
   // --   adds a leaf at a random location to the end
   // --
   // --   parameters:
   // --     - new_data : the data of the new leave
   // --
   // --   return_value: a pointer to the newly added leaf
   // --------------------------------------------------------------------------
   binNode<T>* addLeaf(T new_data);

   // --------------------------------------------------------------------------
   // -- modifyLeaf
   // --   modifies a leaf's information and then resorts the tree
   // --
   // --   parameters:
   // --     - nodePtr : a pointer to the node you want to modify
   // --     - new_data : the new_data to update the node with
   // --------------------------------------------------------------------------
   void modifyLeaf(binNode<T>* nodePtr, T new_data);

   // --------------------------------------------------------------------------
   // -- pickOffRoot
   // --   returns the root data, deletes it from the tree, and reorganizes so
   // -- that it stays sorted
   // --
   // --   return_value: the data of the tree
   // --------------------------------------------------------------------------
   T pickOffRoot();
   
   void checkErrorsSubTree(binNode<T>* rootSubTree);

   // --------------------------------------------------------------------------
   // -- treeNotEmpty
   // --   returns true if the tree has data, false if it is empty
   // --------------------------------------------------------------------------
   bool treeNotEmpty();

   // --------------------------------------------------------------------------
   // -- printTree
   // --   prints the tree.  used for debugging
   // --------------------------------------------------------------------------
   void printTree();
};


#endif /* BINTREE_H */

