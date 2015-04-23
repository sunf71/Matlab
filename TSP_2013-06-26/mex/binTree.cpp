// =============================================================================
// == binNode.cpp
// == --------------------------------------------------------------------------
// == A Binary Tree Node class
// == --------------------------------------------------------------------------
// == Written by Jason Chang 02-16-2008
// =============================================================================

#include "binTree.h"

// --------------------------------------------------------------------------
// -- binNode
// --   constructor; initializes the binary tree node to nothing
// --------------------------------------------------------------------------
template <typename T>
binNode<T>::binNode() : 
   data(), parent(NULL), child1(NULL), child2(NULL)
{
}

// --------------------------------------------------------------------------
// -- binNode
// --   constructor; initializes the binary tree node to contain the data
// --
// --   parameters:
// --     - new_data : the data to put in the new node
// --------------------------------------------------------------------------
template <typename T>
binNode<T>::binNode(T new_data) : 
   data(new_data), parent(NULL), child1(NULL), child2(NULL)

{
}

template <typename T>
binNode<T>::binNode(const binNode<T> &p)
{
   data = p.data;
   parent = p.parent;
   child1 = p.child1;
   child2 = p.child2;
}
   
template <typename T>
binNode<T>& binNode<T>::operator=(const binNode<T> &p)
{
   if (this != &p)
   {
      data = p.data;
      parent = p.parent;
      child1 = p.child1;
      child2 = p.child2;
   }
   return *this;
}

// --------------------------------------------------------------------------
// -- clearNode
// --   deletes dynamic data at node
// --------------------------------------------------------------------------
template <typename T>
void binNode<T>::clearNode()
{
   if (child1 != NULL) delete child1;
   if (child2 != NULL) delete child2;
}

// --------------------------------------------------------------------------
// -- getData
// --   retrieves the data at the node
// --
// --   return_value: the data of the node
// --------------------------------------------------------------------------
template <typename T>
T binNode<T>::getData()
{
   return data;
}



// =============================================================================
// == binTree.h
// == --------------------------------------------------------------------------
// == A Binary Tree class
// == --------------------------------------------------------------------------
// == Written by Jason Chang 12-14-2007
// =============================================================================

// --------------------------------------------------------------------------
// -- deleteSubTree
// --   a helper function used to recursively delete the tree
// --
// --   parameters:
// --     - root : the root of the subtree to delete
// --------------------------------------------------------------------------
template <typename T>
void binTree<T>::deleteSubTree(binNode<T>* rootSubTree)
{
   if (rootSubTree != NULL)
   {
      deleteSubTree(rootSubTree->child1);
      deleteSubTree(rootSubTree->child2);
      (*rootSubTree).clearNode();
   }
}

// --------------------------------------------------------------------------
// -- swapData
// --   a helper function that swaps the the parent and child
// --
// --   parameters:
// --     - orig_parent : the original parent
// --     - orig_child : the original child
// --------------------------------------------------------------------------
template <typename T>
void binTree<T>::swapData(binNode<T>* orig_parent, binNode<T>* orig_child)
{
   // check to see if we are changing root
   if (orig_parent == root)
      root = orig_child;

   // swap pointers
   binNode<T>* orig_sibling;
   if (orig_parent->child1 == orig_child)
      orig_sibling = orig_parent->child2;
   else
      orig_sibling = orig_parent->child1;
   if (orig_sibling != NULL)
      orig_sibling->parent = orig_child;

   binNode<T>* orig_grandparent = orig_parent->parent;
   if (orig_grandparent != NULL)
   {
      if (orig_grandparent->child1 == orig_parent)
         orig_grandparent->child1 = orig_child;
      else
         orig_grandparent->child2 = orig_child;
   }

   binNode<T>* orig_grandchild = orig_child->child1;
   if (orig_grandchild != NULL)
      orig_grandchild->parent = orig_parent;
   orig_parent->child1 = orig_grandchild;
   
   orig_grandchild = orig_child->child2;
   if (orig_grandchild != NULL)
      orig_grandchild->parent = orig_parent;
   orig_parent->child2 = orig_grandchild;

   orig_child->parent = orig_grandparent;
   orig_parent->parent = orig_child;
   orig_child->child1 = orig_parent;
   orig_child->child2 = orig_sibling;
}

// --------------------------------------------------------------------------
// -- binTree
// --   constructor; initializes the binary tree to nothing
// --------------------------------------------------------------------------
template <typename T>
binTree<T>::binTree() :
   root(NULL)
{
}

// --------------------------------------------------------------------------
// -- binTree
// --   destructor
// --------------------------------------------------------------------------
template <typename T>
binTree<T>::~binTree()
{
   if (root != NULL)
   {
      deleteSubTree(root);
      delete root;
   }
}

// --------------------------------------------------------------------------
// -- findLeaf
// --   finds a random leaf
// --
// --   parameters:
// --     - root : the roof of the subtree to start looking from
// --
// --   return_value: a pointer to a random binary tree leaf
// --------------------------------------------------------------------------
template <typename T>
binNode<T>* binTree<T>::findLeafSubTree(binNode<T>* rootSubTree)
{
   bool child1null = ((rootSubTree->child1) == NULL);
   bool child2null = ((rootSubTree->child2) == NULL);

   if (!child1null && !child2null)
   {
      // full node, traverse randomly down
      if (rand()%2)
         return findLeafSubTree(rootSubTree->child1);
      else
         return findLeafSubTree(rootSubTree->child2);
   }
   else if (child1null && !child2null)
   {
      // possible leaf node with child1 empty
      if (rand()%2)
         return rootSubTree;
      else
         return findLeafSubTree(rootSubTree->child2);
   }
   else if (!child1null && child2null)
   {
      // possible leaf node with child2 empty
      if (rand()%2)
         return findLeafSubTree(rootSubTree->child1);
      else
         return rootSubTree;
   }
   else
   {
      // leaf node
      return rootSubTree;
   }
}

// --------------------------------------------------------------------------
// -- addLeaf
// --   adds a leaf at a random location to the end
// --
// --   parameters:
// --     - new_data : the data of the new leave
// --
// --   return_value: a pointer to the newly added leaf
// --------------------------------------------------------------------------
template <typename T>
binNode<T>* binTree<T>::addLeaf(T new_data)
{
   // check if we have an empty tree
   if (root == NULL)
   {
      root = new binNode<T>(new_data);
      return root;
   }
   else
   {
      binNode<T>* new_parent = findLeafSubTree(root);
      binNode<T>* new_leaf = new binNode<T>(new_data);
      
      if (new_parent->child1 == NULL)
         new_parent->child1 = new_leaf;
      else
         new_parent->child2 = new_leaf;
      new_leaf->parent = new_parent;

      sortLeaf(new_leaf);
      return new_leaf;
   }
}

// --------------------------------------------------------------------------
// -- modifyLeaf
// --   modifies a leaf's information and then resorts the tree
// --
// --   parameters:
// --     - nodePtr : a pointer to the node you want to modify
// --     - new_data : the new_data to update the node with
// --------------------------------------------------------------------------
template <typename T>
void binTree<T>::modifyLeaf(binNode<T>* nodePtr, T new_data)
{
   nodePtr->data = new_data;
   sortLeaf(nodePtr);
}

// --------------------------------------------------------------------------
// -- sortLeaf
// --   sorts the tree from the new leaf given as a parameter
// --
// --   parameters:
// --     - new_leaf : the newly added leaf
// --------------------------------------------------------------------------
template <typename T>
void binTree<T>::sortLeaf(binNode<T>* new_leaf)
{
   // check if we are at root
   if (new_leaf->parent == NULL)
   {
      root = new_leaf;
      return;
   }

   else
   {
      binNode<T>* new_parent = new_leaf->parent;
      if (new_leaf->data < new_parent->data)
      {
         // swap the two nodes
         swapData(new_parent, new_leaf);

         // recursively call on the swapped parent, which is new_leaf
         sortLeaf(new_leaf);
      }
   }
}

// --------------------------------------------------------------------------
// -- pickOffRoot
// --   returns the root data, deletes it from the tree, and reorganizes so
// -- that it stays sorted
// --
// --   return_value: the data of the tree
// --------------------------------------------------------------------------
template <typename T>
T binTree<T>::pickOffRoot()
{
   if (root == NULL)
      return T();
   else
   {
      T root_data = root->data;
      root = sortRoot(root);
      return root_data;
   }
}

// --------------------------------------------------------------------------
// -- sortRoot
// --   sorts the tree from the root (ignores the data in the actual root)
// --
// --   parameters:
// --     - new_root : the root to look from
// --
// --   return_value : returns a pointer to the new root
// --------------------------------------------------------------------------
template <typename T>
binNode<T>*  binTree<T>::sortRoot(binNode<T>* new_root)
{
   // make sure that new_root is a valid pointer
   if (new_root == NULL)
      return NULL;
   else
   {
      binNode<T>* new_child1 = new_root->child1;
      binNode<T>* new_child2 = new_root->child2;

      bool child1null = (new_child1 == NULL);
      bool child2null = (new_child2 == NULL);

      // check to see if we have reached a leaf
      if (child1null && child2null)
      {
      /*
         std::cout<<"hello\n";
      
         if (new_root->parent != NULL)
         {
            std::cout<<new_root->parent->parent<<"\t";
            std::cout<<new_root->parent<<"\t";
            std::cout<<new_root->parent->child1<<"\t";
            std::cout<<new_root->parent->child2<<"\n";
         } */     
      
         // reached a leaf... delete it
         // update pointers to this node if it has a parent
         if (new_root->parent != NULL)
         {
            binNode<T>* new_parent = new_root->parent;
            if (new_root == new_parent->child1)
               new_parent->child1 = NULL;
            else
               new_parent->child2 = NULL;
         }
         /*
         std::cout<<root->parent<<"\t";
         std::cout<<root<<"\t";
         std::cout<<root->child1<<"\t";
         std::cout<<root->child2<<"\n";
         
         std::cout<<new_root->parent<<"\t";
         std::cout<<new_root<<"\t";
         std::cout<<new_root->child1<<"\t";
         std::cout<<new_root->child2<<"\n";
         
         if (new_root->parent != NULL)
         {
            std::cout<<new_root->parent->parent<<"\t";
            std::cout<<new_root->parent<<"\t";
            std::cout<<new_root->parent->child1<<"\t";
            std::cout<<new_root->parent->child2<<"\n";
         }
         
         if (root->child1 == new_root->parent || root->child2 == new_root->parent || root == new_root ||
             root == new_root->parent)
            printTree();
         */

//         std::cout<<"DATA: " <<new_root->data<<"\n";
         delete new_root;
//         std::cout<<"bye\n";
         return NULL;
      }
      else if (child2null || (!child1null && (new_child1->data <= new_child2->data)) )
      {
         // swap root with child1
         swapData(new_root, new_child1);
         sortRoot(new_root);
         return new_child1;
      }
      else
      {
         // swap root with child2
         swapData(new_root, new_child2);
         sortRoot(new_root);
         return new_child2;
      }
   }
}

// --------------------------------------------------------------------------
// -- treeNotEmpty
// --   returns true if the tree has data, false if it is empty
// --------------------------------------------------------------------------
template <typename T>
bool binTree<T>::treeNotEmpty()
{
   return (root!=NULL);
}

template <typename T>
void binTree<T>::checkErrorsSubTree(binNode<T>* rootSubTree)
{
   if (rootSubTree != NULL)
   {
      int x;
      if (rootSubTree->child1 != NULL && rootSubTree->child1->parent != rootSubTree)
      {
         std::cout<<"ERROR ON CHILD 1!!!\n";
         std::cin>>x;
         printSubTree(rootSubTree);
         std::cin>>x;
      }
      if (rootSubTree->child2 != NULL && rootSubTree->child2->parent != rootSubTree)
      {
         std::cout<<"ERROR ON CHILD 1!!!\n";
         std::cin>>x;
         printSubTree(rootSubTree);
         std::cin>>x;
      }
      checkErrorsSubTree(rootSubTree->child1);
      checkErrorsSubTree(rootSubTree->child2);
   }
}

// --------------------------------------------------------------------------
// -- printTree
// --   prints the tree.  used for debugging
// --------------------------------------------------------------------------
template <typename T>
void binTree<T>::printSubTree(binNode<T>* rootSubTree)
{
   if (rootSubTree != NULL)
   {
      // print the parent
      std::cout << rootSubTree->data << " ";
      printSubTree(rootSubTree->child1);
      std::cout << rootSubTree->data << " ";
      printSubTree(rootSubTree->child2);
      std::cout << rootSubTree->data << " ";
   }
}

// --------------------------------------------------------------------------
// -- printTree
// --   prints the tree.  used for debugging
// --------------------------------------------------------------------------
template <typename T>
void binTree<T>::printTree()
{
   printSubTree(root);
}

