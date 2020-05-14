/*
 * TreeNode.hpp
 *
 *  Created on: Jan 24, 2016
 *      Author: michel
 */
#ifndef TREENODE_HPP
#define TREENODE_HPP

#include <float.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>

class TreeNode{
	public:
		TreeNode();
		TreeNode(std::string *name);
		~TreeNode();

		TreeNode *leftChild = NULL;
		TreeNode *rightChild = NULL;
		std::string *name;
		float length = FLT_MAX;

		void buildTreeString (std::string *sb);

};

#endif
