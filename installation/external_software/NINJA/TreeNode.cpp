/*
 * TreeNode.cpp
 *
 *  Created on: Jan 24, 2016
 *      Author: michel
 */

#include "TreeNode.hpp"

TreeNode::TreeNode(){
	this->name = new std::string();
	this->name->assign("");
	this->leftChild = NULL;
	this->rightChild = NULL;
}

TreeNode::TreeNode(std::string *name){
	this->name = name;
}
TreeNode::~TreeNode(){
	if (this->rightChild != NULL)
		delete(this->rightChild);
	if (this->leftChild != NULL)
		delete(this->leftChild);
	if( this->name != NULL)
		delete this->name;
}
void TreeNode::buildTreeString (std::string *sb){
	std::string len;
	//printf("name: %s ", this->name);
	//printf("%.5f\n",this->length);
	if (this->length == FLT_MAX){
		len.assign("");
	}else{
		char x[20];
		sprintf(x,"%.5f",this->length);
		len.assign(":");
		len.append(x);
	}
	if (this->leftChild == NULL) { //|| this->rightChild == NULL
			sb->append(*this->name  + len);
	} else {
		sb->append("(");
		this->leftChild->buildTreeString(sb);
		sb->append(",");
		this->rightChild->buildTreeString(sb);
		sb->append(")" + len);
	}
}



