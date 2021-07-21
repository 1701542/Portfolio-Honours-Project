#pragma once

#include <cstdint>
#include "Node.h"

class Voxel
{
public:
	Voxel();
	~Voxel();

	void setNode(unsigned int nodeNum, Node& node);
	Node* getNode(unsigned int nodeNum) { return nodes[nodeNum]; };

	void setPosition(XMFLOAT3 newPosition);
	XMFLOAT3 getPosition() { return position; };

private:

	//	A pointer is used instead of a reference because it can be re-assigned
	/*
	*	|/	  |/
	*	3-----7---
	* |/|	|/|
	* 2-----6---
	* |	1---|-5---
	* |/	|/
	* 0-----4---
	* 
	*/
	Node* nodes[8];
	XMFLOAT3 position;
};

