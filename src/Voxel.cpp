#include "Voxel.h"

Voxel::Voxel()
{

}

Voxel::~Voxel()
{

}

void Voxel::setNode(unsigned int nodeNum, Node& node)
{
	assert(nodeNum < 8 && nodeNum >= 0);
	nodes[nodeNum] = &node;
}

void Voxel::setPosition(XMFLOAT3 newPosition)
{
	position = newPosition;
}