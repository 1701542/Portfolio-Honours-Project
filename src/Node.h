#pragma once

#include "DXF.h"	// include dxframework
#include "cuda_runtime_api.h"

class Node
{
public:
	Node();
	~Node();

	float density;
	XMFLOAT3 position;
	XMINT3 indexPosition;
};

