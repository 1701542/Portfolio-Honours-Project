#ifndef _CUDA_KERNEL_CUH
#define _CUDA_KERNEL_CUH

#include "Voxel.h"

void RunMarchingCubes(float densityCutoff,
	Voxel* voxelArray,
	int voxelArraySize,
	XMFLOAT3* allVertices,
	int* allIndices,
	int& vertCounter,
	int* triTable,
	int* edgeTable,
	double& nodeParseTime,
	double& mallocTime,
	double& memcpyTime,
	double& marchTime,
	double& compactTime,
	double& freeTime);

#endif