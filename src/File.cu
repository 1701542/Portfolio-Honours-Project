#ifndef _CUDA_CU
#define _CUDA_CU

#include "Header.cuh"
#include "cuda_runtime_api.h"
#include "device_launch_parameters.h";
#include "DXF.h"
#include "Voxel.h"

#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>

#include <chrono>
using std::chrono::duration_cast;
using std::chrono::microseconds;
typedef std::chrono::steady_clock the_clock;

//
//	Interpolate between two nodes
__device__ void AddInterpolatedVertex(unsigned int voxelIter,
	unsigned int vertexIter,
	float* devVertexArray,
	float* devNodePositionArray,
	float* devNodeDensityArray,
	unsigned int node1,
	unsigned int node2,
	float devDensityCutoff)
{
	if (abs(devDensityCutoff - devNodeDensityArray[(voxelIter * 8) + node1]) < 0.00001f)
	{
		devVertexArray[(vertexIter * 3)] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3)];
		devVertexArray[(vertexIter * 3) + 1] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 1];
		devVertexArray[(vertexIter * 3) + 2] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 2];
		return;
	}
	if (abs(devDensityCutoff - devNodeDensityArray[(voxelIter * 8) + node2]) < 0.00001f)
	{
		devVertexArray[(vertexIter * 3)] = devNodePositionArray[(voxelIter * (8 * 3)) + (node2 * 3)];
		devVertexArray[(vertexIter * 3) + 1] = devNodePositionArray[(voxelIter * (8 * 3)) + (node2 * 3) + 1];
		devVertexArray[(vertexIter * 3) + 2] = devNodePositionArray[(voxelIter * (8 * 3)) + (node2 * 3) + 2];
		return;
	}
	if (abs(devNodeDensityArray[(voxelIter * 8) + node1] - devNodeDensityArray[(voxelIter * 8) + node2]) < 0.00001f)
	{
		devVertexArray[(vertexIter * 3)] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3)];
		devVertexArray[(vertexIter * 3) + 1] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 1];
		devVertexArray[(vertexIter * 3) + 2] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 2];
		return;
	}

	float mu = (devDensityCutoff - devNodeDensityArray[(voxelIter * 8) + node1]) / (devNodeDensityArray[(voxelIter * 8) + node2] - devNodeDensityArray[(voxelIter * 8) + node1]);
	
	devVertexArray[(vertexIter * 3)] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3)] + mu * (devNodePositionArray[(voxelIter * (8 * 3)) + (node2 * 3)] - devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3)]);
	devVertexArray[(vertexIter * 3) + 1] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 1] + mu * (devNodePositionArray[(voxelIter * (8 * 3)) + (node2 * 3) + 1] - devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 1]);
	devVertexArray[(vertexIter * 3) + 2] = devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 2] + mu * (devNodePositionArray[(voxelIter * (8 * 3)) + (node2 * 3) + 2] - devNodePositionArray[(voxelIter * (8 * 3)) + (node1 * 3) + 2]);
}

__global__ void March(float* devDensityCutoff,
	int* devVoxelArraySize,
	float* devNodeDensityArray,
	float* devNodePositionArray,
	int* devEdgeTable,
	int* devTriTable,
	float* devAllVertices,
	int* devVertCounter)
{
	int voxelIter = (blockIdx.x * blockDim.x) + threadIdx.x;

	if(voxelIter < *devVoxelArraySize)
	{
		unsigned char voxelByte = 0x00;

		if (devNodeDensityArray[(voxelIter * 8) + 1] < *devDensityCutoff) voxelByte |= 1;
		if (devNodeDensityArray[(voxelIter * 8) + 5] < *devDensityCutoff) voxelByte |= 2;
		if (devNodeDensityArray[(voxelIter * 8) + 4] < *devDensityCutoff) voxelByte |= 4;
		if (devNodeDensityArray[(voxelIter * 8) + 0] < *devDensityCutoff) voxelByte |= 8;
		if (devNodeDensityArray[(voxelIter * 8) + 3] < *devDensityCutoff) voxelByte |= 16;
		if (devNodeDensityArray[(voxelIter * 8) + 7] < *devDensityCutoff) voxelByte |= 32;
		if (devNodeDensityArray[(voxelIter * 8) + 6] < *devDensityCutoff) voxelByte |= 64;
		if (devNodeDensityArray[(voxelIter * 8) + 2] < *devDensityCutoff) voxelByte |= 128;

		if (devEdgeTable[voxelByte] != 0)
		{
			float vertexArray[12];

			if (devEdgeTable[voxelByte] & 1)	//	AND operator
				AddInterpolatedVertex(voxelIter, 0, vertexArray, devNodePositionArray, devNodeDensityArray, 1, 5, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 2)
				AddInterpolatedVertex(voxelIter, 1, vertexArray, devNodePositionArray, devNodeDensityArray, 5, 4, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 4)
				AddInterpolatedVertex(voxelIter, 2, vertexArray, devNodePositionArray, devNodeDensityArray, 4, 0, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 8)
				AddInterpolatedVertex(voxelIter, 3, vertexArray, devNodePositionArray, devNodeDensityArray, 0, 1, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 16)
				AddInterpolatedVertex(voxelIter, 4, vertexArray, devNodePositionArray, devNodeDensityArray, 3, 7, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 32)
				AddInterpolatedVertex(voxelIter, 5, vertexArray, devNodePositionArray, devNodeDensityArray, 7, 6, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 64)
				AddInterpolatedVertex(voxelIter, 6, vertexArray, devNodePositionArray, devNodeDensityArray, 6, 2, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 128)
				AddInterpolatedVertex(voxelIter, 7, vertexArray, devNodePositionArray, devNodeDensityArray, 2, 3, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 256)
				AddInterpolatedVertex(voxelIter, 8, vertexArray, devNodePositionArray, devNodeDensityArray, 1, 3, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 512)
				AddInterpolatedVertex(voxelIter, 9, vertexArray, devNodePositionArray, devNodeDensityArray, 5, 7, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 1024)
				AddInterpolatedVertex(voxelIter, 10, vertexArray, devNodePositionArray, devNodeDensityArray, 4, 6, *devDensityCutoff);
			if (devEdgeTable[voxelByte] & 2048)
				AddInterpolatedVertex(voxelIter, 11, vertexArray, devNodePositionArray, devNodeDensityArray, 0, 2, *devDensityCutoff);
		
			for (int vertIter = 0; devTriTable[(voxelByte * 16) + vertIter] != -1; vertIter += 3)
			{
				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 0) * 3) + 0] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter)] * 3) + 0];
				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 0) * 3) + 1] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter)] * 3) + 1];
				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 0) * 3) + 2] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter)] * 3) + 2];

				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 1) * 3) + 0] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter + 1)] * 3) + 0];
				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 1) * 3) + 1] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter + 1)] * 3) + 1];
				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 1) * 3) + 2] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter + 1)] * 3) + 2];

				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 2) * 3) + 0] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter + 2)] * 3) + 0];
				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 2) * 3) + 1] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter + 2)] * 3) + 1];
				devAllVertices[(voxelIter * (3 * 15)) + ((vertIter + 2) * 3) + 2] = vertexArray[(devTriTable[(voxelByte * 16) + (vertIter + 2)] * 3) + 2];

				atomicAdd(&devVertCounter[0], 3);
			}
		}
	}
}

struct is_non_zero
{
	__host__ __device__ bool operator()(const float x)
	{
		return x != 0;
	}
};

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
	double& freeTime)
{
	the_clock::time_point p1 = the_clock::now();

	//
	//	Create and Load intermediate arrays
	float* nodeDensityArray = new float[8 * voxelArraySize];
	float* nodePositionArray = new float[3 * 8 * voxelArraySize];
	
	for (int i = 0; i < voxelArraySize; ++i)
	{
		for (int j = 0; j < 8; ++j)
		{
			nodeDensityArray[(i * 8) + j] = voxelArray[i].getNode(j)->density;

			nodePositionArray[(i * (8 * 3)) + (j * 3) + 0] = voxelArray[i].getNode(j)->position.x;
			nodePositionArray[(i * (8 * 3)) + (j * 3) + 1] = voxelArray[i].getNode(j)->position.y;
			nodePositionArray[(i * (8 * 3)) + (j * 3) + 2] = voxelArray[i].getNode(j)->position.z;
		}
	}

	the_clock::time_point p2 = the_clock::now();
	nodeParseTime = duration_cast<microseconds>(p2 - p1).count();

	float* devDensityCutoff = 0;
	int* devVoxelArraySize = 0;
	float* devNodeDensityArray = 0;
	float* devNodePositionArray = 0;
	int* devEdgeTable = 0;
	int* devTriTable = 0;
	float* devAllVertices = 0;
	int* devVertCounter = 0;


	cudaError_t cudaStatus;

	//
	//	Malloc
	cudaStatus = cudaMallocManaged((void**)&devDensityCutoff, sizeof(float));
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMalloc failed!");

	cudaStatus = cudaMallocManaged((void**)&devVoxelArraySize, sizeof(int));
	if (cudaStatus != cudaStatus)
		fprintf(stderr, "cudaMalloc failed!");

	cudaStatus = cudaMallocManaged((void**)&devNodeDensityArray, 8 * voxelArraySize * sizeof(float));
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMalloc failed!");

	cudaStatus = cudaMallocManaged((void**)&devNodePositionArray, 3 * 8 * voxelArraySize * sizeof(float));
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMalloc failed!");

	cudaStatus = cudaMallocManaged((void**)&devEdgeTable, 256 * sizeof(int));
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMalloc failed!");

	cudaStatus = cudaMallocManaged((void**)&devTriTable, 256 * 16 * sizeof(int));
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMalloc failed!");

	cudaStatus = cudaMallocManaged((void**)&devVertCounter, sizeof(int));
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMalloc failed!");

	cudaStatus = cudaMallocManaged((void**)&devAllVertices, 15 * 3 * voxelArraySize * sizeof(float));
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMalloc failed!");

	the_clock::time_point p3 = the_clock::now();
	mallocTime = duration_cast<microseconds>(p3 - p2).count();

	//
	//	Initialise arrays with values
	cudaStatus = cudaMemcpy(devDensityCutoff, &densityCutoff, sizeof(float), cudaMemcpyHostToDevice);
	if(cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMemcpy failed!");

	cudaStatus = cudaMemcpy(devVoxelArraySize, &voxelArraySize, sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMemcpy failed!");

	cudaStatus = cudaMemcpy(devNodeDensityArray, nodeDensityArray, 8 * voxelArraySize * sizeof(float), cudaMemcpyHostToDevice);
	if(cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMemcpy failed!");

	cudaStatus = cudaMemcpy(devNodePositionArray, nodePositionArray, 3 * 8 * voxelArraySize * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMemcpy failed!");

	cudaStatus = cudaMemcpy(devEdgeTable, edgeTable, 256 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMemcpy failed!");

	cudaStatus = cudaMemcpy(devTriTable, triTable, 256 * 16 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaMemcpy failed!");

	//
	//	Delete intermediate dynamic arrays
	delete nodeDensityArray;
	nodeDensityArray = 0;

	delete nodePositionArray;
	nodePositionArray = 0;

	the_clock::time_point p4 = the_clock::now();
	memcpyTime = duration_cast<microseconds>(p4 - p3).count();

	//
	//	Optimise thread hierarchies
	int numThreads = voxelArraySize % 32 == 0 ? voxelArraySize : ((voxelArraySize / 32) + 1.0f) * 32;
	int numBlocks = 1;

	if (numThreads > 1024)
	{
		numBlocks = numThreads % 1024 == 0 ? (numThreads / 1024) : (numThreads / 1024) + 1;
		numThreads = numThreads / numBlocks;
	}

	dim3 blocks(numBlocks);
	dim3 threads(numThreads);

	//
	//	Run
	March << <blocks, threads>> > (devDensityCutoff,
		devVoxelArraySize,
		devNodeDensityArray,
		devNodePositionArray,
		devEdgeTable,
		devTriTable,
		devAllVertices,
		devVertCounter);

	//
	//	Check error
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));

	//
	//	Sync
	cudaDeviceSynchronize();

	the_clock::time_point p5 = the_clock::now();
	marchTime = duration_cast<microseconds>(p5 - p4).count();

	//
	//	Compact verts and indices
	if (cudaStatus == cudaSuccess)
	{

		thrust::device_vector<float> t_devAllVertices(devAllVertices, devAllVertices + (voxelArraySize * 3 * 15));

		thrust::device_vector<float> t_compactAllVertices(voxelArraySize * 3 * 15, 0);

		thrust::copy_if(thrust::device, t_devAllVertices.begin(), t_devAllVertices.end(), t_compactAllVertices.begin(),is_non_zero());

		thrust::copy(t_compactAllVertices.begin(), t_compactAllVertices.end(), devAllVertices);

		cudaStatus = cudaMemcpy(allVertices, devAllVertices, 15 * 3 * voxelArraySize * sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
			fprintf(stderr, "cudaMemcpy failed!");

		cudaStatus = cudaMemcpy(&vertCounter, devVertCounter,sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
			fprintf(stderr, "cudaMemcpy failed!");
		
		thrust::sequence(thrust::host, allIndices, allIndices + vertCounter, 0);
	}

	the_clock::time_point p6 = the_clock::now();
	compactTime = duration_cast<microseconds>(p6 - p5).count();

	//
	//	Free
	cudaFree(devDensityCutoff);
	cudaFree(devVoxelArraySize);
	cudaFree(devNodeDensityArray);
	cudaFree(devNodePositionArray);
	cudaFree(devEdgeTable);
	cudaFree(devTriTable);
	cudaFree(devAllVertices);
	cudaFree(devVertCounter);

	the_clock::time_point p7 = the_clock::now();
	freeTime = duration_cast<microseconds>(p7 - p6).count();
}

#endif