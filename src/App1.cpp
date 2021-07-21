// Lab1.cpp
// Lab 1 example, simple coloured triangle mesh
#include "Header.cuh"
#include "App1.h"
#include "OpenSimplexNoise.hpp"
#include <fstream>
#include <iostream>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::microseconds;
typedef std::chrono::steady_clock the_clock;

App1::App1()
{
	m_BasicShader = nullptr;
}

void App1::init(HINSTANCE hinstance, HWND hwnd, int screenWidth, int screenHeight, Input *in, bool VSYNC, bool FULL_SCREEN)
{
	BaseApplication::init(hinstance, hwnd, screenWidth, screenHeight, in, VSYNC, FULL_SCREEN);

	textureMgr->loadTexture( L"brick", L"res/Brick.png" );

	m_Isosurface = new Isosurface(renderer->getDevice(), renderer->getDeviceContext());
	m_BasicShader = new BasicShader(renderer->getDevice(), hwnd);

	ambientColour = XMFLOAT4(0.25f, 0.25f, 0.25f, 1.0f);
	diffuseColour = XMFLOAT4(0.75f, 0.75f, 0.75f, 1.0f);
	direction = XMFLOAT3(0.3f, 0.7f, 0.0f);

	light = new Light;
	light->setAmbientColour(ambientColour.x, ambientColour.y, ambientColour.z, ambientColour.w);
	light->setDiffuseColour(diffuseColour.x, diffuseColour.y, diffuseColour.z, diffuseColour.w);
	light->setDirection(direction.x, direction.y, direction.z);

	camera->setPosition( 0.0f, 2.0f, -15.0f );
	camera->setRotation( 0.0f, 0.0f, 0.0f );

	//	Init random values
	randomRange = 100;
	minValue = 10;

	terrainSize = 10;
	voxelStackSize = 10;
	densityCutoff = 50;
	useGPU = true;

	time_taken = 0.0f;
	gpuCompactTime = 0.0f;
	gpuMallocTime = 0.0f;
	gpuMarchTime = 0.0f;

	//	Init ImGui testing values
	testIterations = 5;
	testStepSize = 1;
	testUseGPU = false;
	testDensity = 50.0f;
	testStackSize = 10;
	testSurface = 50.0f;
	testScale = 0.1f;

	testDensityRange = XMFLOAT2(0.0f, 100.0f);
	testStackSizeRange = XMINT2(1,100);
	testSurfaceRange = XMFLOAT2(0,100);
	testScaleRange = XMFLOAT2(0.1f,0.5f);


	InitialiseStack();

	NoisifyDensities();

	MarchingCubes();
}


App1::~App1()
{
	// Run base application deconstructor
	BaseApplication::~BaseApplication();

	if (m_BasicShader)
	{
		delete m_BasicShader;
		m_BasicShader = 0;
	}
	if (m_Isosurface) {
		delete m_Isosurface;
		m_Isosurface = 0;
	}
}

//Render a bunch of instanced cubes
void App1::BuildCubeInstances()
{
	//	Render voxels
	XMFLOAT3* pos = new XMFLOAT3[voxelStackSize * voxelStackSize * voxelStackSize];

	int instanceCount = 0;

	for (int i = 0; i < voxelVector.size(); ++i)
	{
		pos[instanceCount] = voxelVector[i].getPosition();
		instanceCount++;
	}

	delete[] pos;
	pos = 0;
}

void App1::InitialiseStack()
{
	XMFLOAT3 position = XMFLOAT3(0.0f, 0.0f, 0.0f);
	
	if (!nodeVector.empty())
		nodeVector.clear();
	if (!voxelVector.empty())
		voxelVector.clear();

	//	Reserve enough contiguous memory for vectors
	nodeVector.reserve((voxelStackSize + 1) * (voxelStackSize + 1) * (voxelStackSize + 1));
	voxelVector.reserve(voxelStackSize * voxelStackSize * voxelStackSize);

	float scale = terrainSize / voxelStackSize;

	for (size_t depthIter = 0; depthIter < voxelStackSize; ++depthIter)			//	Depth
	{
		for (size_t heightIter = 0; heightIter < voxelStackSize; ++heightIter)	//	Row
		{
			for (size_t widthIter = 0; widthIter < voxelStackSize; ++widthIter)	//	Column
			{
				Voxel* newVoxel = new Voxel();
				Node* newNode = new Node();

				newVoxel->setPosition(position);

				if (depthIter == 0)	//	First depth
				{
					if (heightIter == 0)	//	Voxels in the first row
					{
						if (widthIter == 0)	//	First Voxel
						{
							//	Create 8 new node isntances

							XMFLOAT3 nodePosition = XMFLOAT3(position.x - (1.0f * scale), position.y - (1.0f * scale), position.z - (1.0f * scale));

							for (int nodeDepthIter = 0; nodeDepthIter < 2; ++nodeDepthIter)
							{
								for (int nodeHeightIter = 0; nodeHeightIter < 2; ++nodeHeightIter)
								{
									for (int nodeWidthIter = 0; nodeWidthIter < 2; ++nodeWidthIter)
									{
										newNode->density = RandomDensity();
										newNode->position = nodePosition;

										nodeVector.push_back(*newNode);
									
										nodePosition.z += (2.0f * scale);
									}
									nodePosition.y += (2.0f * scale);
									nodePosition.z = position.z - (1.0f * scale);
								}
								nodePosition.x += (2.0f * scale);
								nodePosition.y = position.y - (1.0f * scale);
							}

							//	Assign nodes to the voxel
							for (int nodeIter = 0; nodeIter < 8; ++nodeIter)
								newVoxel->setNode(nodeIter, nodeVector[nodeIter]);
						}
						else	//	Voxels after the first column still on the first row, copies 4 nodes from the previous voxel and creates 4 new instances
						{
							
							newVoxel->setNode(0, *voxelVector[voxelVector.size() - 1].getNode(1));	//	Node 0	

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x - (1.0f * scale), position.y - (1.0f * scale), position.z + (1.0f * scale));	//	Node 1
							nodeVector.push_back(*newNode);
							newVoxel->setNode(1, nodeVector[nodeVector.size() - 1]);

							newVoxel->setNode(2, *voxelVector[voxelVector.size() - 1].getNode(3));	//	Node 2

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x - (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));	//	Node 3
							nodeVector.push_back(*newNode);
							newVoxel->setNode(3, nodeVector[nodeVector.size() - 1]);

							newVoxel->setNode(4, *voxelVector[voxelVector.size() - 1].getNode(5));	//	Node 4

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y - (1.0f * scale), position.z + (1.0f * scale));	//	Node 5
							nodeVector.push_back(*newNode);
							newVoxel->setNode(5, nodeVector[nodeVector.size() - 1]);

							newVoxel->setNode(6, *voxelVector[voxelVector.size() - 1].getNode(7));	//	Node 6

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));	//	Node 7
							nodeVector.push_back(*newNode);
							newVoxel->setNode(7, nodeVector[nodeVector.size() - 1]);
						}
					}
					else	//	Voxels after the first row
					{
						if (widthIter == 0)	//	Voxel in the first column, reuses 4 nodes from voxel below it
						{
							//	3D arrays in 1D arrays are traversed by using (depth * width^2) + (width * row) + column

							float voxelStackSizeSquared = pow(voxelStackSize, 2.0f);

							newVoxel->setNode(0, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(2));
							newVoxel->setNode(1, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(3));

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x - (1.0f * scale), position.y + (1.0f * scale), position.z - (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(2, nodeVector[nodeVector.size() - 1]);

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x - (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(3, nodeVector[nodeVector.size() - 1]);

							newVoxel->setNode(4, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(6));
							newVoxel->setNode(5, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(7));

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z - (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(6, nodeVector[nodeVector.size() - 1]);

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(7, nodeVector[nodeVector.size() - 1]);
						}
						else	//	Voxels in subsequent columns, copies 6 nodes (4 from behind, 2 from below)
						{
							float voxelStackSizeSquared = pow(voxelStackSize, 2.0f);

							newVoxel->setNode(0, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(1));	//	Node 0 copies node 1 from the voxel behind
							newVoxel->setNode(1, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(3));	//	Node 1 copies node 3 from the voxel below
							newVoxel->setNode(2, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(3));	//	Node 2 copies node 3 from the voxel behind

							//	Create a new node instance
							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x - (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(3, nodeVector[nodeVector.size() - 1]);

							newVoxel->setNode(4, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(5));	//	Behind
							newVoxel->setNode(5, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(7));	//	Below
							newVoxel->setNode(6, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(7));	//	Behind

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(7, nodeVector[nodeVector.size() - 1]);
						}
					}
				}
				else	//	Voxels after the first depth iteration
				{
					if (heightIter == 0)	//	The first row of voxels
					{
						if (widthIter == 0)	//	The first column of voxels i.e. first voxel
						{
							//	Copies 4 nodes from the voxel in previous depth iteration at same row and column, creates 4 new node instances

							float voxelStackSizeSquared = pow(voxelStackSize, 2.0f);

							newVoxel->setNode(0, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(4));
							newVoxel->setNode(1, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(5));
							newVoxel->setNode(2, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(6));
							newVoxel->setNode(3, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(7));

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y - (1.0f * scale), position.z - (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(4, nodeVector[nodeVector.size() - 1]);

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y - (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(5, nodeVector[nodeVector.size() - 1]);

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z - (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(6, nodeVector[nodeVector.size() - 1]);

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(7, nodeVector[nodeVector.size() - 1]);
						}
						else	//	Subsequent columns of voxels
						{
							//	Copies 2 nodes from voxel in previous column and 4 from voxel in previous depth, creates 2 new node instances
							/*
							* 4 nodes copied from previous depth instead of 4 from previous
							* column because of continuity: (0,1,2,3 = 4,5,6,7) vs (0,2,4,6 = 1,3,5,7)
							*/

							float voxelStackSizeSquared = pow(voxelStackSize, 2.0f);

							newVoxel->setNode(0, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(4));	//	Previous depth
							newVoxel->setNode(1, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(5));
							newVoxel->setNode(2, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(6));
							newVoxel->setNode(3, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(7));

							newVoxel->setNode(4, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(5));	//	Behind

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y - (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(5, nodeVector[nodeVector.size() - 1]);

							newVoxel->setNode(6, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(7));	//	Behind

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(7, nodeVector[nodeVector.size() - 1]);
						}
					}
					else	//	Subsequent rows of voxels
					{
						if (widthIter == 0)	//	First column of voxels i.e. first voxel
						{
							//	Copies 4 nodes from voxel in previous depth and 2 nodes from the voxel in previous row, creates 2 new node instances

							float voxelStackSizeSquared = pow(voxelStackSize, 2.0f);

							newVoxel->setNode(0, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(4));	//	Previous depth
							newVoxel->setNode(1, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(5));
							newVoxel->setNode(2, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(6));
							newVoxel->setNode(3, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(7));

							newVoxel->setNode(4, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(6));	//	Below
							newVoxel->setNode(5, *voxelVector[((depthIter) * voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(7));	//	Below

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z - (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(6, nodeVector[nodeVector.size() - 1]);

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(7, nodeVector[nodeVector.size() - 1]);

						}
						else	//	Subsequent columns of voxels
						{
							//	Copies 4 nodes from voxel in previous depth, 2 nodes from voxel in previous column, and 1 node from voxel in previous row, creates 1 new node instance

							float voxelStackSizeSquared = pow(voxelStackSize, 2.0f);

							newVoxel->setNode(0, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(4));	//	Previous depth
							newVoxel->setNode(1, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(5));
							newVoxel->setNode(2, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(6));
							newVoxel->setNode(3, *voxelVector[((depthIter - 1) * voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter)].getNode(7));

							newVoxel->setNode(4, *voxelVector[((depthIter)* voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(5));	//	Behind
							newVoxel->setNode(5, *voxelVector[((depthIter)* voxelStackSizeSquared) + (voxelStackSize * (heightIter - 1)) + (widthIter)].getNode(7));	//	Below
							newVoxel->setNode(6, *voxelVector[((depthIter)* voxelStackSizeSquared) + (voxelStackSize * (heightIter)) + (widthIter - 1)].getNode(7));	//	Behind

							newNode->density = RandomDensity();
							newNode->position = XMFLOAT3(position.x + (1.0f * scale), position.y + (1.0f * scale), position.z + (1.0f * scale));
							nodeVector.push_back(*newNode);
							newVoxel->setNode(7, nodeVector[nodeVector.size() - 1]);

						}
					}
				}

				voxelVector.push_back(*newVoxel);

				delete newVoxel;
				delete newNode;

				position.z += (2.0f * scale);
			}
			position.y += (2.0f * scale);
			position.z = 0.0f;
		}
		position.x += (2.0f * scale);
		position.y = 0.0f;
	}

}

void App1::NoisifyDensities()
{
	//	Returns double between -1 and 1
	OpenSimplexNoise* simplexNoise = new OpenSimplexNoise(63129);

	std::vector<float> noises;
	float max = 50.0f, min = 50.0f;

	for (int i = 0; i < nodeVector.size(); ++i)
	{
		if (nodeVector[i].position.y / nodeVector.back().position.y * 100.0f <= surfaceLevel)
		{
			float newDensity = simplexNoise->Evaluate(nodeVector[i].position.x * noiseScale, nodeVector[i].position.y * noiseScale, nodeVector[i].position.z * noiseScale);
		
			noises.push_back(((newDensity)+1.0f) * 50.0f);

			if (noises[noises.size() - 1] > max)
				max = noises[noises.size() - 1];
			if (noises[noises.size() - 1] < min)
				min = noises[noises.size() - 1];

			newDensity = (newDensity + 1.0f) * 50.0f;	//	Scales between 0 and 100

			nodeVector[i].density = newDensity;
		}
		else
		{
			nodeVector[i].density = densityCutoff;
		}
	}

	delete simplexNoise;
	simplexNoise = 0;
}

float App1::RandomDensity()
{
	return  rand() % 100;
}

void App1::MarchingCubes()
{
	int vertCounter = 0;
	int triCounter = 0;
	XMFLOAT3* allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];	//	Maximum of 5 triangles per voxel. Each triangle needs 3 vertices.
	int* allIndices = new int[(5 * 3) * voxelVector.size()];

	if (useGPU)
	{
		the_clock::time_point start = the_clock::now();

		RunMarchingCubes(densityCutoff, voxelVector.data(), voxelVector.size(), allVertices, allIndices, vertCounter, *triTable, edgeTable, gpuNodeParseTime, gpuMallocTime, gpuMemcpyTime, gpuMarchTime, gpuCompactTime, gpuFreeTime);
		
		the_clock::time_point end = the_clock::now();
	
		time_taken = duration_cast<microseconds>(end - start).count();

	}
	else
	{
		the_clock::time_point start = the_clock::now();

		for (int voxelIter = 0; voxelIter < voxelVector.size(); ++voxelIter)
		{
			XMFLOAT3 vertexArray[12];
			unsigned char voxelByte = 0x00;
		
			if (voxelVector[voxelIter].getNode(1)->density < densityCutoff) voxelByte |= 1;
			if (voxelVector[voxelIter].getNode(5)->density < densityCutoff) voxelByte |= 2;	//	If density is below cutoff set bit number 2 to 1
			if (voxelVector[voxelIter].getNode(4)->density < densityCutoff) voxelByte |= 4;
			if (voxelVector[voxelIter].getNode(0)->density < densityCutoff) voxelByte |= 8;
			if (voxelVector[voxelIter].getNode(3)->density < densityCutoff) voxelByte |= 16;
			if (voxelVector[voxelIter].getNode(7)->density < densityCutoff) voxelByte |= 32;
			if (voxelVector[voxelIter].getNode(6)->density < densityCutoff) voxelByte |= 64;
			if (voxelVector[voxelIter].getNode(2)->density < densityCutoff) voxelByte |= 128;

			//	If voxel is entirely inside or outside of surface
			if (edgeTable[voxelByte] == 0)
				continue;

			//	If the binary value of element voxelByte in edgeTable has a 1 in the same position as 1=0x00000001 perform linear interpolation and store the value of that vertex
			//	Had to translate indexing convention for vertices so that they match with the edges
			if (edgeTable[voxelByte] & 1)	//	AND operator
				vertexArray[0] = VertexInterp(voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(5)->density);
			if (edgeTable[voxelByte] & 2)
				vertexArray[1] = VertexInterp(voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(4)->density);
			if (edgeTable[voxelByte] & 4)
				vertexArray[2] = VertexInterp(voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(0)->density);
			if (edgeTable[voxelByte] & 8)
				vertexArray[3] = VertexInterp(voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(1)->density);
			if (edgeTable[voxelByte] & 16)
				vertexArray[4] = VertexInterp(voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(3)->density, voxelVector[voxelIter].getNode(7)->density);
			if (edgeTable[voxelByte] & 32)
				vertexArray[5] = VertexInterp(voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(7)->density, voxelVector[voxelIter].getNode(6)->density);
			if (edgeTable[voxelByte] & 64)
				vertexArray[6] = VertexInterp(voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(6)->density, voxelVector[voxelIter].getNode(2)->density);
			if (edgeTable[voxelByte] & 128)
				vertexArray[7] = VertexInterp(voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(2)->density, voxelVector[voxelIter].getNode(3)->density);
			if (edgeTable[voxelByte] & 256)
				vertexArray[8] = VertexInterp(voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(3)->density);
			if (edgeTable[voxelByte] & 512)
				vertexArray[9] = VertexInterp(voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(7)->density);
			if (edgeTable[voxelByte] & 1024)
				vertexArray[10] = VertexInterp(voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(6)->density);
			if (edgeTable[voxelByte] & 2048)
				vertexArray[11] = VertexInterp(voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(2)->density);
		
			for (int i = 0; triTable[voxelByte][i] != -1; i += 3)
			{
				allVertices[(triCounter * 3)] =	vertexArray[triTable[voxelByte][i]];
				allVertices[(triCounter * 3) + 1] = vertexArray[triTable[voxelByte][i + 1]];
				allVertices[(triCounter * 3) + 2] = vertexArray[triTable[voxelByte][i + 2]];

				allIndices[(triCounter * 3)] = triCounter * 3;
				allIndices[(triCounter * 3) + 1] = (triCounter * 3) + 1;
				allIndices[(triCounter * 3) + 2] = (triCounter * 3) + 2;
			
				triCounter++;
				vertCounter += 3;
			}
		}

		the_clock::time_point end = the_clock::now();

		time_taken = duration_cast<microseconds>(end - start).count();
	}

	

	m_Isosurface->initBuffers(renderer->getDevice(), allVertices, allIndices, vertCounter);
	delete allVertices;
	delete allIndices;
	allVertices = 0;
	allIndices = 0;
}

XMFLOAT3 App1::VertexInterp(XMFLOAT3 p1, XMFLOAT3 p2, float valp1, float valp2)
{
	double mu;
	XMFLOAT3 p;

	if (abs(densityCutoff - valp1) < 0.00001)
		return p1;
	if (abs(densityCutoff - valp2) < 0.00001)
		return p2;
	if (abs(valp1 - valp2) < 0.00001)
		return p1;
	mu = (densityCutoff - valp1) / (valp2 - valp1);
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);

	return p;
}

XMFLOAT3 App1::VertexInterp(float otherDensityCut, XMFLOAT3 p1, XMFLOAT3 p2, float valp1, float valp2)
{
	double mu;
	XMFLOAT3 p;

	if (abs(otherDensityCut - valp1) < 0.00001)
		return p1;
	if (abs(otherDensityCut - valp2) < 0.00001)
		return p2;
	if (abs(valp1 - valp2) < 0.00001)
		return p1;
	mu = (otherDensityCut - valp1) / (valp2 - valp1);
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);

	return p;
}

bool App1::frame()
{
	bool result;

	result = BaseApplication::frame();
	if (!result)
	{
		return false;
	}
	// Render the graphics.
	result = render();
	if (!result)
	{
		return false;
	}

	return true;
}

bool App1::render()
{
	XMMATRIX worldMatrix, viewMatrix, projectionMatrix;

	// Clear the scene. (default blue colour)
	renderer->beginScene(0.39f, 0.58f, 0.92f, 1.0f);

	// Generate the view matrix based on the camera's position.
	camera->update();

	// Get the world, view, projection, and ortho matrices from the camera and Direct3D objects.
	worldMatrix = renderer->getWorldMatrix();
	viewMatrix = camera->getViewMatrix();
	projectionMatrix = renderer->getProjectionMatrix();

	m_Isosurface->sendData(renderer->getDeviceContext());
	m_BasicShader->setShaderParameters(renderer->getDeviceContext(), worldMatrix, viewMatrix, projectionMatrix, textureMgr->getTexture(L"brick"), light);
	m_BasicShader->render(renderer->getDeviceContext(), m_Isosurface->getIndexCount());

	// Render GUI
	gui();

	// Swap the buffers
	renderer->endScene();

	return true;
}

void App1::gui()
{
	// Force turn off unnecessary shader stages.
	renderer->getDeviceContext()->GSSetShader(NULL, NULL, 0);
	renderer->getDeviceContext()->HSSetShader(NULL, NULL, 0);
	renderer->getDeviceContext()->DSSetShader(NULL, NULL, 0);

	// Build UI
	ImGui::Text("FPS: %.2f", timer->getFPS());
	ImGui::Text( "Camera Pos: (%.2f, %.2f, %.2f)", camera->getPosition().x, camera->getPosition().y, camera->getPosition().z );
	ImGui::Checkbox("Wireframe mode", &wireframeToggle);
	
	ImGui::Checkbox("Use GPU", &useGPU);

	if(ImGui::DragFloat("Density Cutoff", &densityCutoff, 0.5f, 0.0f, 100.0f))
	{
		MarchingCubes();
	}
	if (ImGui::DragFloat("Scale", &noiseScale, 0.01f, 0.0f, 1.0f))
	{
		NoisifyDensities();
		MarchingCubes();
	}
	if (ImGui::DragFloat("TerrainSize", &terrainSize, 0.5f, 0.0f, 50.0f) ||
		ImGui::DragInt("StackSize", &voxelStackSize, 1, 0, 50))
	{
		InitialiseStack();
		NoisifyDensities();
		MarchingCubes();
	}
	if (ImGui::DragFloat("SurfaceLevel", &surfaceLevel, 0.5f, 0.0f, 100.0f))
	{
		NoisifyDensities();
		MarchingCubes();
	}

	ImGui::Text("Time taken: \t%f", time_taken);
	ImGui::Text("NodeParseGPU time: \t%f", gpuNodeParseTime);
	ImGui::Text("MallocGPU time: \t%f", gpuMallocTime);
	ImGui::Text("MemcpyGPU time: \t%f", gpuMemcpyTime);
	ImGui::Text("MarchGPU time: \t%f", gpuMarchTime);
	ImGui::Text("CompactGPU time: \t%f", gpuCompactTime);
	ImGui::Text("FreeGPU time: \t%f", gpuFreeTime);


	ImGui::NewLine();
	ImGui::Separator();

	//
	//	Testing Variables
	ImGui::Text("Constant Test Variables");
	ImGui::InputInt("Iterations", &testIterations);
	ImGui::InputFloat("StepSize", &testStepSize);
	ImGui::Checkbox("Use GPU?", &testUseGPU);
	ImGui::InputFloat("Density", &testDensity);
	ImGui::InputInt("Stack Size", &testStackSize);
	ImGui::InputFloat("Surface Y", &testSurface);
	ImGui::InputFloat("Noise Scale", &testScale);

	ImGui::NewLine();
	ImGui::Separator();

	//
	//	Testing Ranges and Buttons
	ImGui::Text("Test Ranges");

	ImGui::InputFloat2("Density Range", &testDensityRange.x);
	if (ImGui::Button("Test Densities"))
		DensityTest(testIterations, testStepSize, testUseGPU, testDensityRange, testStackSize, testSurface, testScale);
	
	ImGui::NewLine();

	ImGui::InputInt2("Stack Range", &testStackSizeRange.x);
	if (ImGui::Button("Test Stack Sizes"))
		StackSizeTest(testIterations, testStepSize, testUseGPU, testDensity, testStackSizeRange, testSurface,  testScale);

	ImGui::NewLine();
	
	ImGui::InputFloat2("Surface Range", &testSurfaceRange.x);
	if (ImGui::Button("Test Surface Levels"))
		SurfaceTest(testIterations, testStepSize, testUseGPU, testDensity, testStackSize, testSurfaceRange, testScale);

	ImGui::NewLine();

	ImGui::InputFloat2("Scale Range", &testScaleRange.x);
	if (ImGui::Button("Test Scales"))
		ScaleTest(testIterations, testStepSize, testUseGPU, testDensity, testStackSize, testSurface, testScaleRange);

	ImGui::Render();
	ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
}

void App1::DensityTest(unsigned int iterations,
	float stepSize,
	bool useGPU,
	XMFLOAT2 mcDensityRange,
	unsigned int stackSize,
	float mcSurface,
	float mcScale) 
{
	std::ofstream myFile;
	double test_time;

	XMFLOAT3* allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];	//	Maximum of 5 triangles per voxel. Each triangle needs 3 vertices.
	int* allIndices = new int[(5 * 3) * voxelVector.size()];
	int vertCounter = 0;

	voxelStackSize = stackSize;
	surfaceLevel = mcSurface;
	noiseScale = mcScale;
	densityCutoff = mcDensityRange.x;

	InitialiseStack();
	NoisifyDensities();

	the_clock::time_point start, end;
	std::string fileString = "DensityTest ";
	useGPU ? fileString += "GPU " : fileString += "CPU ";
	fileString += std::to_string(mcDensityRange.x) + " " + std::to_string(mcDensityRange.y) + " " + std::to_string(iterations) + ".csv";

	myFile.open(fileString.c_str());

	myFile << "Number of Rows," << iterations * int(((mcDensityRange.y + stepSize) - mcDensityRange.x) / stepSize) << std::endl;
	myFile << "Iterations," << iterations << std::endl;
	myFile << "Step Size," << stepSize << std::endl;
	myFile << "Used GPU?," << useGPU << std::endl;
	myFile << "Density Range," << mcDensityRange.x << "," << mcDensityRange.y << std::endl;
	myFile << "Stack Size," << stackSize << std::endl;
	myFile << "Surface Level," << mcSurface << std::endl;
	myFile << "Scale," << mcScale << std::endl;

	if (useGPU)
	{
		myFile << "Density Cutoff,Iteration,Node Parse,CudaMalloc Time,CudaMemcpy Time,March Time,Compact Time,Free Time,Total Time" << std::endl;

		for (float densityIter = mcDensityRange.x; densityIter <= mcDensityRange.y; densityIter += stepSize)
		{
			for (int i = 0; i < iterations; ++i)
			{
				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;

				start = the_clock::now();

				RunMarchingCubes(densityIter, voxelVector.data(), voxelVector.size(), allVertices, allIndices, vertCounter, *triTable, edgeTable, gpuNodeParseTime, gpuMallocTime, gpuMemcpyTime, gpuMarchTime, gpuCompactTime, gpuFreeTime);

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();


				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << densityIter << "," << i << "," << gpuNodeParseTime << "," << gpuMallocTime << "," << gpuMemcpyTime << "," << gpuMarchTime << "," << gpuCompactTime << "," << gpuFreeTime << "," << test_time << std::endl;
			}
		}
	}
	else
	{
		myFile << "Density Cutoff,Iteration,Time" << std::endl;

		int triCounter = 0;
		XMFLOAT3 vertexArray[12];
		unsigned char voxelByte = 0x00;

		for (float densityIter = mcDensityRange.x; densityIter <= mcDensityRange.y; densityIter += stepSize)
		{
			for (int i = 0; i < iterations; ++i)
			{

				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;
				triCounter = 0;

				start = the_clock::now();

				for (int voxelIter = 0; voxelIter < voxelVector.size(); ++voxelIter)
				{
					voxelByte = 0x00;

					if (voxelVector[voxelIter].getNode(1)->density < densityIter) voxelByte |= 1;
					if (voxelVector[voxelIter].getNode(5)->density < densityIter) voxelByte |= 2;	//	If density is below cutoff set bit number 2 to 1
					if (voxelVector[voxelIter].getNode(4)->density < densityIter) voxelByte |= 4;
					if (voxelVector[voxelIter].getNode(0)->density < densityIter) voxelByte |= 8;
					if (voxelVector[voxelIter].getNode(3)->density < densityIter) voxelByte |= 16;
					if (voxelVector[voxelIter].getNode(7)->density < densityIter) voxelByte |= 32;
					if (voxelVector[voxelIter].getNode(6)->density < densityIter) voxelByte |= 64;
					if (voxelVector[voxelIter].getNode(2)->density < densityIter) voxelByte |= 128;

					//	If voxel is entirely inside or outside of surface
					if (edgeTable[voxelByte] == 0)
						continue;

					//	If the binary value of element voxelByte in edgeTable has a 1 in the same position as 1=0x00000001 perform linear interpolation and store the value of that vertex
					//	Had to translate indexing convention for vertices so that they match with the edges
					if (edgeTable[voxelByte] & 1)	//	AND operator
						vertexArray[0] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(5)->density);
					if (edgeTable[voxelByte] & 2)
						vertexArray[1] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(4)->density);
					if (edgeTable[voxelByte] & 4)
						vertexArray[2] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(0)->density);
					if (edgeTable[voxelByte] & 8)
						vertexArray[3] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(1)->density);
					if (edgeTable[voxelByte] & 16)
						vertexArray[4] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(3)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 32)
						vertexArray[5] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(7)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 64)
						vertexArray[6] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(6)->density, voxelVector[voxelIter].getNode(2)->density);
					if (edgeTable[voxelByte] & 128)
						vertexArray[7] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(2)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 256)
						vertexArray[8] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 512)
						vertexArray[9] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 1024)
						vertexArray[10] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 2048)
						vertexArray[11] = VertexInterp(densityIter, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(2)->density);

					for (int i = 0; triTable[voxelByte][i] != -1; i += 3)
					{
						allVertices[(triCounter * 3)] = vertexArray[triTable[voxelByte][i]];
						allVertices[(triCounter * 3) + 1] = vertexArray[triTable[voxelByte][i + 1]];
						allVertices[(triCounter * 3) + 2] = vertexArray[triTable[voxelByte][i + 2]];

						allIndices[(triCounter * 3)] = triCounter * 3;
						allIndices[(triCounter * 3) + 1] = (triCounter * 3) + 1;
						allIndices[(triCounter * 3) + 2] = (triCounter * 3) + 2;

						triCounter++;
						vertCounter += 3;
					}
				}

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();
				
				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << densityIter << "," << i << "," << test_time << std::endl;
			}
		}
	}

	myFile.close();

	if (allVertices)
	{
		delete allVertices;
		allVertices = 0;
	}
	if (allIndices)
	{
		delete allIndices;
		allIndices = 0;
	}
}

void App1::StackSizeTest(unsigned int iterations,
	int stepSize,
	bool useGPU,
	float mcDensity,
	XMINT2 stackSizeRange,
	float mcSurface,
	float mcScale)
{
	std::ofstream myFile;
	double test_time;

	XMFLOAT3* allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];	//	Maximum of 5 triangles per voxel. Each triangle needs 3 vertices.
	int* allIndices = new int[(5 * 3) * voxelVector.size()];
	int vertCounter = 0;

	the_clock::time_point start, end;
	std::string fileString = "StackSizeTest ";
	useGPU ? fileString += "GPU " : fileString += "CPU ";
	fileString += std::to_string(stackSizeRange.x) + " " + std::to_string(stackSizeRange.y) + " " + std::to_string(iterations) + ".csv";

	myFile.open(fileString.c_str());

	myFile << "Number of Rows," << iterations * int(((stackSizeRange.y + stepSize) - stackSizeRange.x) / stepSize) << std::endl;
	myFile << "Iterations," << iterations << std::endl;
	myFile << "Step Size," << stepSize << std::endl;
	myFile << "Used GPU?," << useGPU << std::endl;
	myFile << "Density Cutoff," << mcDensity << std::endl;
	myFile << "Stack Size Range," << stackSizeRange.x << "," << stackSizeRange.y << std::endl;
	myFile << "Surface Level," << mcSurface << std::endl;
	myFile << "Scale," << mcScale << std::endl;

	if (useGPU)
	{
		myFile << "Stack Size,Iteration,Node Parse,CudaMalloc Time,CudaMemcpy Time,March Time,Compact Time,Free Time,Total Time" << std::endl;

		//
		//	Init voxel vars
		densityCutoff = mcDensity;
		surfaceLevel = mcSurface;
		noiseScale = mcScale;

		for (int stackSizeIter = stackSizeRange.x; stackSizeIter <= stackSizeRange.y; stackSizeIter += stepSize)
		{
			voxelStackSize = stackSizeIter;

			InitialiseStack();
			NoisifyDensities();

			for (int i = 0; i < iterations; ++i)
			{
				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;

				start = the_clock::now();

				RunMarchingCubes(mcDensity, voxelVector.data(), voxelVector.size(), allVertices, allIndices, vertCounter, *triTable, edgeTable, gpuNodeParseTime, gpuMallocTime, gpuMemcpyTime, gpuMarchTime, gpuCompactTime, gpuFreeTime);

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();


				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << stackSizeIter << "," << i << "," << gpuNodeParseTime << "," << gpuMallocTime << "," << gpuMemcpyTime << "," << gpuMarchTime << "," << gpuCompactTime << "," << gpuFreeTime << "," << test_time << std::endl;
			}
		}
	}
	else
	{
		myFile << "Stack Size,Iteration,Time" << std::endl;

		int triCounter = 0;
		XMFLOAT3 vertexArray[12];
		unsigned char voxelByte = 0x00;

		//
		//	Init voxel vars
		densityCutoff = mcDensity;
		surfaceLevel = mcSurface;
		noiseScale = mcScale;

		for (int stackSizeIter = stackSizeRange.x; stackSizeIter <= stackSizeRange.y; stackSizeIter += stepSize)
		{
			voxelStackSize = stackSizeIter;

			InitialiseStack();
			NoisifyDensities();

			for (int i = 0; i < iterations; ++i)
			{
				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;
				triCounter = 0;

				start = the_clock::now();

				for (int voxelIter = 0; voxelIter < voxelVector.size(); ++voxelIter)
				{
					voxelByte = 0x00;

					if (voxelVector[voxelIter].getNode(1)->density < mcDensity) voxelByte |= 1;
					if (voxelVector[voxelIter].getNode(5)->density < mcDensity) voxelByte |= 2;	//	If density is below cutoff set bit number 2 to 1
					if (voxelVector[voxelIter].getNode(4)->density < mcDensity) voxelByte |= 4;
					if (voxelVector[voxelIter].getNode(0)->density < mcDensity) voxelByte |= 8;
					if (voxelVector[voxelIter].getNode(3)->density < mcDensity) voxelByte |= 16;
					if (voxelVector[voxelIter].getNode(7)->density < mcDensity) voxelByte |= 32;
					if (voxelVector[voxelIter].getNode(6)->density < mcDensity) voxelByte |= 64;
					if (voxelVector[voxelIter].getNode(2)->density < mcDensity) voxelByte |= 128;

					//	If voxel is entirely inside or outside of surface
					if (edgeTable[voxelByte] == 0)
						continue;

					//	If the binary value of element voxelByte in edgeTable has a 1 in the same position as 1=0x00000001 perform linear interpolation and store the value of that vertex
					//	Had to translate indexing convention for vertices so that they match with the edges
					if (edgeTable[voxelByte] & 1)	//	AND operator
						vertexArray[0] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(5)->density);
					if (edgeTable[voxelByte] & 2)
						vertexArray[1] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(4)->density);
					if (edgeTable[voxelByte] & 4)
						vertexArray[2] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(0)->density);
					if (edgeTable[voxelByte] & 8)
						vertexArray[3] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(1)->density);
					if (edgeTable[voxelByte] & 16)
						vertexArray[4] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(3)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 32)
						vertexArray[5] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(7)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 64)
						vertexArray[6] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(6)->density, voxelVector[voxelIter].getNode(2)->density);
					if (edgeTable[voxelByte] & 128)
						vertexArray[7] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(2)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 256)
						vertexArray[8] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 512)
						vertexArray[9] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 1024)
						vertexArray[10] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 2048)
						vertexArray[11] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(2)->density);

					for (int i = 0; triTable[voxelByte][i] != -1; i += 3)
					{
						allVertices[(triCounter * 3)] = vertexArray[triTable[voxelByte][i]];
						allVertices[(triCounter * 3) + 1] = vertexArray[triTable[voxelByte][i + 1]];
						allVertices[(triCounter * 3) + 2] = vertexArray[triTable[voxelByte][i + 2]];

						allIndices[(triCounter * 3)] = triCounter * 3;
						allIndices[(triCounter * 3) + 1] = (triCounter * 3) + 1;
						allIndices[(triCounter * 3) + 2] = (triCounter * 3) + 2;

						triCounter++;
						vertCounter += 3;
					}
				}

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();

				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << stackSizeIter << "," << i << "," << test_time << std::endl;
			}
		}
	}

	myFile.close();

	if (allVertices)
	{
		delete allVertices;
		allVertices = 0;
	}
	if (allIndices)
	{
		delete allIndices;
		allIndices = 0;
	}
}

void App1::SurfaceTest(unsigned int iterations,
	float stepSize,
	bool useGPU,
	float mcDensity,
	unsigned int stackSize,
	XMFLOAT2 mcSurfaceRange,
	float mcScale)
{
	std::ofstream myFile;
	double test_time;

	XMFLOAT3* allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];	//	Maximum of 5 triangles per voxel. Each triangle needs 3 vertices.
	int* allIndices = new int[(5 * 3) * voxelVector.size()];
	int vertCounter = 0;

	//
	//	Init voxel vars
	voxelStackSize = stackSize;
	noiseScale = mcScale;
	densityCutoff = mcDensity;

	InitialiseStack();

	the_clock::time_point start, end;
	std::string fileString = "SurfaceTest ";
	useGPU ? fileString += "GPU " : fileString += "CPU ";
	fileString += std::to_string(mcSurfaceRange.x) + " " + std::to_string(mcSurfaceRange.y) + " " + std::to_string(iterations) + ".csv";

	myFile.open(fileString.c_str());

	myFile << "Number of Rows," << iterations * int(((mcSurfaceRange.y + stepSize) - mcSurfaceRange.x) / stepSize) << std::endl;
	myFile << "Iterations," << iterations << std::endl;
	myFile << "Step Size," << stepSize << std::endl;
	myFile << "Used GPU?," << useGPU << std::endl;
	myFile << "Density Cutoff," << mcDensity << std::endl;
	myFile << "Stack Size," << stackSize << std::endl;
	myFile << "Surface Range," << mcSurfaceRange.x << "," << mcSurfaceRange.y << std::endl;
	myFile << "Scale," << mcScale << std::endl;

	if (useGPU)
	{
		myFile << "Surface Level,Iteration,Node Parse,CudaMalloc Time,CudaMemcpy Time,March Time,Compact Time,Free Time,Total Time" << std::endl;


		for (float surfaceIter = mcSurfaceRange.x; surfaceIter <= mcSurfaceRange.y; surfaceIter += stepSize)
		{
			surfaceLevel = surfaceIter;

			NoisifyDensities();


			for (int i = 0; i < iterations; ++i)
			{
				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;

				start = the_clock::now();

				RunMarchingCubes(mcDensity, voxelVector.data(), voxelVector.size(), allVertices, allIndices, vertCounter, *triTable, edgeTable, gpuNodeParseTime, gpuMallocTime, gpuMemcpyTime, gpuMarchTime, gpuCompactTime, gpuFreeTime);

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();


				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << surfaceIter << "," << i << "," << gpuNodeParseTime << "," << gpuMallocTime << "," << gpuMemcpyTime << "," << gpuMarchTime << "," << gpuCompactTime << "," << gpuFreeTime << "," << test_time << std::endl;
			}
		}
	}
	else
	{
		myFile << "Surface Level,Iteration,Time" << std::endl;

		int triCounter = 0;
		XMFLOAT3 vertexArray[12];
		unsigned char voxelByte = 0x00;

		for (float surfaceIter = mcSurfaceRange.x; surfaceIter <= mcSurfaceRange.y; surfaceIter += stepSize)
		{
			surfaceLevel = surfaceIter;

			NoisifyDensities();


			for (int i = 0; i < iterations; ++i)
			{
				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;
				triCounter = 0;

				start = the_clock::now();

				for (int voxelIter = 0; voxelIter < voxelVector.size(); ++voxelIter)
				{
					voxelByte = 0x00;

					if (voxelVector[voxelIter].getNode(1)->density < mcDensity) voxelByte |= 1;
					if (voxelVector[voxelIter].getNode(5)->density < mcDensity) voxelByte |= 2;	//	If density is below cutoff set bit number 2 to 1
					if (voxelVector[voxelIter].getNode(4)->density < mcDensity) voxelByte |= 4;
					if (voxelVector[voxelIter].getNode(0)->density < mcDensity) voxelByte |= 8;
					if (voxelVector[voxelIter].getNode(3)->density < mcDensity) voxelByte |= 16;
					if (voxelVector[voxelIter].getNode(7)->density < mcDensity) voxelByte |= 32;
					if (voxelVector[voxelIter].getNode(6)->density < mcDensity) voxelByte |= 64;
					if (voxelVector[voxelIter].getNode(2)->density < mcDensity) voxelByte |= 128;

					//	If voxel is entirely inside or outside of surface
					if (edgeTable[voxelByte] == 0)
						continue;

					//	If the binary value of element voxelByte in edgeTable has a 1 in the same position as 1=0x00000001 perform linear interpolation and store the value of that vertex
					//	Had to translate indexing convention for vertices so that they match with the edges
					if (edgeTable[voxelByte] & 1)	//	AND operator
						vertexArray[0] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(5)->density);
					if (edgeTable[voxelByte] & 2)
						vertexArray[1] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(4)->density);
					if (edgeTable[voxelByte] & 4)
						vertexArray[2] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(0)->density);
					if (edgeTable[voxelByte] & 8)
						vertexArray[3] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(1)->density);
					if (edgeTable[voxelByte] & 16)
						vertexArray[4] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(3)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 32)
						vertexArray[5] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(7)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 64)
						vertexArray[6] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(6)->density, voxelVector[voxelIter].getNode(2)->density);
					if (edgeTable[voxelByte] & 128)
						vertexArray[7] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(2)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 256)
						vertexArray[8] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 512)
						vertexArray[9] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 1024)
						vertexArray[10] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 2048)
						vertexArray[11] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(2)->density);

					for (int i = 0; triTable[voxelByte][i] != -1; i += 3)
					{
						allVertices[(triCounter * 3)] = vertexArray[triTable[voxelByte][i]];
						allVertices[(triCounter * 3) + 1] = vertexArray[triTable[voxelByte][i + 1]];
						allVertices[(triCounter * 3) + 2] = vertexArray[triTable[voxelByte][i + 2]];

						allIndices[(triCounter * 3)] = triCounter * 3;
						allIndices[(triCounter * 3) + 1] = (triCounter * 3) + 1;
						allIndices[(triCounter * 3) + 2] = (triCounter * 3) + 2;

						triCounter++;
						vertCounter += 3;
					}
				}

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();

				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << surfaceIter << "," << i << "," << test_time << std::endl;
			}
		}
	}

	myFile.close();

	if (allVertices)
	{
		delete allVertices;
		allVertices = 0;
	}
	if (allIndices)
	{
		delete allIndices;
		allIndices = 0;
	}
}

void App1::ScaleTest(unsigned int iterations,
	float stepSize,
	bool useGPU,
	float mcDensity,
	unsigned int stackSize,
	float mcSurface,
	XMFLOAT2 mcScaleRange)
{
	std::ofstream myFile;
	double test_time;

	XMFLOAT3* allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];	//	Maximum of 5 triangles per voxel. Each triangle needs 3 vertices.
	int* allIndices = new int[(5 * 3) * voxelVector.size()];
	int vertCounter = 0;

	//
	//	Init voxel vars
	voxelStackSize = stackSize;
	densityCutoff = mcDensity;
	surfaceLevel = mcSurface;

	InitialiseStack();

	the_clock::time_point start, end;
	std::string fileString = "ScaleTest ";
	useGPU ? fileString += "GPU " : fileString += "CPU ";
	fileString += std::to_string(mcScaleRange.x) + " " + std::to_string(mcScaleRange.y) + " " + std::to_string(iterations) + ".csv";

	myFile.open(fileString.c_str());

	myFile << "Number of Rows," << iterations * int(((mcScaleRange.y + stepSize) - mcScaleRange.x) / stepSize) << std::endl;
	myFile << "Iterations," << iterations << std::endl;
	myFile << "Step Size," << stepSize << std::endl;
	myFile << "Used GPU?," << useGPU << std::endl;
	myFile << "Density Cutoff," << mcDensity << std::endl;
	myFile << "Stack Size," << stackSize << std::endl;
	myFile << "Surface," << mcSurface <<std::endl;
	myFile << "Scale Range,"  <<  mcScaleRange.x << "," << mcScaleRange.y << std::endl;

	if (useGPU)
	{
		myFile << "Noise Scale,Iteration,Node Parse,CudaMalloc Time,CudaMemcpy Time,March Time,Compact Time,Free Time,Total Time" << std::endl;


		for (float scaleIter = mcScaleRange.x; scaleIter <= mcScaleRange.y; scaleIter += stepSize)
		{
			noiseScale = scaleIter;

			NoisifyDensities();


			for (int i = 0; i < iterations; ++i)
			{
				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;

				start = the_clock::now();

				RunMarchingCubes(mcDensity, voxelVector.data(), voxelVector.size(), allVertices, allIndices, vertCounter, *triTable, edgeTable, gpuNodeParseTime, gpuMallocTime, gpuMemcpyTime, gpuMarchTime, gpuCompactTime, gpuFreeTime);

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();


				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << scaleIter << "," << i << "," << gpuNodeParseTime << "," << gpuMallocTime << "," << gpuMemcpyTime << "," << gpuMarchTime << "," << gpuCompactTime << "," << gpuFreeTime << "," << test_time << std::endl;
			}
		}
	}
	else
	{
		myFile << "Surface Level,Iteration,Time" << std::endl;

		int triCounter = 0;
		XMFLOAT3 vertexArray[12];
		unsigned char voxelByte = 0x00;

		for (float scaleIter = mcScaleRange.x; scaleIter <= mcScaleRange.y; scaleIter += stepSize)
		{
			noiseScale = scaleIter;

			NoisifyDensities();


			for (int i = 0; i < iterations; ++i)
			{
				allVertices = new XMFLOAT3[(5 * 3) * voxelVector.size()];
				allIndices = new int[(5 * 3) * voxelVector.size()];
				vertCounter = 0;
				triCounter = 0;

				start = the_clock::now();

				for (int voxelIter = 0; voxelIter < voxelVector.size(); ++voxelIter)
				{
					voxelByte = 0x00;

					if (voxelVector[voxelIter].getNode(1)->density < mcDensity) voxelByte |= 1;
					if (voxelVector[voxelIter].getNode(5)->density < mcDensity) voxelByte |= 2;	//	If density is below cutoff set bit number 2 to 1
					if (voxelVector[voxelIter].getNode(4)->density < mcDensity) voxelByte |= 4;
					if (voxelVector[voxelIter].getNode(0)->density < mcDensity) voxelByte |= 8;
					if (voxelVector[voxelIter].getNode(3)->density < mcDensity) voxelByte |= 16;
					if (voxelVector[voxelIter].getNode(7)->density < mcDensity) voxelByte |= 32;
					if (voxelVector[voxelIter].getNode(6)->density < mcDensity) voxelByte |= 64;
					if (voxelVector[voxelIter].getNode(2)->density < mcDensity) voxelByte |= 128;

					//	If voxel is entirely inside or outside of surface
					if (edgeTable[voxelByte] == 0)
						continue;

					//	If the binary value of element voxelByte in edgeTable has a 1 in the same position as 1=0x00000001 perform linear interpolation and store the value of that vertex
					//	Had to translate indexing convention for vertices so that they match with the edges
					if (edgeTable[voxelByte] & 1)	//	AND operator
						vertexArray[0] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(5)->density);
					if (edgeTable[voxelByte] & 2)
						vertexArray[1] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(4)->density);
					if (edgeTable[voxelByte] & 4)
						vertexArray[2] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(0)->density);
					if (edgeTable[voxelByte] & 8)
						vertexArray[3] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(1)->density);
					if (edgeTable[voxelByte] & 16)
						vertexArray[4] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(3)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 32)
						vertexArray[5] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(7)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 64)
						vertexArray[6] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(6)->density, voxelVector[voxelIter].getNode(2)->density);
					if (edgeTable[voxelByte] & 128)
						vertexArray[7] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(2)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 256)
						vertexArray[8] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(1)->position, voxelVector[voxelIter].getNode(3)->position, voxelVector[voxelIter].getNode(1)->density, voxelVector[voxelIter].getNode(3)->density);
					if (edgeTable[voxelByte] & 512)
						vertexArray[9] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(5)->position, voxelVector[voxelIter].getNode(7)->position, voxelVector[voxelIter].getNode(5)->density, voxelVector[voxelIter].getNode(7)->density);
					if (edgeTable[voxelByte] & 1024)
						vertexArray[10] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(4)->position, voxelVector[voxelIter].getNode(6)->position, voxelVector[voxelIter].getNode(4)->density, voxelVector[voxelIter].getNode(6)->density);
					if (edgeTable[voxelByte] & 2048)
						vertexArray[11] = VertexInterp(mcDensity, voxelVector[voxelIter].getNode(0)->position, voxelVector[voxelIter].getNode(2)->position, voxelVector[voxelIter].getNode(0)->density, voxelVector[voxelIter].getNode(2)->density);

					for (int i = 0; triTable[voxelByte][i] != -1; i += 3)
					{
						allVertices[(triCounter * 3)] = vertexArray[triTable[voxelByte][i]];
						allVertices[(triCounter * 3) + 1] = vertexArray[triTable[voxelByte][i + 1]];
						allVertices[(triCounter * 3) + 2] = vertexArray[triTable[voxelByte][i + 2]];

						allIndices[(triCounter * 3)] = triCounter * 3;
						allIndices[(triCounter * 3) + 1] = (triCounter * 3) + 1;
						allIndices[(triCounter * 3) + 2] = (triCounter * 3) + 2;

						triCounter++;
						vertCounter += 3;
					}
				}

				end = the_clock::now();

				test_time = duration_cast<microseconds>(end - start).count();

				delete allVertices;
				allVertices = 0;
				delete allIndices;
				allIndices = 0;

				myFile << scaleIter << "," << i << "," << test_time << std::endl;
			}
		}
	}

	myFile.close();

	if (allVertices)
	{
		delete allVertices;
		allVertices = 0;
	}
	if (allIndices)
	{
		delete allIndices;
		allIndices = 0;
	}
}