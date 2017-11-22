//
//  CompFab.cpp
//  voxelizer
//
//
//

#include "../include/CompFab.h"
#include <iostream>

using namespace CompFab;





CompFab::Vec3Struct::Vec3Struct()
{
    m_x = m_y = m_z = 0.0;
}

CompFab::Vec3Struct::Vec3Struct(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

void CompFab::Vec3Struct::normalize() {
    
    double magnitude = sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
    
    if(magnitude > EPSILON)
    {
        m_x /= magnitude;
        m_y /= magnitude;
        m_z /= magnitude;
    }
}

//Data Types
CompFab::Vec3iStruct::Vec3iStruct()
{
    m_x = m_y = m_z = 0.0;
}

CompFab::Vec3iStruct::Vec3iStruct(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

CompFab::Vec2fStruct::Vec2fStruct()
{
    m_x = m_y = 0.0;
}

CompFab::Vec2fStruct::Vec2fStruct(double x, double y)
{
    m_x = x;
    m_y = y;
}

CompFab::RayStruct::RayStruct()
{
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_direction[0] = 1.0;
    m_direction[1] = m_direction[2] = 0.0;
}

CompFab::RayStruct::RayStruct(Vec3 &origin, Vec3 &direction)
{
    m_origin = origin;
    m_direction = direction;
}

CompFab::TriangleStruct::TriangleStruct(Vec3 &v1, Vec3 &v2,Vec3 &v3)
{
    m_v1 = v1;
    m_v2 = v2;
    m_v3 = v3;
}

CompFab::Vec3 CompFab::operator-(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[0] - v2[0];
    v3[1] = v1[1] - v2[1];
    v3[2] = v1[2] - v2[2];

    return v3;
}

CompFab::Vec3 CompFab::operator+(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[0] + v2[0];
    v3[1] = v1[1] + v2[1];
    v3[2] = v1[2] + v2[2];
    
    return v3;
}


//Cross Product
Vec3 CompFab::operator%(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return v3;
}

//Dot Product
double CompFab::operator*(const Vec3 &v1, const Vec3 &v2)
{
    return v1.m_x*v2.m_x + v1.m_y*v2.m_y+v1.m_z*v2.m_z;
}


//Grid structure for Voxels
CompFab::VoxelGridStruct::VoxelGridStruct(Vec3 lowerLeft, unsigned int dimX, unsigned int dimY, unsigned int dimZ, double spacing)
{
    m_lowerLeft = lowerLeft;
    m_dimX = dimX;
    m_dimY = dimY;
    m_dimZ = dimZ;
    m_size = dimX*dimY*dimZ;
    m_spacing = spacing;
    
    //Allocate Memory
    m_insideArray = new bool[m_size];
    
    for(unsigned int ii=0; ii<m_size; ++ii)
    {
        m_insideArray[ii] = false;
    }
    
}

CompFab::VoxelGridStruct::~VoxelGridStruct()
{
    delete[] m_insideArray;
}

MeshTree* CompFab::MeshTree::build(std::vector<Triangle> meshList)
{       
    bbox = BoundingBox();
    bbox.x_min = bbox.x_max = meshList[0].m_v1.m_x;
    bbox.y_min = bbox.y_max = meshList[0].m_v1.m_y;
    bbox.z_min = bbox.z_max = meshList[0].m_v1.m_z;
    for (unsigned int i=0; i < meshList.size(); i++) {
        double x_min_new = std::min(meshList[i].m_v1.m_x, 
                           std::min(meshList[i].m_v2.m_x,
                                    meshList[i].m_v3.m_x));
        double x_max_new = std::max(meshList[i].m_v1.m_x, 
                           std::max(meshList[i].m_v2.m_x,
                                    meshList[i].m_v3.m_x));
        double y_min_new = std::min(meshList[i].m_v1.m_y, 
                           std::min(meshList[i].m_v2.m_y,
                                    meshList[i].m_v3.m_y));
        double y_max_new = std::max(meshList[i].m_v1.m_y, 
                           std::max(meshList[i].m_v2.m_y,
                                    meshList[i].m_v3.m_y));
        double z_min_new = std::min(meshList[i].m_v1.m_z, 
                           std::min(meshList[i].m_v2.m_z,
                                    meshList[i].m_v3.m_z));
        double z_max_new = std::max(meshList[i].m_v1.m_z, 
                           std::max(meshList[i].m_v2.m_z,
                                    meshList[i].m_v3.m_z));

        bbox.update(x_min_new, y_min_new, z_min_new, 
                    x_max_new, y_max_new, z_max_new);
    }

    return build(meshList, bbox);
}


MeshTree* CompFab::MeshTree::build(std::vector<Triangle> meshList, BoundingBox bbox)
{    
    MeshTree* meshTree = new MeshTree();
    meshTree->left = NULL;
    meshTree->right = NULL;
    meshTree->bbox = bbox;
    meshTree->meshs = meshList; 
    if (meshList.size() == 0) {
        meshTree->isLeaf = true;
        return meshTree;
    }
    else if (meshList.size() <= 10) {
        meshTree->isLeaf = true;
        return meshTree;
    }
    else {
        meshTree->isLeaf = false;
        int axisNum = bbox.maxAxis();
        double splitPoint;
        BoundingBox bboxLeft = bbox, bboxRight = bbox; //might have problem
        std::vector<Triangle> leftList, rightList;
        if (axisNum == 1) {
            splitPoint = (bbox.x_max + bbox.x_min) / 2;
            for (unsigned int i=0; i < meshList.size(); i++) {
                if ((meshList[i].m_v1.m_x + meshList[i].m_v2.m_x + meshList[i].m_v3.m_x) / 3 <= splitPoint) {
                    leftList.push_back(meshList[i]);
                }
                else {
                    rightList.push_back(meshList[i]);
                }
            }
            bboxLeft.x_max = bboxRight.x_min = splitPoint;
        }
        else if(axisNum == 2) {
            splitPoint = (bbox.y_max + bbox.y_min) / 2;
            for (unsigned int i=0; i < meshList.size(); i++) {
                if ((meshList[i].m_v1.m_y + meshList[i].m_v2.m_y + meshList[i].m_v3.m_y) / 3 <= splitPoint) {
                    leftList.push_back(meshList[i]);
                }
                else {
                    rightList.push_back(meshList[i]);
                }
            }
            bboxLeft.y_max = bboxRight.y_min = splitPoint;
        }
        else {
            splitPoint = (bbox.z_max + bbox.z_min) / 2;
            for (unsigned int i=0; i < meshList.size(); i++) {
                if ((meshList[i].m_v1.m_z + meshList[i].m_v2.m_z + meshList[i].m_v3.m_z) / 3 <= splitPoint) {
                    leftList.push_back(meshList[i]);
                }
                else {
                    rightList.push_back(meshList[i]);
                }
            }
            bboxLeft.z_max = bboxRight.z_min = splitPoint;
        }

        //log of creating tree
        // std::cout << "=================="  << std::endl;
        // std::cout << leftList.size() << std::endl;
        // std::cout << bboxLeft.x_min << "|" << bboxLeft.x_max << std::endl;
        // std::cout << bboxLeft.y_min << "|" << bboxLeft.y_max << std::endl;
        // std::cout << bboxLeft.z_min << "|" << bboxLeft.z_max << std::endl;
        // std::cout << rightList.size() << std::endl;
        // std::cout << bboxRight.x_min << "|" << bboxRight.x_max << std::endl;
        // std::cout << bboxRight.y_min << "|" << bboxRight.y_max << std::endl;
        // std::cout << bboxRight.z_min << "|" << bboxRight.z_max << std::endl;

        meshTree->left = build(leftList, bboxLeft);
        meshTree->right = build(rightList, bboxRight);

        return meshTree;
    }
}
CompFab::BoundingBox::BoundingBoxStruct(){
    x_min = x_max = y_min = y_max = z_min = z_max = 0;
}


void CompFab::BoundingBox::update(double x_min_new, double y_min_new, double z_min_new, 
                                  double x_max_new, double y_max_new, double z_max_new)
{
    if (x_min_new < x_min) {
        x_min = x_min_new;
        }
    if (y_min_new < y_min) {
        y_min = y_min_new;
        }
    if (z_min_new < z_min) {
        z_min = z_min_new;
        }
    if (x_max_new > x_max) {
        x_max = x_max_new;
        }
    if (y_max_new > y_max) {
        y_max = y_max_new;
        }
    if (z_max_new > z_max) {
        z_max = z_max_new;
    }
}

int CompFab::BoundingBox::maxAxis() 
{
    double x_diff = x_max - x_min;
    double y_diff = y_max - y_min;
    double z_diff = z_max - z_min;
    if (x_diff >= y_diff && x_diff >= z_diff) {
        return 1;
    }
    else if(y_diff >= x_diff && y_diff >= z_diff) {
        return 2;
    }
    return 3;
}




