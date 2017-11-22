//Computational Fabrication Assignment #1
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "../include/CompFab.h"
#include "../include/Mesh.h"

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */
    CompFab::Vec3 normal = (triangle.m_v2 - triangle.m_v1) % (triangle.m_v3 - triangle.m_v1);
    normal.normalize();
    double D = - (normal * triangle.m_v1);
    double isIntersect = normal * ray.m_direction;
    if (isIntersect < EPSILON && isIntersect > -EPSILON) {
        return 0;
    }
    double t = -((normal * ray.m_origin) + D) / isIntersect;
    if (t < 0) {
        return 0;
    }
    CompFab::Vec3 P = ray.m_origin + CompFab::Vec3 (ray.m_direction[0] * t, 
                                                    ray.m_direction[1] * t, 
                                                    ray.m_direction[2] * t);
    CompFab::Vec3 R = P - triangle.m_v1;
    CompFab::Vec3 Q1 = triangle.m_v2 - triangle.m_v1;
    CompFab::Vec3 Q2 = triangle.m_v3 - triangle.m_v1;
    double Q1Q1 = Q1 * Q1;
    double Q1Q2 = Q1 * Q2;
    double Q2Q2 = Q2 * Q2;
    double RQ1 = R * Q1;
    double RQ2 = R * Q2;
    double det = (Q1Q1 * Q2Q2 - Q1Q2 * Q1Q2);
    if (det < EPSILON) {
        return 0;
    }
    double w1 = (Q2Q2 * RQ1 - Q1Q2 * RQ2) / det;
    double w2 = (Q1Q1 * RQ2 - Q1Q2 * RQ1) / det;
    if (1 - w1 - w2 > 0 && w1 > 0 && w2 > 0) {
        return 1;
    }

    return 0;

}

bool rayBoxIntersection(CompFab::Ray &ray, CompFab::BoundingBox &bbox) 
{
    CompFab::Vec3 min(bbox.x_min, bbox.y_min, bbox.z_min);
    CompFab::Vec3 max(bbox.x_max, bbox.y_max, bbox.z_max);

    double txmin = (min.m_x - ray.m_origin[0]) / ray.m_direction[0];
    double txmax = (max.m_x - ray.m_origin[0]) / ray.m_direction[0];
    double tymin = (min.m_y - ray.m_origin[1]) / ray.m_direction[1];
    double tymax = (max.m_y - ray.m_origin[1]) / ray.m_direction[1];
    double tzmin = (min.m_z - ray.m_origin[2]) / ray.m_direction[2];
    double tzmax = (max.m_z - ray.m_origin[2]) / ray.m_direction[2];

    if (txmin > txmax) std::swap(txmin, txmax);
    if (tymin > tymax) std::swap(tymin, tymax);
    if (tzmin > tzmax) std::swap(tzmin, tzmax);

    if (txmin > tymax || tymin > txmax) {
        return false;
    }
    txmin = std::max(txmin, tymin);
    tymax = std::min(txmax, tymax);

    if(txmin > tzmax || tzmin > txmax) {
        return false;
    }

    return true;
}

//Tree version of counting intersections
int numSurfaceIntersectionsTree(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir, CompFab::MeshTree & meshTree) {
    dir.normalize();
    CompFab::Ray ray(voxelPos, dir);
    int count = 0;
    if (rayBoxIntersection(ray, meshTree.bbox)) {
        if (meshTree.isLeaf) {
            for (unsigned int i=0; i < meshTree.meshs.size(); i++) {
                count += rayTriangleIntersection(ray, meshTree.meshs[i]);
            }
        }
        else {
            count = numSurfaceIntersectionsTree(voxelPos, dir, *meshTree.left) +
                    numSurfaceIntersectionsTree(voxelPos, dir, *meshTree.right);
        }
    }
    return count;
}

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */
    dir.normalize();
    CompFab::Ray ray(voxelPos, dir);
    for (unsigned int i=0; i < g_triangleList.size(); i++) {
        numHits += rayTriangleIntersection(ray, g_triangleList[i]);
    }
    return numHits;
}

bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return true;
   
}

void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}


int main(int argc, char **argv)
{
    clock_t t1,t2;
    t1=clock();
    unsigned int dim = 64; //dimension of voxel grid (e.g. 32x32x32)
    std::cout<<"=========="<<std::endl;
    //Load OBJ
    if(argc < 3)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        return 0;
    }
    
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], dim);
    

    
    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction;
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */
    CompFab::MeshTree* meshTree = CompFab::MeshTree::MeshTree().build(g_triangleList);
    std::cout << "Tree Built" << std::endl;

    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    bool useTree = false;
    bool multipleRay = false;
    CompFab::Vec3 ll = g_voxelGrid->m_lowerLeft;
    double spacing = g_voxelGrid->m_spacing;
    for (int ii = 0; ii < nx; ii++) {
        voxelPos[0] = ll[0]+ ii * spacing;
        for (int jj = 0; jj < ny; jj++) {
            voxelPos[1] = ll[1] + jj * spacing;
            for (int kk = 0; kk < nz; kk++) {
                voxelPos[2] = ll[2] + kk * spacing;
                if (!multipleRay) {
                    direction = CompFab::Vec3((double)(rand() % 100 - 50),(double)(rand() % 100 - 50),(double)(rand() % 100 - 50));
                    if (!useTree && numSurfaceIntersections(voxelPos, direction) % 2 == 1) {
                        g_voxelGrid->isInside(ii,jj,kk) = true;
                    }
                    if (useTree && numSurfaceIntersectionsTree(voxelPos, direction, *meshTree) % 2 == 1) {
                        g_voxelGrid->isInside(ii,jj,kk) = true;
                    }
                }
                if (multipleRay) {
                    int count = 0;
                    for (int n = 0; n < 9; n++) {
                        direction = CompFab::Vec3((double)(rand() % 100 - 50),(double)(rand() % 100 - 50),(double)(rand() % 100 - 50));
                        if (!useTree && numSurfaceIntersections(voxelPos, direction) % 2 == 1) {
                            count ++;
                        }
                        if (useTree && numSurfaceIntersectionsTree(voxelPos, direction, *meshTree) % 2 == 1) {
                            count ++;
                        }
                    }
                    if (count >= 5) {
                        g_voxelGrid->isInside(ii,jj,kk) = true;
                    }
                }
            }
        }
    }


    
    //Write out voxel data as obj
    saveVoxelsToObj(argv[2]);
    t2=clock();
    float diff ((float)t2-(float)t1);
    if (useTree) {
        std::cout<<"Using Tree: True"<<std::endl;
    }
    else {
        std::cout<<"Using Tree: False"<<std::endl;
    }
    if (multipleRay) {
        std::cout<<"Multiple Ray: True"<<std::endl;
    }
    else {
        std::cout<<"Multiple Ray: False"<<std::endl;
    }
    std::cout<<"Dimension: "<<dim<<std::endl;
    std::cout<<"Runtime in seconds: "<<diff / CLOCKS_PER_SEC<<std::endl;
    std::cout<<"=========="<<std::endl;
    
    delete g_voxelGrid;
}