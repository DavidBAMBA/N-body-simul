#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>


const int N = 1000;
const double theta = 0.5;

const double G = 4 * (3.14159216)*(3.14159216); 
const double dt = 0.01; 
const int nSteps = 1200;
const double epsilon = 1e-7;


typedef struct Vec3 {
    double x, y, z;
} Vec3;

typedef struct Particle{
    Vec3 position;
    Vec3 velocity;
    double mass;
} Particle;

typedef struct Node {
    Vec3 centerOfMass;
    double totalMass;
    Vec3 bbox[2];
    struct Node* children[8];
    Particle* particle;
} Node;


Node* Init_Node(Vec3 min, Vec3 max);
Node* createChildNodeForOctant(Node* parent, int octantIndex);
int OctantIndex(Node* node, Particle* particle);
void Init_Particles(Particle* particles, int N, Vec3 min, Vec3 max);
void updateMassCenterOfMass(Node* node, Particle* particle);
void subdivideNode(Node* node);
void insertParticle(Node* node, Particle* particle);
void CenterOfMass(Node* node);
void Force(Node* node, Particle* particle, Vec3* f);
void Euler(Particle* particle, Vec3 force, double dt);
void freeNode(Node* node);
double randomDouble(double min, double max);
void writePositionsToCSV(Particle* particles, int numParticles, const char* filename);


int main() {

    Particle particles[N];
    Vec3 force[N];
    Vec3 rootMin = {-70, -70, -70};
    Vec3 rootMax = {70, 70, 70};

    // Init quantities
    Init_Particles(particles, N, rootMin, rootMax);
    Node* rootNode = Init_Node(rootMin, rootMax);

    // Initial insertion of particles into the tree
    for (size_t ii = 0; ii < N; ii++) {
        insertParticle(rootNode, &particles[ii]);
    }

    char filename[nSteps];

    for (int step = 0; step < nSteps; step++) {
        // Write positions
        sprintf(filename, "positions_%d.csv", step);
        writePositionsToCSV(particles, N, filename);

        // clean forces
        memset(force, 0, sizeof(force));

        // calculate force of each particle
        for (size_t ii = 0; ii < N; ii++) {

            Vec3 f = {0, 0, 0}; 
            Force(rootNode, &particles[ii], &f);
            force[ii] = f; 
        }

        // Update vel and acc
        for (size_t ii = 0; ii < N; ii++) {
            Euler(&particles[ii], force[ii], dt);
        }

        // Clean last node
        freeNode(rootNode);
        rootNode = Init_Node(rootMin, rootMax);

        for (size_t ii = 0; ii < N; ii++) {
            insertParticle(rootNode, &particles[ii]); 
        }

        CenterOfMass(rootNode);


    }


    // CLEAN
    freeNode(rootNode);

    return 0;
}


Node* Init_Node(Vec3 min, Vec3 max){

    Node* node = (Node*)malloc(sizeof(Node));
    if (node == NULL)
    {
        fprintf(stderr, "Error al asignar memoria para el nodo\n");
        exit(EXIT_FAILURE);
    }

    //Init node parameters
    node->bbox[0] = min; 
    node->bbox[1] = max;
    node->centerOfMass = (Vec3){0,0,0};
    node->totalMass = 0;
    
    for(int ii = 0; ii< 8 ; ++ii){
        node->children[ii]= NULL;
    }

    node->particle = NULL;

    return node;
}


Node* createChildNodeForOctant(Node* parent, int octantIndex) {
    Vec3 center = {
        (parent->bbox[0].x + parent->bbox[1].x) / 2,
        (parent->bbox[0].y + parent->bbox[1].y) / 2,
        (parent->bbox[0].z + parent->bbox[1].z) / 2
    };

    Vec3 min = parent->bbox[0];
    Vec3 max = parent->bbox[1];

    // octant limits
    switch (octantIndex) {
        case 0: // Octant -x, -y, -z
            max.x = center.x;
            max.y = center.y;
            max.z = center.z;
            break;
        case 1: // Octant -x, -y, +z
            max.x = center.x;
            max.y = center.y;
            min.z = center.z;
            break;
        case 2: // Octant -x, +y, -z
            max.x = center.x;
            min.y = center.y;
            max.z = center.z;
            break;
        case 3: // Octant -x, +y, +z
            max.x = center.x;
            min.y = center.y;
            min.z = center.z;
            break;
        case 4: // Octant +x, -y, -z
            min.x = center.x;
            max.y = center.y;
            max.z = center.z;
            break;
        case 5: // Octant +x, -y, +z
            min.x = center.x;
            max.y = center.y;
            min.z = center.z;
            break;
        case 6: // Octant +x, +y, -z
            min.x = center.x;
            min.y = center.y;
            max.z = center.z;
            break;
        case 7: // Octant +x, +y, +z
            min.x = center.x;
            min.y = center.y;
            min.z = center.z;
            break;
    }

    return Init_Node(min, max);
}


double randomDouble(double min, double max) {
    return min + (rand() / (RAND_MAX / (max - min)));
}


void Init_Particles(Particle* particles, int N, Vec3 min, Vec3 max) {
    // seed
    srand(time(NULL));

    for (size_t ii = 0; ii < N; ii++) {
        // random pos in domain
        particles[ii].position.x = randomDouble(min.x+20, max.x-20);
        particles[ii].position.y = randomDouble(min.y+20, max.y-20);
        particles[ii].position.z = randomDouble(min.z+20, max.z-20);

        particles[ii].velocity.x = 0.01; 
        particles[ii].velocity.y = 0.01; 
        particles[ii].velocity.z = 0.01; 

        // mass
        particles[ii].mass = 1.0; 
    }
}


int OctantIndex(Node* node, Particle* particle) {
    Vec3 center = {
        (node->bbox[0].x + node->bbox[1].x) / 2,
        (node->bbox[0].y + node->bbox[1].y) / 2,
        (node->bbox[0].z + node->bbox[1].z) / 2
    };
    
    // Determine Octant
    if (particle->position.x < center.x) {
        if (particle->position.y < center.y) {
            if (particle->position.z < center.z) {
                return 0; // Octant -x, -y, -z
            } else {
                return 1; // Octant -x, -y, +z
            }
        } else {
            if (particle->position.z < center.z) {
                return 2; // Octant -x, +y, -z
            } else {
                return 3; // Octant -x, +y, +z
            }
        }
    } else {
        if (particle->position.y < center.y) {
            if (particle->position.z < center.z) {
                return 4; // Octant +x, -y, -z
            } else {
                return 5; // Octant +x, -y, +z
            }
        } else {
            if (particle->position.z < center.z) {
                return 6; // Octant +x, +y, -z
            } else {
                return 7; // Octant +x, +y, +z
            }
        }
    }
}


void updateMassCenterOfMass(Node* node, Particle* particle){
    if (node->totalMass == 0){
        //Empty node
        node->centerOfMass = particle->position;
    }
    else { //Average CM
        node->centerOfMass.x = (node->centerOfMass.x * node->totalMass + particle->position.x * particle->mass) / (node->totalMass + particle->mass);
        node->centerOfMass.y = (node->centerOfMass.y * node->totalMass + particle->position.y * particle->mass) / (node->totalMass + particle->mass);
        node->centerOfMass.z = (node->centerOfMass.z * node->totalMass + particle->position.z * particle->mass) / (node->totalMass + particle->mass);
    }
    node->totalMass += particle->mass;
}


void insertParticle(Node* node, Particle* particle) {
    // Empty Node
    if (node->particle == NULL && node->totalMass == 0) {
        node->particle = particle;
        updateMassCenterOfMass(node, particle);
        return;
    }

    //If the node is a leaf but already contains a particle, subdivide the node.
    if (node->particle != NULL && node->totalMass != 0) {
        subdivideNode(node);  
    }

    // Insert Child
    int newParticleIndex = OctantIndex(node, particle);
    if (node->children[newParticleIndex] == NULL) {
        node->children[newParticleIndex] = createChildNodeForOctant(node, newParticleIndex);
    }

    insertParticle(node->children[newParticleIndex], particle);

    // Update the parent node Center of Mass
    updateMassCenterOfMass(node, particle);
}


void subdivideNode(Node* node) 
{
    Vec3 center = {
        (node->bbox[0].x + node->bbox[1].x) / 2,
        (node->bbox[0].y + node->bbox[1].y) / 2,
        (node->bbox[0].z + node->bbox[1].z) / 2
    };

    // Create child nodes
    for (int ii = 0; ii < 8; ii++) {
        if (node->children[ii] == NULL) {
            node->children[ii] = createChildNodeForOctant(node, ii);
        }
    }
    
}


double distance(Vec3 ri, Vec3 rj)
{
    return sqrt( (rj.x - ri.x)*(rj.x - ri.x) + (rj.y-ri.y)*(rj.y-ri.y) + (rj.z-ri.z)*(rj.z-ri.z) );
}


void CenterOfMass(Node* node)
{
    node->centerOfMass = (Vec3){0,0,0};
    node->totalMass = 0;

    int count = 0;

    //Count particles
    if (node->particle !=NULL){
        node->centerOfMass.x = node->particle->position.x * node->particle->mass;
        node->centerOfMass.y = node->particle->position.y * node->particle->mass;
        node->centerOfMass.z = node->particle->position.z * node->particle->mass;
    }

    // Run over child nodes
    for (size_t ii = 0; ii < count; ii++){
        if (node->children[ii] != NULL){

            // Recursive call
            CenterOfMass(node->children[ii]);

            // Sum of Mass
            node->totalMass += node->children[ii]->totalMass;

            if (node->children[ii]->totalMass > 0){
                node->centerOfMass.x += node->children[ii]->centerOfMass.x * node->children[ii]->totalMass;
                node->centerOfMass.y += node->children[ii]->centerOfMass.y * node->children[ii]->totalMass;
                node->centerOfMass.z += node->children[ii]->centerOfMass.z * node->children[ii]->totalMass;
                count++;
            } 
        }
    }

    // if there are mass, adjust the CM
    if (node->totalMass > 0 && count > 0 ){
        node->centerOfMass.x /=node->totalMass;
        node->centerOfMass.y /=node->totalMass;
        node->centerOfMass.z /=node->totalMass;
    }
     
}


void Force(Node* node, Particle* particle, Vec3* f ){
    // Emty node
    if (node->totalMass == 0) return;
    
    Vec3 dir;
    double r_ij = distance(particle->position, node->centerOfMass);

    if (r_ij < epsilon) return;

    // Avoid divergence
    if (r_ij == epsilon) return;

    // size of node
    double nodeSize = node->bbox[1].x - node->bbox[0].x;
    double h = nodeSize / r_ij;

    // If the node is far enough away or is a leaf, treat the node as a point mass
    if (h < theta || node->children[0] == NULL ){

        dir.x = node->centerOfMass.x - particle->position.x;
        dir.y = node->centerOfMass.y - particle->position.y;
        dir.z = node->centerOfMass.z - particle->position.z;

        double F = G * particle->mass * node->totalMass / (r_ij*r_ij*r_ij);
        f->x += F*dir.x;
        f->y += F*dir.y;
        f->z += F*dir.z;
    }
    else{
        // Otherwise calculate force recursively on subnodes
        for (size_t ii = 0; ii < 8; ii++){
            if (node->children[ii] != NULL){
                Force(node->children[ii], particle, f);
            }    
        }   
    }

}


void Euler(Particle* particle, Vec3 force, double dt) {
    // Calculate Aceleration
    Vec3 acceleration = {force.x / particle->mass, 
                         force.y / particle->mass, 
                         force.z / particle->mass};

    // Update Velocity
    particle->velocity.x += 0.5*acceleration.x * dt;
    particle->velocity.y += 0.5*acceleration.y * dt;
    particle->velocity.z += 0.5*acceleration.z * dt;

    // Update Positions
    particle->position.x += 0.5*particle->velocity.x * dt;
    particle->position.y += 0.5*particle->velocity.y * dt;
    particle->position.z += 0.5*particle->velocity.z * dt;
}


void freeNode(Node* node) {
    // Empty node
    if (node == NULL) {
        return;
    }
    
    // Run over nodes
    for (size_t ii = 0; ii < 8; ii++) {
        if (node->children[ii] != NULL) {
            freeNode(node->children[ii]);
            node->children[ii] = NULL; 
        }
    }

    free(node);
}


void writePositionsToCSV(Particle* particles, int numParticles, const char* filename) {
    FILE* file = fopen(filename, "w"); 
    if (file == NULL) {
        printf("Error al abrir el archivo.\n");
        return;
    }

    fprintf(file, "X,Y,Z\n");

    for (int i = 0; i < numParticles; i++) {
        fprintf(file, "%f,%f,%f\n", particles[i].position.x, particles[i].position.y, particles[i].position.z);
    }

    fclose(file);
}



