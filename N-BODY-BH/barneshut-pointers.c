#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>


const int N = 1000;
const double theta = 0.5;
//const double M_PI = 3.14159216;


const double G = 4 * 3.14159216 * 3.14159216; 
const double dt = 0.001; 
const int nSteps = 1200;
const double epsilon = 1e-7;


typedef struct Vec3 {
    double x, y, z;
} Vec3;

typedef struct {
    double  *r;
    double  *v;
    double  *m;
} body;

typedef struct Node {
    double *centerOfMass;
    double *totalMass;
    double *bbox[2];
    struct Node *children[8];
    body* bd;
} Node;


Node* Init_Node(double min, double max);
Node* createChildNodeForOctant(Node* node, int octantIndex);
int OctantIndex(Node* node, body* body);
void Init_bd(body* bd, int N, double min, double max);
void updateMassCenterOfMass(Node* nd, body* bd);
void subdivideNode(Node* nd);
void insertbody(Node* node, body* bd);
void CenterOfMass(Node* node);
void Force(Node* node, body* bd, Vec3* f);
void Euler(body* bd, Vec3 force, double dt);
void freeNode(Node* node);
double randomDouble(double min, double max);
void writePositionsToCSV(body* bd, int numbd, const char* filename);
void Init_bodys2(body* bd, int N, Vec3 center);



int main() {

    body bd[N];
    double *force = (double *)malloc(3 * N * sizeof(double));
    double rootMin = {-70, -70, -70};
    double rootMax = {70, 70, 70};
    double centerg = {0,0,0};

    // Init quantities
    //Init_bodys2(bd, N, centerg);
    Init_bd(bd, N, rootMin, rootMax);
    Node* rootNode = Init_Node(rootMin, rootMax);

    // Initial insertion of bodys into the tree
    for (size_t ii = 0; ii < N; ii++) {
        insertbody(rootNode, &bd[ii]);
    }

    char filename[nSteps];

    for (int step = 0; step < nSteps; step++) {
        // Write positions
        if (step % 1 == 0)
        {
            sprintf(filename, "positions_%d.csv", step);
            writePositionsToCSV(bd, N, filename);
        }

        // clean forces
        memset(force, 0, sizeof(force));

        // calculate force of each body
        for (size_t ii = 0; ii < N; ii++) {

            Vec3 f = {0, 0, 0}; 
            Force(rootNode, &bd[ii], &f);
            force[ii] = f; 
        }

        // Update vel and acc
        for (size_t ii = 0; ii < N; ii++) {
            Euler(&bd[ii], force[ii], dt);
        }

        // Clean last node
        freeNode(rootNode);
        rootNode = Init_Node(rootMin, rootMax);

        for (size_t ii = 0; ii < N; ii++) {
            insertbody(rootNode, &bd[ii]); 
        }

        CenterOfMass(rootNode);


    }
    


    // CLEAN
    freeNode(rootNode);

    return 0;
}


Node* Init_Node(double *min, double *max){

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

    node->bd = NULL;

    return node;
}


Node* createChildNodeForOctant(Node* bd, int octantIndex) {
    Vec3 center = {
        (bd->bbox[0].x + bd->bbox[1].x) / 2,
        (bd->bbox[0].y + bd->bbox[1].y) / 2,
        (bd->bbox[0].z + bd->bbox[1].z) / 2
    };

    Vec3 min = bd->bbox[0];
    Vec3 max = bd->bbox[1];

    // Octant index
    if (octantIndex & 4) { 
        min.x = center.x;
    } else {
        max.x = center.x;
    }

    if (octantIndex & 2) { 
        min.y = center.y;
    } else {
        max.y = center.y;
    }

    if (octantIndex & 1) { 
        min.z = center.z;
    } else {
        max.z = center.z;
    } 


    return Init_Node(min, max);
}


double randomDouble(double min, double max) {
    return min + (rand() / (RAND_MAX / (max - min)));
}


void Init_bodys(body *bd, int N, Vec3 min, Vec3 max) {
    // seed
    srand(time(NULL));

    bd->r = (double *)malloc(3 * N * sizeof(double));
    bd->v = (double *)malloc(3 * N * sizeof(double));
    bd->m = (double *)malloc(N * sizeof(double));

    for (size_t ii = 0; ii < N; ii++) {
        // random pos in domain
        bd->r[3*ii]     = randomDouble(min.x + 20, max.x - 20);  // x component
        bd->r[3*ii + 1] = randomDouble(min.y + 20, max.y - 20);  // y component
        bd->r[3*ii + 2] = randomDouble(min.z + 20, max.z - 20);  // z component

        bd->v[3*ii]     = 0;
        bd->v[3*ii + 1] = 0;
        bd->v[3*ii + 2] = 0;

        // mass
        bd->m[ii] = 1.0; 
    }
}


void Init_bodys2(body* bd, int N, Vec3 center) {
    srand(time(NULL));
    double galaxy_radius = 50;
    double disk_thickness = 1.5;
    double galaxy_mass = N;

    bd->r = (double *)malloc(3 * N * sizeof(double));
    bd->v = (double *)malloc(3 * N * sizeof(double));
    bd->v = (double *)malloc(N * sizeof(double));



    for (size_t ii = 0; ii < N; ii++) {
        // Distribución radial desde el centro de la galaxia
        double radius = ((double)rand() / RAND_MAX) * galaxy_radius;
        double theta = ((double)rand() / RAND_MAX) * 2 * 3.14159216; // Ángulo aleatorio en radianes

        // Altura aleatoria para simular el espesor del disco de la galaxia
        double height = ((double)rand() / RAND_MAX) * disk_thickness - disk_thickness / 2;

        // Conversión de coordenadas polares a cartesianas
        bd->r[3*ii ]    = center.x + radius * cos(theta);
        bd->r[3*ii + 1] = center.y + radius * sin(theta);
        bd->r[3*ii + 2] = center.x + height;


        // Velocidad inicial: orbital simple para simular rotación galáctica
        double orbital_velocity = 1.5*sqrt((G * galaxy_mass) / radius); // Fórmula simplificada
        bd->v[3*ii]     = -orbital_velocity * sin(theta);
        bd->v[3*ii + 1] = orbital_velocity * cos(theta);
        bd->v[3*ii + 2] = 0.0; // Sin movimiento inicial en la dirección Z

        // Masa
        bd->m[3*ii] = 1.0; // Considera ajustar según la distribución de masa deseada
    }
}

  
int OctantIndex(Node* node, body* bd) {
    // Calcula el centro del nodo
    Vec3 center;
    center.x = (node->bbox[0].x + node->bbox[1].x) / 2;
    center.y = (node->bbox[0].y + node->bbox[1].y) / 2;
    center.z = (node->bbox[0].z + node->bbox[1].z) / 2;

    int index = 0;

    //bd->r = (double *)malloc(3 * sizeof(double));
    //bd->v = (double *)malloc(3 * sizeof(double));

    // X-axis
    if (bd->r[0] >= center.x) {
        index += 4;
    }
    // Y-xis
    if (bd->r[1] >= center.y) {
        index += 2;
    }
    // Z-axis
    if (bd->r[2] >= center.z) {
        index += 1;
    }

    return index;
}



void updateMassCenterOfMass(Node* node, body* bd){
    if (node->totalMass == 0){
        //Empty node
        node->centerOfMass = bd->r;
    }
    else { //Average CM
        node->centerOfMass.x = (node->centerOfMass.x * node->totalMass + bd->r.x * body->mass) / (node->totalMass + body->mass);
        node->centerOfMass.y = (node->centerOfMass.y * node->totalMass + body->position.y * body->mass) / (node->totalMass + body->mass);
        node->centerOfMass.z = (node->centerOfMass.z * node->totalMass + body->position.z * body->mass) / (node->totalMass + body->mass);
    }
    node->totalMass += bd->ma;
}


void insertbody(Node* node, body* bd) {
    // Empty Node
    if (node->bd == NULL && node->totalMass == 0) {
        node->bd = bd;
        updateMassCenterOfMass(node, bd);
        return;
    }

    //If the node is a leaf but already contains a body, subdivide the node.
    if (node->bd != NULL && node->totalMass != 0) {
        subdivideNode(node);  
    }

    // Insert Child
    int newbodyIndex = OctantIndex(node, bd);
    if (node->children[newbodyIndex] == NULL) {
        node->children[newbodyIndex] = createChildNodeForOctant(node, newbodyIndex);
    }

    insertbody(node->children[newbodyIndex], bd);

    // Update the parent node Center of Mass
    updateMassCenterOfMass(node, bd);
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

    //Count bodys
    if (node->bd !=NULL){
        node->centerOfMass.x = node->bd->position.x * node->bd->mass;
        node->centerOfMass.y = node->bd->position.y * node->bd->mass;
        node->centerOfMass.z = node->bd->position.z * node->bd->mass;
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


void Force(Node* node, body* bd, Vec3* f ){
    // Emty node
    if (node->totalMass == 0) return;
    
    Vec3 dir;
    double r_ij = distance(bd->position, node->centerOfMass);

    if (r_ij < epsilon) return;

    // Avoid divergence
    if (r_ij == epsilon) return;

    // size of node
    double nodeSize = node->bbox[1].x - node->bbox[0].x;
    double h = nodeSize / r_ij;

    // If the node is far enough away or is a leaf, treat the node as a point mass
    if (h < theta || node->children[0] == NULL ){

        dir.x = node->centerOfMass.x - bd->position.x;
        dir.y = node->centerOfMass.y - bd->position.y;
        dir.z = node->centerOfMass.z - bd->position.z;

        double F = -G * bd->mass * node->totalMass / (r_ij*r_ij*r_ij);
        f->x += F*dir.x;
        f->y += F*dir.y;
        f->z += F*dir.z;
    }
    else{
        // Otherwise calculate force recursively on subnodes
        for (size_t ii = 0; ii < 8; ii++){
            if (node->children[ii] != NULL){
                Force(node->children[ii], bd, f);
            }    
        }   
    }

}


void Euler(body* bd, Vec3 force, double dt) {
    // Calculate Aceleration
    Vec3 acceleration = {force.x / bd->mass, 
                         force.y / bd->mass, 
                         force.z / bd->mass};

    // Update Velocity
    bd->velocity.x += 0.5*acceleration.x * dt;
    bd->velocity.y += 0.5*acceleration.y * dt;
    bd->velocity.z += 0.5*acceleration.z * dt;

    // Update Positions
    bd->position.x += 0.5*bd->velocity.x * dt;
    bd->position.y += 0.5*bd->velocity.y * dt;
    bd->position.z += 0.5*bd->velocity.z * dt;
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


void writePositionsToCSV(body* bd, int numbd, const char* filename) {
    FILE* file = fopen(filename, "w"); 
    if (file == NULL) {
        printf("Error al abrir el archivo.\n");
        return;
    }

    fprintf(file, "X,Y,Z\n");

    for (int i = 0; i < numbd; i++) {
        fprintf(file, "%f,%f,%f\n", bd[i].position.x, bd[i].position.y, bd[i].position.z);
    }

    fclose(file);
}

