/**
 * Solar System with test particles
 *
 * This example integrates all planets of the Solar
 * System and 10000 test particles. The initial data comes 
 * from the NASA HORIZONS system and was saved to
 * a binary file beforehand. The integrator used is WHFast
 * with a 4 day timestep. Note that close encounters are
 * not resolved. The OpenMP speedup you get depends on the 
 * compiler and CPU that you are using. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

#define FILENAME_SIZE 64
#define MAX_LINE 256
#define AST_NUM 10




void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation_from_binary("ss-2020-03-03.bin");
    // Setup constants
    r->dt           = 0.5/365.25*2.*M_PI;        // 0.5 days
    r->heartbeat    = heartbeat;
    r->N_active     = r->N;
    r->ri_ias15.min_dt       = 0.01;

    // Read Synthetic Population
    FILE *fp;
    char buffer[MAX_LINE];

 

    fp = fopen("asteroids.cat", "r");
    if (fp == NULL){
        perror("Error");
        return(-1);
    } 
    for (int i=0; i<AST_NUM; i++){
    
        fgets(buffer, MAX_LINE, fp);
        double a = atof(buffer);
        fgets(buffer, MAX_LINE, fp);
        double e = atof(buffer);
        fgets(buffer, MAX_LINE, fp);
        double i = atof(buffer);
        fgets(buffer, MAX_LINE, fp);
        double node = atof(buffer);
        fgets(buffer, MAX_LINE, fp);
        double omega = atof(buffer);
        fgets(buffer, MAX_LINE, fp);
        double M = atof(buffer);
        double f = reb_tools_M_to_f(e, M);
        struct reb_particle p = reb_tools_orbit_to_particle(1.,r->particles[0],0.,a,e,i,node,omega,f);
        reb_add(r, p);
    }
    fclose(fp);

    reb_move_to_com(r);
    reb_integrate(r, 365.25*2.*M_PI);
}

void heartbeat(struct reb_simulation* r){
    printf("%f\n",r->t);
    if (reb_output_check(r, 1.00)){
        printf("yee XD");
        char* file = "testoutput.txt";
        reb_output_ascii(r, file);
        
    }
}

