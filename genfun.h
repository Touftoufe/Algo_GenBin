/************************************************************************************************
 * L3 EEA - Projet Techniques Scientifiques
 * Theme : Algorithme Genetique
 * Auteurs:
 *      - AIT MOUHOUB Farid
 *      - AOUCI Sofiane
 *      - JAQUET E.Prince
 *      - LAVAL Hugo
 ************************************************************************************************/


/* Constants intervals and conversion rules
 * lambda_0 E [6555, 6570] coded in 19 bit integer ~ [0, 524287] <=> [6540, 6592.4287]
 * lambda_L E [7, 15.191] coded in 13 bit integer ~ [0, 8191] <=> [7, 15.191]
 * Y_0      E [0.2, 0.25] coded in 16 bit integer ~ [0, 65535] <=> [0.2, 0.265535]
 * Y_max    E [0.72, 0.78] coded in 16 bit integer ~ [0, 65535] <=> [0.72, 0.785535]
 */


#ifndef GENFUN_H
#define GENFUN_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "nongenfun.h"

#define pi 3.14159265358979323846

//lower bounds of our constants sets
#define inf_lam_0 6540
#define inf_lam_L 7
#define inf_Y0 0.2
#define inf_Ymax 0.72

//length of the profil.txt file
#define data_l 1024

#define inc 100         //number of chromosomes
#define X 0.15          //Crossover rate
#define Mu 0.015        //Mutation rate
#define c_length 64     //Cromosome length (bits)
#define ferror 0.204    //max overall error
#define max_gen 5000    //maximun number of generations (itterations)


void init(uint64_t* chromo);
int decode(uint64_t chromo, double* Y0, double* lambda_L, double* lambda_0, double* Ymax);
uint64_t code(double Y0, double lambda_L, double lambda_0, double Ymax);
double fitting(uint64_t c, double lambda);
void evaluate(uint64_t *chromo, double data[][3], double *eval);
void survival(uint64_t **chromo, double* eval);
void crossover(uint64_t *chromo);
void mutation(uint64_t *chromo);
double best_chromo(double *eval, int *i_sol);


#endif // GENFUN_H
