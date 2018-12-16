/************************************************************************************************
 * L3 EEA - Projet Techniques Scientifiques
 * Theme : Algorithme Genetique
 * Auteurs:
 *      - AIT MOUHOUB Farid
 *      - AOUCI Sofiane
 *      - JAQUET E.Prince
 *      - LAVAL Hugo
 ************************************************************************************************/

#ifndef NONGENFUN_H
#define NONGENFUN_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define octave "C:\\Octave\\Octave-4.4.1\\bin"


void plot(int soft, double Y0, double lambda_L, double lambda_0, double Ymax);
double electronic_density(double lambda_L);

double sqr(double x);


#endif // NONGENFUN_H
