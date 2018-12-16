/************************************************************************************************
 * L3 EEA - Projet Techniques Scientifiques
 * Theme : Algorithme Genetique
 * Auteurs:
 *      - AIT MOUHOUB Farid
 *      - AOUCI Sofiane
 *      - JAQUET E.Prince
 *      - LAVAL Hugo
 ************************************************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "nongenfun.h"
#include "genfun.h"

int main()
{
    srand((unsigned int) time(0)); //initialisation of the rand() function
    time_t duration = time(0);  //time at which the program has been executed to calculate how much time it takes to find the solution
    double data[data_l][3], eval[inc] = {0};
    double Y0 = 0, lambda_L = 0, lambda_0 = 0, Ymax = 0; //Our constants
    int i = 0, i_sol = 0;
    uint64_t *chromo = calloc(inc,sizeof(uint64_t));

    FILE *profile = fopen("profil.txt","r");

    if(profile == NULL) exit(-1);

    for(i = 0; i < data_l; i++)
    {
        fscanf(profile,"%lf %lf %lf",&data[i][0],&data[i][1],&data[i][2]); //retreiving the data from profile.txt
    }
    fclose(profile);

    init(chromo); // initialisation of the chromosomes randomly (it's the first generation)
    i = 0;
    do
    {
        evaluate(chromo,data,eval); //We evaluate the chromosomes
        survival(&chromo,eval); //We build the new generation from the old one
        crossover(chromo);  //make the crossover operation
        mutation(chromo);   //make the mutation
        printf("%lf     %d\n",best_chromo(eval, &i_sol),i++); //print the error of the best chromosome of the new generation and the current itteration
        //if the number of itterations reaches 5000, the program is stuck and is very slowly converging to the solution, so we
        //reinitialise the chromosomes using init()
        if(i > max_gen)
        {
            init(chromo); 
            i = 0;
        }
    }
    while(best_chromo(eval, &i_sol) > ferror); //do all this until we go under an error of 0.204

    printf("Algo duration : %d s\n",(unsigned int)(time(0) - duration)); //display the execution time

    decode(chromo[i_sol],&Y0,&lambda_L,&lambda_0,&Ymax);
    printf("Ymax = %lf, lambda_0 = %lf, lambda_L = %lf, Y0 = %lf\n",Ymax,lambda_0,lambda_L,Y0); //display our constants
    
    //plot the curves using matlab or octave
    printf("Which software do you want to use ? 0 - Matlab, 1 - Octave\n");
    scanf("%d",&i);
    plot(i, Y0, lambda_L, lambda_0, Ymax);

    return 0;
}
