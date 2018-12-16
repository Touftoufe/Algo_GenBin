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
    srand((unsigned int) time(0));
    time_t duration = time(0);
    double data[data_l][3], eval[inc] = {0};;
    double Y0 = 0, lambda_L = 0, lambda_0 = 0, Ymax = 0;
    int i = 0, i_sol = 0;
    uint64_t *chromo = calloc(inc,sizeof(uint64_t));

    FILE *profile = fopen("profil.txt","r");

    if(profile == NULL) exit(-1);

    for(i = 0; i < data_l; i++)
    {
        fscanf(profile,"%lf %lf %lf",&data[i][0],&data[i][1],&data[i][2]);
    }
    fclose(profile);

    init(chromo);
    i = 0;
    do
    {
        evaluate(chromo,data,eval);
        survival(&chromo,eval);
        crossover(chromo);
        mutation(chromo);
        printf("%lf     %d\n",best_chromo(eval, &i_sol),i++);
    }
    while(best_chromo(eval, &i_sol) > ferror);

    printf("Algo duration : %d s\n",(unsigned int)(time(0) - duration));

    decode(chromo[i_sol],&Y0,&lambda_L,&lambda_0,&Ymax);
    printf("Ymax = %lf, lambda_0 = %lf, lambda_L = %lf, Y0 = %lf\n",Ymax,lambda_0,lambda_L,Y0);
    
    printf("Which software do you want to use ? 0 - Matlab, 1 - Octave\n");
    scanf("%d",&i);
    plot(i, Y0, lambda_L, lambda_0, Ymax);


    return 0;
}
