/************************************************************************************************
 * L3 EEA - Projet Techniques Scientifiques
 * Theme : Algorithme Genetique
 * Auteurs:
 *      - AIT MOUHOUB Farid
 *      - AOUCI Sofiane
 *      - JAQUET E.Prince
 *      - LAVAL Hugo
 ************************************************************************************************/

#include "nongenfun.h"

/**
 * @brief sqr: calculates the scare of a number
 * @param x
 * @return x^2
 */
double sqr(double x){
    return x*x;
}

/**
 * @brief plot: uses matlab or octave to plot the curves. The script is run using the system command lines
 * @param soft: an integer, if it's 0, we use matlab otherwise, we use octave
 * @param Y0
 * @param lambda_L
 * @param lambda_0
 * @param Ymax
 */
void plot(int soft, double Y0, double lambda_L, double lambda_0, double Ymax)
{
    char script[200], dir[100];
    getcwd(dir,sizeof(dir));
    if(soft) sprintf(script,"%s\\octave --persist --eval \"cd %s; c = [%lf, %lf, %lf, %lf]; plot_result\"",octave,dir,Ymax,lambda_0,lambda_L,Y0);
    else sprintf(script,"matlab -nosplash -nodesktop -r \"cd %s; c = [%lf, %lf, %lf, %lf]; plot_result\"",dir,Ymax,lambda_0,lambda_L,Y0);
    system(script);

}
