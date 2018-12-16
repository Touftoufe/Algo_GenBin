/************************************************************************************************
 * L3 EEA - Projet Techniques Scientifiques
 * Theme : Algorithme Genetique
 * Auteurs:
 *      - AIT MOUHOUB Farid
 *      - AOUCI Sofiane
 *      - JAQUET E.Prince
 *      - LAVAL Hugo
 ************************************************************************************************/

#include "genfun.h"


void init(uint64_t* chromo)
{
    int i = 0;
    for(i = 0; i < inc; i++)
    {
        double Y0 = (rand() % (1<<16))* 1e-6 + inf_Y0;
        double lambda_L = (rand() % (1<<13))* 1e-3  + inf_lam_L;
        double lambda_0 = (rand() % (1<<19))* 1e-4 + inf_lam_0;
        double Ymax = (rand() % (1<<16))* 1e-6 + inf_Ymax;
        chromo[i] = code(Y0,lambda_L,lambda_0,Ymax);
    }
}

uint64_t code(double Y0, double lambda_L, double lambda_0, double Ymax)
{
    uint64_t bin;
    uint16_t a = (uint16_t) ( (Y0       - inf_Y0)    * 1e6 + 1 );
    uint16_t b = (uint16_t) ( (lambda_L - inf_lam_L) * 1e3 + 1 );
    uint32_t c = (uint32_t) ( (lambda_0 - inf_lam_0) * 1e4 + 1 );
    uint16_t d = (uint16_t) ( (Ymax     - inf_Ymax)  * 1e6 + 1 );
    bin = ((uint64_t)a << (16 + 19 + 13)) + ((uint64_t)b << (16 + 19)) + ((uint64_t)c << 16) + ((uint64_t)d << 0);

    return bin;
}

int decode(uint64_t chromo, double *Y0, double *lambda_L, double *lambda_0, double *Ymax)
{
    uint16_t a =  (chromo & ((((uint64_t)0b1<<16) - 1) << (16 + 19 + 13)) ) >> (16 + 19 + 13);
    uint16_t b =  (chromo & ((((uint64_t)0b1<<13) - 1) << (16 + 19)) ) >> (16 + 19);
    uint32_t c =  (chromo & ((((uint64_t)0b1<<19) - 1) <<  16) ) >>  16;
    uint16_t d =  (chromo & ((((uint64_t)0b1<<16) - 1) << 0) ) >> 0;

    *Y0       = a * 1e-6 + inf_Y0;
    *lambda_L = b * 1e-3 + inf_lam_L;
    *lambda_0 = c * 1e-4 + inf_lam_0;
    *Ymax     = d * 1e-6 + inf_Ymax;

    return 1;
}

double fitting(uint64_t chromo, double lambda)
{
    double Y0, lambda_L, lambda_0, Ymax;
    decode(chromo,&Y0,&lambda_L,&lambda_0,&Ymax);

    return Y0 + ((Ymax - Y0) * sqr(lambda_L) / 4)/(sqr(lambda-lambda_0) + sqr(lambda_L/2));
}

void evaluate(uint64_t* chromo,double data[][3], double* eval)
{
    int i = 0, j = 0;
    for(i = 0;i < inc; i++)
    {
        eval[i] = 0;
        for(j = 0; j < data_l; j++)
        {
            eval[i] += sqr(data[j][1] - fitting(chromo[i],data[j][0]));
        }
    }
}

void survival(uint64_t **chromo, double *eval)
{
    uint64_t *temp_c = calloc(inc,sizeof(uint64_t));
    double c_prob[inc] = {0}, s = 0;
    int i = 0;
    for(i = 0; i < inc; i++)
    {
        c_prob[i] = (double)1/(double)(eval[i]+1);
        s += c_prob[i];
    }

    c_prob[0] /= s;

    for(i = 1; i < inc; i++)
    {
        c_prob[i] = c_prob[i] / s + c_prob[i - 1];
    }

    for(i = 0; i < inc; i++)
    {
        int j = 0;
        s = rand()/(double)RAND_MAX;

        if(s <= c_prob[0] && s > 0)
        {
            temp_c[i] = (*chromo)[0];
        }
        else
        {
            j++;
            while(j<inc)
            {
                if(s <= c_prob[j] && s > c_prob[j-1])
                {
                    temp_c[i] = (*chromo)[j];
                    j = inc;
                }
                j++;
            }
        }
    }

    free(*chromo);
    *chromo = temp_c;
}

void crossover(uint64_t *chromo)
{
    int parent[inc] = {0}, i = 0, n = 0;
    for(i = 0; i < inc; i++)
    {
        if(rand()/(double)RAND_MAX < X) parent[n++] = i;
    }

    for(i = 0; i < n; i++)
    {
        int R = rand() % (c_length - 1) + 1;
        uint32_t filter = (1 << R) - 1, temp;
        temp = chromo[parent[i]] & filter;
        chromo[parent[i]] = (chromo[parent[i]] & ~filter) + (chromo[parent[(i+1)%n]] & filter);
        chromo[parent[(i+1)%n]] = (chromo[parent[(i+1)%n]] & ~filter) + temp;
    }

}

void mutation(uint64_t *chromo)
{
    int i = 0, R = 0;
    for(i = 0; i < (int)(Mu*inc*c_length); i++)
    {
        R = rand()%(c_length*inc);
        chromo[R/c_length] ^= (uint64_t)1 << (R % c_length);
    }
}

double best_chromo(double *eval, int *i_sol)
{
    double min = eval[0];
    *i_sol = 0;
    int i = 0;
    for(i = 1; i < inc; i++)
    {
        if(min > eval[i])
        {
            min = eval[i];
            *i_sol = i;
        }
    }

    return min;
}
