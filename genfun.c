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

/**
 * @brief init: initialisation of the first gen chromosomes
 *        randomly generate our constants within their respective intervals and code them into 64 bits integers
 *        all stored in the chromo vector
 * @param chromo: a (inc = 100) length vector
 */
void init(uint64_t* chromo)
{
    int i = 0;
    for(i = 0; i < inc; i++)
    {
        double Y0 = (rand() % (1<<16)) * 1e-6 + inf_Y0;
        double lambda_L = (rand() % (1<<13)) * 1e-3  + inf_lam_L;
        double lambda_0 = (rand() % (1<<19)) * 1e-4 + inf_lam_0;
        double Ymax = (rand() % (1<<16)) * 1e-6 + inf_Ymax;
        chromo[i] = code(Y0,lambda_L,lambda_0,Ymax);
    }
}

/**
 * @brief code: 64 bits integer coding function
 *        the constants are converted into integers (see Constants intervals and conversion rules on genfun.h)
 *        and then all stored in one 64 bits integer
 * @param Y0
 * @param lambda_L
 * @param lambda_0
 * @param Ymax
 * @return 64 bit integer containing the 4 constants [Y0 (16 bits), lambda_L (13 bits), lambda_0 (19 bits), Ymax (16 bits)]
 */
uint64_t code(double Y0, double lambda_L, double lambda_0, double Ymax)
{
    uint16_t a = (uint16_t) ( (Y0       - inf_Y0)    * 1e6 + 1 );
    uint16_t b = (uint16_t) ( (lambda_L - inf_lam_L) * 1e3 + 1 );
    uint32_t c = (uint32_t) ( (lambda_0 - inf_lam_0) * 1e4 + 1 );
    uint16_t d = (uint16_t) ( (Ymax     - inf_Ymax)  * 1e6 + 1 );

    return ((uint64_t)a << (16 + 19 + 13)) + ((uint64_t)b << (16 + 19)) + ((uint64_t)c << 16) + ((uint64_t)d << 0);
}

/**
 * @brief decode: 64 bits integer decoding function
 *        Retreive the information (4 constants) from the 64 bits chromosome using masks and convert them into doubles
 *        (see Constants intervals and conversion rules on genfun.h)
 * @param chromo: the 64 bit chromosome
 * @param Y0
 * @param lambda_L
 * @param lambda_0
 * @param Ymax
 * @return 1
 */
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

/**
 * @brief fitting: calculates the image of lambda under the lorentzien function using the constants in a given chromosome
 *        chromo is convertend into the four constants usinf decode and the values are replaced in the fonction
 *        the image of lambda is then calculated and returned
 * @param chromo: 64 bit chromosome containing the four constants
 * @param lambda
 * @return j(lambda) = Y0 + B / ( (lambda - lambda_0)² + (lambda_L / 2)² ); where B =  (Ymax - Y0) * (lambda_L / 2)²
 */
double fitting(uint64_t chromo, double lambda)
{
    double Y0, lambda_L, lambda_0, Ymax;
    decode(chromo,&Y0,&lambda_L,&lambda_0,&Ymax);

    return Y0 + ((Ymax - Y0) * sqr(lambda_L / 2))/(sqr(lambda-lambda_0) + sqr(lambda_L/2));
}

/**
 * @brief evaluate: function that evaluates the chromosomes
 *        for each chromosome [i], the diffrence between the experimental data points and those calculated using
 *        fitting() function are squared and summed (we get the error). They are stored in the eval vecor
 * @param chromo: vector of chromosomes
 * @param data: data retreived from the profil.txt
 * @param eval: vector of the scores, the results of the evaluation
 */
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

/**
 * @brief crossover: makes the crossover operation
 * @param chromo: vector of chromosomes
 */
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

/**
 * @brief mutation: makes the mutation operation
 *        the Mutation rate 'Mu' represents how many bits will be mutated:
 *          - There is 'inc * c_length' bits (100 * 64 bits)
 *          - 'Mu * (inc * c_length)' bits are randomly chosen to be mutated
 * @param chromo: vector of chromosomes
 */
void mutation(uint64_t *chromo)
{
    int i = 0, R = 0;
    for(i = 0; i < (int)(Mu*inc*c_length); i++)
    {
        R = rand()%(c_length*inc); // random number between 0 and total number of bits 'inc * c_length'

        /* 'R/c_length' is the index of the cromosome and 'R % c_length' is the position of the bit to mutate
         * ex: lets take R = 197. So the bit to change is in the 197/64 = 3 rd chromosome
         * the position of the bit in the 3rd chromo is 197 % 64 = 5;
         * We use a mask (1 << 5 = 0b100000) to toggle the bit using the XOR (^) operator
         * */
        chromo[R/c_length] ^= (uint64_t)1 << (R % c_length);
    }
}

/**
 * @brief best_chromo: returns the error of the best chromosome (the smallest error) and put its index on i_sol
 * @param eval: vector of the evaluations (errors)
 * @param i_sol: index of the best chromosome
 * @return smalest error (smalest value on eval)
 */
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
