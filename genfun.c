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
 *        randomly generate our constants within their respective intervals (see Constants intervals and conversion rules on genfun.h)
 *        and code them into 64 bits integers all stored in the chromo vector
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
 *        chromo is converted into the four constants using decode() and the values are replaced in the function
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

/**
 * @brief survival: the Selection, this function chooses which chromosomes will be part of the next generation
 *        it's done randomly, using a chromosomes' probability to survive. The smallest is the error for a given
 *        chromosome, the bigger is its chance to survive
 * @param chromo: double pointer to the chromo vector (because we need to change the adress it points)
 * @param eval: evaluations vector (errors)
 */
void survival(uint64_t **chromo, double *eval)
{
    uint64_t *temp_c = calloc(inc,sizeof(uint64_t)); //allocated space for the new generation of chromosomes
    double c_prob[inc] = {0}; //Cumulative Probability vector
    double s = 0; //sum of all the chromosomes' errors (that we get using evaluate())
    int i = 0;
    for(i = 0; i < inc; i++)
    {
        c_prob[i] = (double)1/(double)(eval[i]+1); //since we're trying to minimize eval, we also trying to maximise 1/eval
        s += c_prob[i]; //calculate the sum of 1/eval[i]
    }

    c_prob[0] /= s; //so the first probability is (1/eval[0]) / sum(1/eval[i])

    for(i = 1; i < inc; i++)
    {
        c_prob[i] = c_prob[i] / s + c_prob[i - 1]; //so we calculate the cumulative probabilities
    }

    /*now begins the selection, we randomly generate numbers between 0 and 1 stored in 's'
     * the idea is:
     * if 0 < s < c_prob[0] then the 1st chromosome 'chromo[0]' is selected to be in the next gen and is strored in temp_c
     * if c_prob[0] < s < c_prob[1] the the second one is selected 'chromo[1]'
     * if c_prob[j - 1] < s < c_prob[j] the the 'chromo[j]' is selected
     * and we do this until we have all 'inc' number of chromosomes
     */
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

    /*in order to gain some performance and avoide the copy of the entire vector, we just free the old 'chromo' vector
     * and put the 'temp_c' vector adress in 'chromo' to point at it
     */
    free(*chromo);
    *chromo = temp_c;
}

/**
 * @brief crossover: makes the crossover operation
 *        a number of chromosomes are chosen to be crossed according to the cross rate 'X'
 *        the crossing position is randomly chosen
 * @param chromo: vector of chromosomes
 */
void crossover(uint64_t *chromo)
{
    //'parent' is the vector of the chromosomes' indexes chosen to be crossed and 'n' is their number
    int parent[inc] = {0}, i = 0, n = 0;

    /* a random number between 0 and 1 is generated and compared to the 'X' rate, if it is smaller, the
     * chromosome 'i' is selected to be a parent, its index is added to 'parent' and 'n' is incremented
     */
    for(i = 0; i < inc; i++)
    {
        if(rand()/(double)RAND_MAX < X) parent[n++] = i;
    }

    /* The parents are corossed 2-by-2 so the 1st X 2nd, 3nd X 4th, ... (n-1)th X (n)th
     * the position 'R' of the cross is chosen randomly for every pair of chromosomes R = 1:63
     * then we use a mask to get the 4 parts using the AND (&) operator from the parents and bring them together
     * to create 2 new chromosomes that will replace their parents on the 'chromo' vector
     *
     * ex: let's say we have two 4 bits integers we want to cross: c1 = 1101 and c2 = 0110 and R = 2 means we have
     * to part the chromosomes between the bit 1 and 2. mask = (1 << R) - 1 = 0b0100 - 1 = 0b0011
     * so C1 = c1 & mask + c2 & ~mask = 0b1101 & 0b0011 + 0b0110 & 0b1100 = 0b0001 + 0b0100 = 0b0101
     * and C2 = c1 & ~mask + c2 & mask = 0b1101 & 0b1100 + 0b0110 & 0b0011 = 0b1100 + 0b0010 = 0b1110
     */
    for(i = 0; i < n; i+=2)
    {
        int R = rand() % (c_length - 1) + 1;
        uint32_t mask = (1 << R) - 1, temp;
        temp = chromo[parent[i]] & mask;
        chromo[parent[i]] = (chromo[parent[i]] & ~mask) + (chromo[parent[i+1]] & mask);
        chromo[parent[i+1]] = (chromo[parent[i+1]] & ~mask) + temp;
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
         * ex: lets take R = 197. So the bit to change is in the chromosome f index 197/64 = 3
         * the position of the bit in the chromo[3] is 197 % 64 = 5;
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
