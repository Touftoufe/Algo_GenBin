/************************************************************************************************
 * L3 EEA - Projet Techniques Scientifiques
 * Theme : Algorithme Genetique
 * Auteurs:
 *      - AIT MOUHOUB Farid
 *      - AOUCI Sofiane
 *      - JAQUET E.Prince
 *      - LAVAL Hugo
 ************************************************************************************************/

/*
 * lambda_0 E [6555, 6570] codee en 19 bit chacun donc [0, 524287] <=> [6555, 6607,4287]
 * lambda_L E [7, 15.191] codee en 13 bit chacun donc [0, 8191] <=> [7, 15.191]
 * Y_0 E [0.2, 0.25] codee en 16 bit donc [0, 65535] <=> [0.2, 0.265535]
 * Y_max E [0.72, 0.785535] codee en 16 bit donc [0, 65535]
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

//#define s_abs(x) ((x) > (0) ? (x):(-x))
//#define sqr(x) (x * x)

double s_abs(double x){
    if(x > 0) return x;
    else return -x;
}

double sqr(double x){
    return x*x;
}

#define pi 3.14159265358979323846

#define inf_lam_0 6555
#define inf_lam_L 7
#define inf_Y0 0.2
#define inf_Ymax 0.72
#define data_l 1024

#define inc 100
#define X 0.15
#define Mu 0.01
#define c_length 64

void init(uint64_t* chromo);
int decode(uint64_t chromo, double* Y0, double* lambda_L, double* lambda_0, double* Ymax);
uint64_t code(double Y0, double lambda_L, double lambda_0, double Ymax);
double fitting(uint64_t c, double lambda);
void evaluate(uint64_t *chromo, double data[][3], double *eval);
void survival(uint64_t **chromo, double* eval);
void crossover(uint64_t *chromo);
void mutation(uint64_t *chromo);
double soso(double *eval, int *i_sol);
uint64_t grayInverse(uint64_t n);

int main()
{
    srand((unsigned int)time(0));
    FILE *profile = fopen("profil.txt","r");
    double data[data_l][3];
    int i = 0;
    if(profile == NULL) exit(-1);

    for(i = 0; i < data_l; i++)
    {
        fscanf(profile,"%lf %lf %lf",&data[i][0],&data[i][1],&data[i][2]);
    }
    fclose(profile);

    uint64_t *chromo = calloc(inc,sizeof(uint64_t));
    int i_sol = 0;
    double eval[inc] = {0};
    init(chromo);
    i = 0;
    do
    {
        evaluate(chromo,data,eval);
        survival(&chromo,eval);
        crossover(chromo);
        mutation(chromo);
        //printf("%d..... %d  %d  %d  %d  %d  %d\n\n",s_abs(-i),eval[0],eval[1],eval[2],eval[3],eval[4],eval[5]);
        printf("%lf\n",soso(eval, &i_sol));
    }
    while(soso(eval, &i_sol) > 0.204);

    double a = 0, b = 0, c = 0, d = 0;
    decode(chromo[i_sol],&a,&b,&c,&d);
    printf("%lf, %lf, %lf, %lf\n",d,c,b,a);
    //printf("-------- %d --------- ",grayInverse(15 ^ (15 >> 1)));

    return 0;
}

void init(uint64_t* chromo)
{
    int i = 0;
    for(i = 0; i < inc; i++)
    {
        double a = (rand() % 65535)* 1e-6 + inf_Y0;
        double b = (rand() % 8191)* 1e-3  + inf_lam_L;
        double c = (rand() % 524287)* 1e-4 + inf_lam_0;
        double d = (rand() % 65535)* 1e-6 + inf_Ymax;
        chromo[i] = code(a,b,c,d);
    }
}

uint64_t code(double Y0, double lambda_L, double lambda_0, double Ymax)
{
    uint64_t bin;
    uint16_t a =  (Y0 - inf_Y0) * 1e6 + 1 ;
    uint16_t b =  (lambda_L - inf_lam_L) *1e3 + 1 ;
    uint32_t c =  (lambda_0 - inf_lam_0) * 1e4 + 1 ;
    uint16_t d =  (Ymax - inf_Ymax) * 1e6 + 1 ;
    bin = ((uint64_t)a << (16 + 19 + 13)) + ((uint64_t)b << (16 + 19)) + ((uint64_t)c << 16) + ((uint64_t)d << 0);
    return bin;// ^ (bin >> 1);
}

int decode(uint64_t chromo, double *Y0, double *lambda_L, double *lambda_0, double *Ymax)
{
    //chromo = grayInverse(chromo);
    uint16_t a =  (chromo & ((((uint64_t)0b1<<16) - 1) << (16 + 19 + 13)) ) >> (16 + 19 + 13);
    uint16_t b =  (chromo & ((((uint64_t)0b1<<13) - 1) << (16 + 19)) ) >> (16 + 19);
    uint32_t c =  (chromo & ((((uint64_t)0b1<<19) - 1) <<  16) ) >>  16;
    uint16_t d =  (chromo & ((((uint64_t)0b1<<16) - 1) << 0) ) >> 0;

    *Y0 = a * 1e-6 + inf_Y0;
    *lambda_L = b * 1e-3 + inf_lam_L;
    *lambda_0 = c * 1e-4 + inf_lam_0;
    *Ymax = d * 1e-6 + inf_Ymax;
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
        s = (double)rand()/(double)RAND_MAX;

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
        if((double)rand()/(double)RAND_MAX < X) parent[n++] = i;
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

double soso(double *eval, int *i_sol)
{
    double min = eval[0];
    int i = 0;
    for(i = 1; i < inc; i++)
    {
        //printf("%lf     %lf\n",min,eval[i]);
        if(min > eval[i])
        {
            min = eval[i];
            *i_sol = i;
        }
    }
    return min;
}

uint64_t grayInverse(uint64_t n) {
        uint64_t translate, idiv;
        translate = 1;
        while(1) {
            idiv = n >> translate;
            n ^= idiv;
            if (idiv <= 1 || translate == 32)
                return n;
            translate <<= 1; // double le nb de shifts la prochaine fois
        }
  }
