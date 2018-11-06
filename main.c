#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#define inc 6
#define X 0.8
#define Mu 0.5
#define c_length 20
#define abs(x) ((x < 0) ? (-1*x):(x))

int decode(uint32_t chromo, int* a, int* b, int* c, int* d);
int fitting(uint32_t c);
uint32_t code(int a, int b, int c, int d);
void init(uint32_t* chromo);
void evaluate(uint32_t *chromo, int* eval);
void survival(uint32_t **chromo, int* eval);
void crossover(uint32_t *chromo);
void mutation(uint32_t *chromo);
int best(int *eval, int *i_sol);

int main()
{
    srand((unsigned int)time(0));

    uint32_t *chromo = calloc(inc,sizeof(uint32_t));
    int eval[inc] = {0}, i = 0, i_sol = 0;
    init(chromo);

    do
    {
        evaluate(chromo,eval);
        survival(&chromo,eval);
        crossover(chromo);
        mutation(chromo);
        printf("%d..... %d  %d  %d  %d  %d  %d\n\n",abs(-i),eval[0],eval[1],eval[2],eval[3],eval[4],eval[5]);
    }
    while(i++ < 10000 && best(eval,&i_sol) > 0);


    int a = 0, b = 0, c = 0, d = 0;
    decode(chromo[i_sol],&a,&b,&c,&d);
    printf("sol = %d, %d, %d, %d\n",a,b,c,d);
    return 0;
}

void init(uint32_t* chromo)
{
    chromo[0] = code(12,5,23,8);
    chromo[1] = code(2,21,18,3);
    chromo[2] = code(10,4,13,14);
    chromo[3] = code(20,1,10,6);
    chromo[4] = code(1,4,13,19);
    chromo[5] = code(20,5,17,1);
}

uint32_t code(int a, int b, int c, int d)
{
    return (a << 15) + (b << 10) + (c << 5) + (d << 0);
}

int decode(uint32_t chromo, int* a, int* b, int* c, int* d)
{
    *a = (chromo & (0b11111<<15)) >> 15;
    *b = (chromo & (0b11111<<10)) >> 10;
    *c = (chromo & (0b11111<<5)) >> 5;
    *d = (chromo & (0b11111<<0)) >> 0;
    return 1;
}

int fitting(uint32_t chromo)
{
    int a,b,c,d;
    decode(chromo,&a,&b,&c,&d);
    return abs(a + 2*b + 3*c + 4*d - 30);
}

void evaluate(uint32_t* chromo, int* eval)
{
    int i = 0;
    for(i = 0;i < inc; i++)
    {
        eval[i] = abs(fitting(chromo[i]));
    }
}

void survival(uint32_t **chromo, int* eval)
{
    uint32_t *temp_c = calloc(inc,sizeof(uint32_t));
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

void crossover(uint32_t *chromo)
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

void mutation(uint32_t *chromo)
{
    int i = 0, R = 0;
    for(i = 0; i < (int)(Mu*inc*c_length); i++)
    {
        R = rand()%(c_length/inc);
        chromo[R/c_length] ^= 1 << (R % c_length);
    }
}

int best(int *eval,int *i_sol)
{
    int min = abs(eval[0]),i = 0;
    for(i = 1; i < inc; i++)
    {
        if(min > abs(eval[i]))
        {
            min = eval[i];
            *i_sol = i;
        }
    }
    return min;
}
