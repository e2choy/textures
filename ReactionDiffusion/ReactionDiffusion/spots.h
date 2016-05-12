#pragma once

#define MAX_SIZE 200
extern int interval;
extern int value_switch;
extern int stripe_flag;

#define  STRIPES  1
#define  SPOTS    2
#define  LEOPARD  3  
extern float a[MAX_SIZE][MAX_SIZE];
extern float b[MAX_SIZE][MAX_SIZE];
extern float c[MAX_SIZE][MAX_SIZE];
extern float d[MAX_SIZE][MAX_SIZE];
extern float e[MAX_SIZE][MAX_SIZE];

extern float da[MAX_SIZE][MAX_SIZE];
extern float db[MAX_SIZE][MAX_SIZE];
extern float dc[MAX_SIZE][MAX_SIZE];
extern float dd[MAX_SIZE][MAX_SIZE];
extern float de[MAX_SIZE][MAX_SIZE];

extern float ai[MAX_SIZE][MAX_SIZE];
extern float p1, p2, p3;
extern float diff1, diff2;
extern float arand;
extern float a_steady;
extern float b_steady;
extern float beta_init;
extern float beta_rand;
extern float speed;
extern int sim;
extern int xsize;
extern int ysize;

void do_stripes();
void do_spots();
void compute();
void semi_equilibria();
void show( float values[MAX_SIZE][MAX_SIZE] );
void multiplicative_help();
void turing();
void rescale_values( float values[MAX_SIZE][MAX_SIZE], float min_final, float max_final );