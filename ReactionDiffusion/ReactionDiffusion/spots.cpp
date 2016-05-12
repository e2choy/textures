/*

Make spots and stripes with reaction-diffusion.

The spot-formation system is described in the article:

"A Model for Generating Aspects of Zebra and Other Mammailian
Coat Patterns"
Jonathan B. L. Bard
Journal of Theoretical Biology, Vol. 93, No. 2, pp. 363-385
(November 1981)

The stripe-formation system is described in the book:

Models of Biological Pattern Formation
Hans Meinhardt
Academic Press, 1982


Permission is granted to modify and/or distribute this program so long
as the program is distributed free of charge and this header is retained
as part of the program.

Copyright (c) Greg Turk, 1991

*/

#include <stdio.h>
#include <math.h>
#include "stdafx.h"
#include "spots.h"
#include "leopard.h"

using namespace std;
int xsize;
int ysize;

/******************************************************************************
Pick a random number between min and max.
******************************************************************************/
float frand( float min, float max )
{
    float f = (float )rand() / RAND_MAX;
    return min + f * (max - min);
}

Mat m_image;
float frand();

/* simulation variables */
int interval = 200;
int value_switch = 1;
int stripe_flag = 0;



float a[MAX_SIZE][MAX_SIZE];
float b[MAX_SIZE][MAX_SIZE];
float c[MAX_SIZE][MAX_SIZE];
float d[MAX_SIZE][MAX_SIZE];
float e[MAX_SIZE][MAX_SIZE];

float da[MAX_SIZE][MAX_SIZE];
float db[MAX_SIZE][MAX_SIZE];
float dc[MAX_SIZE][MAX_SIZE];
float dd[MAX_SIZE][MAX_SIZE];
float de[MAX_SIZE][MAX_SIZE];

float ai[MAX_SIZE][MAX_SIZE];
float p1, p2, p3;
float diff1, diff2;
float arand;
float a_steady;
float b_steady;
float beta_init;
float beta_rand;
float speed = 1.0;
int sim = 1;


/******************************************************************************
Main routine.
******************************************************************************/

int main( int argc, char* argv[] )
{
    char *s;

    /* parse the command line options */
    xsize = MAX_SIZE;
    ysize = MAX_SIZE;
    
    int mode = 0;   //0 spots, 1 stripes, 2 leopard
    m_image = Mat( xsize, ysize, CV_8UC4, cv::Scalar( 0, 0, 0, 0 ) );

    if ( mode == 0 ){
        do_spots();
    }
    else if ( mode == 1 ){
        do_stripes();
    }
    else if ( mode == 2 ){
        do_leopard();
    }
    cv::imshow( "reaction-diffusion", m_image );
    cv::waitKey();
}


/******************************************************************************
Run Meinhardt's stripe-formation system.
******************************************************************************/

void do_stripes()
{
    float scale = 0.10f;
    p1 = 0.04f * scale;
    p2 = 0.06f * scale;
    p3 = 0.04f * scale;

    diff1 = 0.009f;
    diff2 = 0.2f;

    arand = 0.02f;

    sim = STRIPES;
    value_switch = 1;
    compute();
}


/******************************************************************************
Run Turing reaction-diffusion system.
******************************************************************************/
void do_spots()
{
    beta_init = 12.0f;
    beta_rand = 0.1f;

    a_steady = 4;
    b_steady = 4;

    diff1 = 0.25f;
    diff2 = 0.0625f;

    p1 = 0.02f; //set this for large spots
    //p1 = 0.2f;
    p2 = 0.0f;
    p3 = 0.0f;

    sim = SPOTS;
    value_switch = 2;
    compute();
}



/******************************************************************************
Diffuse and react.
******************************************************************************/
void compute()
{
    int k;
    //int iterations = 99999999;
    int iterations = 100000; // 99999999;

    /* calculate semistable equilibria */
    semi_equilibria();

    /* start things diffusing */
    for ( k = 0; k < iterations; k++ ) {
        if ( k % interval == 0 ) {
            printf( "iteration %d\n", k );
        }
        /* perform reaction and diffusion */
        switch ( sim ) {
        case STRIPES:
            multiplicative_help();
            break;
        case SPOTS:
            turing();
            break;
        default: break;
        }
    }

    switch ( value_switch ) {
    case 1:
        show( a );
        break;
    case 2:
        show( b );
        break;
    case 3:
        show( c );
        break;
    case 4:
        show( d );
        break;
    case 5:
        show( e );
        break;
    default:
        printf( "bad switch in compute: %d\n", value_switch );
        break;
    }
}


/******************************************************************************
Create stripes with what Hans Meinhardt calls a two-species balance.
******************************************************************************/

void multiplicative_help()
{
    int i, j;
    int iprev, inext, jprev, jnext;
    float aval, bval, cval, dval, eval;
    float ka, kc, kd;
    float temp1, temp2;
    float dda, ddb;
    float ddd, dde;

    /* compute change in each cell */
    for ( i = 0; i < xsize; i++ ) {
        ka = -p1 - 4 * diff1;
        kc = -p2;
        kd = -p3 - 4 * diff2;
        iprev = (i + xsize - 1) % xsize;
        inext = (i + 1) % xsize;
        for ( j = 0; j < ysize; j++ ) {
            jprev = (j + ysize - 1) % ysize;
            jnext = (j + 1) % ysize;

            aval = a[i][j];
            bval = b[i][j];
            cval = c[i][j];
            dval = d[i][j];
            eval = e[i][j];

            temp1 = 0.01 * aval * aval * eval * ai[i][j];
            temp2 = 0.01 * bval * bval * dval;

            dda = a[i][jprev] + a[i][jnext] + a[iprev][j] + a[inext][j];
            ddb = b[i][jprev] + b[i][jnext] + b[iprev][j] + b[inext][j];
            ddd = d[i][jprev] + d[i][jnext] + d[iprev][j] + d[inext][j];
            dde = e[i][jprev] + e[i][jnext] + e[iprev][j] + e[inext][j];

            da[i][j] = aval * ka + diff1 * dda + temp1 / cval;
            db[i][j] = bval * ka + diff1 * ddb + temp2 / cval;
            dc[i][j] = cval * kc + temp1 + temp2;
            dd[i][j] = dval * kd + diff2 * ddd + p3 * aval;
            de[i][j] = eval * kd + diff2 * dde + p3 * bval;
        }
    }

    /* affect change */

    for ( i = 0; i < xsize; i++ )
        for ( j = 0; j < ysize; j++ ) {
        a[i][j] += (speed * da[i][j]);
        b[i][j] += (speed * db[i][j]);
        c[i][j] += (speed * dc[i][j]);
        d[i][j] += (speed * dd[i][j]);
        e[i][j] += (speed * de[i][j]);
        }
}


/******************************************************************************
Turing's reaction-diffusion equations.
******************************************************************************/
void turing()
{
    int i, j;
    int iprev, inext, jprev, jnext;
    float aval, bval;
    float ka;
    float dda, ddb;
    float Diff1, Diff2;

    Diff1 = diff1 / 2.0;
    Diff2 = diff2 / 2.0;
    ka = p1 / 16.0;

    /* compute change in each cell */
    for ( i = 0; i < xsize; i++ ) {
        iprev = (i + xsize - 1) % xsize;
        inext = (i + 1) % xsize;
        for ( j = 0; j < ysize; j++ ) {
            jprev = (j + ysize - 1) % ysize;
            jnext = (j + 1) % ysize;
            aval = a[i][j];
            bval = b[i][j];
            dda = a[i][jprev] + a[i][jnext] + a[iprev][j] + a[inext][j] - 4 * aval;
            ddb = b[i][jprev] + b[i][jnext] + b[iprev][j] + b[inext][j] - 4 * bval;
            da[i][j] = ka * (16 - aval * bval) + Diff1 * dda;
            db[i][j] = ka * (aval * bval - bval - c[i][j]) + Diff2 * ddb;
        }
    }

    /* affect change */
    for ( i = 0; i < xsize; i++ )
        for ( j = 0; j < ysize; j++ ) {
        a[i][j] += (speed * da[i][j]);
        b[i][j] += (speed * db[i][j]);
        if ( b[i][j] < 0 )
            b[i][j] = 0;
        }
}


/******************************************************************************
Calculate semi-stable equilibria.
******************************************************************************/
void semi_equilibria()
{
    int i, j;
    float ainit, binit;
    float cinit, dinit, einit;

    ainit = binit = cinit = dinit = einit = 0;

    /* figure the values */
    switch ( sim ) {
    case STRIPES:
        for ( i = 0; i < xsize; i++ ) {

            ainit = p2 / (2 * p1);
            binit = ainit;
            cinit = 0.02 * ainit * ainit * ainit / p2;
            dinit = ainit;
            einit = ainit;

            for ( j = 0; j < ysize; j++ ) {
                a[i][j] = ainit;
                b[i][j] = binit;
                c[i][j] = cinit;
                d[i][j] = dinit;
                e[i][j] = einit;
                ai[i][j] = 1 + frand( -0.5 * arand, 0.5 * arand );
            }
        }
        break;

    case SPOTS:
    case LEOPARD:
        for ( i = 0; i < xsize; i++ )
            for ( j = 0; j < ysize; j++ ) {
            a[i][j] = a_steady;
            b[i][j] = b_steady;
            c[i][j] = beta_init + frand( -beta_rand, beta_rand );
            }
        break;

    default:
        printf( "bad case in semi_equilibria\n" );
        break;
    }
}


/******************************************************************************
Switch for picking array to rescale.
******************************************************************************/

void do_rescale( int index, float min, float max )
{
    switch ( index ) {
    case 1:
        rescale_values( a, min, max );
        break;
    case 2:
        rescale_values( b, min, max );
        break;
    case 3:
        rescale_values( c, min, max );
        break;
    case 4:
        rescale_values( d, min, max );
        break;
    case 5:
        rescale_values( e, min, max );
        break;
    default:
        printf( "bad switch in do_rescale: %d\n", index );
        break;
    }
}


/******************************************************************************
Rescale values in array.

Entry:
values    - array to rescale
min_final - minimum value to map to
max_final - maximum value to map to
******************************************************************************/

void rescale_values( float values[MAX_SIZE][MAX_SIZE], float min_final, float max_final )
{
    int i, j;
    float val;
    float min = 1e20;
    float max = -1e20;

    /* find minimum and maximum values */

    for ( i = 0; i < xsize; i++ )
        for ( j = 0; j < ysize; j++ ) {
        if ( values[i][j] < min )
            min = values[i][j];
        if ( values[i][j] > max )
            max = values[i][j];
        }

    if ( min == max ) {
        min = max - .001;
        max = min + .002;
    }

    /* rescale the values */

    for ( i = 0; i < xsize; i++ )
        for ( j = 0; j < ysize; j++ ) {
        val = (values[i][j] - min) / (max - min);
        val = min_final + val * (max_final - min_final);
        values[i][j] = val;
        }
}


/******************************************************************************
Display the activator.
******************************************************************************/

void show( float values[MAX_SIZE][MAX_SIZE] )
{
    int i, j;
    float output;
    float min = 1e20;
    float max = -1e20;

    /* find minimum and maximum values */

    for ( i = 0; i < xsize; i++ )
        for ( j = 0; j < ysize; j++ ) {
        if ( values[i][j] < min )
            min = values[i][j];
        if ( values[i][j] > max )
            max = values[i][j];
        }

    if ( min == max ) {
        min = max - 1;
        max = min + 2;
    }

    printf( "min max diff: %f %f %f\n", min, max, max - min );

    /* display the values */

    for ( i = 0; i < xsize; i++ ){
        for ( j = 0; j < ysize; j++ ) {
            output = (values[i][j] - min) / (max - min);
            if ( output > 0.5f )
                output = 255.0f;
            else
                output = 0.0f;
            m_image.at<Vec4b>( j, i ) = Vec4b( output, output, output, 255 );
        }
    }

    cv::imwrite( "output_image.png", m_image );
}



