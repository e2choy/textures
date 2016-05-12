#include "stdafx.h"
#include "leopard.h"
#include "spots.h"
#include <cstdio>
//////////////////////////////////////////////////////////////////////////
void do_leopard(){
    beta_init = 12.0f;
    //beta_rand = 0.1f;
    beta_rand = 0.2f;

    a_steady = 4;
    b_steady = 4;

    diff1 = 0.25f;
    diff2 = 0.0625f;

    p1 = 0.02f;     //set this for larger spots
    //p1 = 0.05;    //default leopard spot
    //p1 = 0.1f;    //default 0.2
    p2 = 0.0f;
    p3 = 0.0f;

    sim = LEOPARD;
    value_switch = 2;
    compute_leopard();
}
//////////////////////////////////////////////////////////////////////////
static Mat frozenMat;
void run_second_stage();
void turingFrozen( );
void compute_leopard()
{
    int k;
    //int iterations = 99999999;
    int iterations = 50000; // 99999999;

    /* calculate semistable equilibria */
    semi_equilibria();

    /* start things diffusing */
    cout << "stage1" << endl;
    for ( k = 0; k < iterations; k++ ) {
        if ( k % interval == 0 ) {
            printf( "iteration %d\n", k );
        }
        turing();
    }

    //iterate through to detect frozen
    cout << ysize << "," << xsize << endl;
    frozenMat = Mat( ysize, xsize, CV_32FC1 );
    for ( int i = 0; i < xsize; ++i ){
        for ( int j = 0; j < ysize; ++j ){
            float aval = a[i][j];
            float bval = b[i][j];
            if ( bval < 4.0000f && bval >= 0.0f){
                frozenMat.at<float>( j, i ) = 1.0f;
            }
            else{
                frozenMat.at<float>( j, i ) = 0.0f;
            }
        }
    }
    cv::imshow( "frozenMat", frozenMat );

    //iterate through second simulation
    run_second_stage();

    //freeze the initial iteration
    //chemical concentrations of b
    show( b );
}

void run_second_stage(){
    cout << "stage2" << endl;
    //change the value of p1
    p1 = 4.0f * p1;
    semi_equilibria();
    //set the chemical concentrations of a and b
    for ( int i = 0; i < xsize; ++i ){
        for ( int j = 0; j < ysize; ++j ){
            if ( frozenMat.at<float>( j, i ) > 0.0f ){
                a[i][j] = 4.0f;
                b[i][j] = 4.0f;
            }
        }
    }

    int k;
    int iterations = 50000; // 99999999;
    for ( k = 0; k < iterations; k++ ) {
        if ( k % interval == 0 ) {
            printf( "iteration %d\n", k );
        }
        turingFrozen();
    }

    //iterate again
    Mat test = Mat( ysize, xsize, CV_32FC1 );
    for ( int i = 0; i < xsize; ++i ){
        for ( int j = 0; j < ysize; ++j ){
            float aval = a[i][j];
            float bval = b[i][j];
            if ( bval < 4.0000f && bval >= 0.0f && frozenMat.at<float>(j,i) == 0.0f){
                test.at<float>( j, i ) = 1.0f;
            }
            else{
                test.at<float>( j, i ) = 0.0f;
            }
        }
    }
    cv::imshow( "testMat", test );

    //merge the frozen with the test
    Mat finalImage = Mat( ysize, xsize, CV_8UC4 );
    for ( int i = 0; i < xsize; ++i ){
        for ( int j = 0; j < ysize; ++j ){
            if ( frozenMat.at<float>( j, i ) == 1.0f ){
                finalImage.at<Vec4b>( j, i ) = Vec4b( 128, 128, 128, 255 );
            }
            else if ( test.at<float>( j, i ) == 1.0f ){
                finalImage.at<Vec4b>( j, i ) = Vec4b( 0, 0, 0, 255 );
            }
            else{
                finalImage.at<Vec4b>( j, i ) = Vec4b( 255, 255, 255, 255 );
            }
        }
    }
    cv::imshow( "finalImage", finalImage );
    cv::imwrite( "leopard.png", finalImage );

}

void turingFrozen( ){
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
    for ( i = 0; i < xsize; i++ ){
        for ( j = 0; j < ysize; j++ ) {
            if ( frozenMat.at<float>( j, i ) == 0.0f ){
                a[i][j] += (speed * da[i][j]);
                b[i][j] += (speed * db[i][j]);
                if ( b[i][j] < 0 ){ b[i][j] = 0; }
            }
        }
    }
}