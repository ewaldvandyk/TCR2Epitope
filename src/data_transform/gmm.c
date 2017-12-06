#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>


//gsl_vectors and matrices to be filled with the data.
static int dim, numMix;
 
static gsl_vector *g_gaussCoeffs, *g_mixPriors;
static gsl_matrix *g_mus, *g_sigInv;
static gsl_vector *x, *x_sqr, *g_vec, *gaussP;


/* replaces the python function of the same name in Python */
void c_vec_2_posterior(double *invec, double *outvec)
{
    //local vars
    int i,bI; //loop variables
    double pwr, postNorm;    //temporary variable
    
    //copy data from p_vec into gsl vector
    for(i=0;i<dim;i++) 
    {
        gsl_vector_set(g_vec,i,invec[i]);
    }
    
          
    //initialise zero-filled matrix that will contain the final data
    /*** do the math ***/
    for(bI=0; bI<numMix; bI++)
    {
        gsl_vector_memcpy(x,g_vec);//copy contents of vec into vector x (a reset if you will)
        
        gsl_vector_view view1=gsl_matrix_row(g_mus,bI); //store vector slice of matrix
        gsl_vector_sub(x,&view1.vector); //subtract mus from x; store results in x
        
        gsl_matrix_view view2=gsl_matrix_submatrix(g_sigInv,bI*dim,0,dim, dim); //store piece of matrix
        gsl_blas_dgemv(CblasTrans,1.,&view2.matrix,x,0.,x_sqr); //multiply sigInv submatrix to x
        
        gsl_blas_ddot(x_sqr,x,&pwr); //multiply x_sqr to x, store result in pwr
        
        pwr=exp(pwr*-0.5);
        
        gsl_vector_set(gaussP,bI,gsl_vector_get(g_gaussCoeffs,bI)*pwr);         //fill a position in GaussP with a multiplied value
    }
       
    
    gsl_blas_ddot(gaussP,g_mixPriors,&postNorm);
    
    postNorm=1./postNorm;
    
    gsl_vector_mul(gaussP,g_mixPriors);
    
    gsl_vector_scale(gaussP,postNorm);
    
    /********************/
    
    //copy final data into outvec
    for(i=0; i<numMix; i++) {
        outvec[i]=gsl_vector_get(gaussP,i);
//         printf("%d: %e\n", i, outvec[i]);
    }
    return;
    
    
}


//function that receives the fixed data from python and stores it in gsl vectors and matrices
void init(int p_numMix,int p_dim)
{
    //allocate memory for static model vectors
    g_gaussCoeffs=gsl_vector_calloc(p_numMix);
    g_mixPriors=gsl_vector_calloc(p_numMix);
    
    g_mus=gsl_matrix_calloc(p_numMix,p_dim);
    g_sigInv=gsl_matrix_calloc(p_numMix*p_dim,p_dim);
//     printf("g_sigInv length = %ld\n", (long)(p_dim*p_dim)*p_numMix);
    
    //allocate memory for static non-model vectors
    x=gsl_vector_alloc(p_dim); 
    x_sqr=gsl_vector_alloc(p_dim);
    g_vec=gsl_vector_calloc(p_dim);
    gaussP=gsl_vector_calloc(p_numMix);
    
    numMix=p_numMix;
    dim=p_dim;
    return;
}

void set_mixPriors(double *p_mixPriors){
    int i;

    for (i = 0; i < numMix; i++){
        gsl_vector_set(g_mixPriors,i,p_mixPriors[i]);
    }
    return;
}

void set_gaussCoeffs(double *p_gaussCoeffs){
    int i;
    
    for (i = 0; i < numMix; i++){
        gsl_vector_set(g_gaussCoeffs,i,p_gaussCoeffs[i]);
    }
    return;
}

void set_mus(double *p_mus){
    int i,j;
    
    for(i=0; i<numMix;i++){
        for(j=0;j<dim;j++){
            gsl_matrix_set(g_mus,i,j,p_mus[j+i*dim]);
        }
    }
    return;
}

void set_sig_inv(double *p_sigInv){
    int i,j;
    
    for(i=0; i<numMix*dim;i++){
        for(j=0;j<dim;j++){
          gsl_matrix_set(g_sigInv,i,j,p_sigInv[j+i*dim]);             
        }    
    }
    return;
}



//function that should be called to clean up c-allocated memory
void clean(){

    gsl_vector_free(g_gaussCoeffs);
    gsl_vector_free(g_mixPriors);
    
    gsl_matrix_free(g_mus);    
    gsl_matrix_free(g_sigInv);
    
    gsl_vector_free(x);
    gsl_vector_free(x_sqr);
    gsl_vector_free(g_vec);
    gsl_vector_free(gaussP);
    return;
}
