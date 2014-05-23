#include <stdio.h>
#define MATHLIB_STANDALONE 
#include "Rmath.h"

main(){
   double shape1,shape2,prob;
   
   printf("Enter first shape parameter: ");
   scanf("%lf",&shape1);

   printf("Enter second shape parameter: ");
   scanf("%lf",&shape2);

   printf("Enter probability level: ");
   scanf("%lf",&prob);

   printf("Critical value is %lf\n",qbeta(prob,shape1,shape2,1,0));
}

