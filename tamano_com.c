#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

unsigned int get_new_community_size(const unsigned int num_comunities){
	float cc=(float)num_comunities+1.0; //
	float aa= 0.800000;//  2.0003158  0.3999369 6.985290e-01
	float bb = 2.254545;//  0.3152941  7.1506099 5.362819e-05
	float max_number= num_comunities; 
	float minimum_number= 1.0;
	float myrand =  ((float)rand()/(float)(RAND_MAX));
	float num =  myrand * (max_number + 1 - minimum_number) + minimum_number;
	//int num= ( rand() % (max_number + 1 - minimum_number)) + minimum_number;   //queremos que sea un numero entre 1 y num_comunities
	
    float newsize =ceil( -bb/log((num-aa)/cc));
	if (newsize==1) newsize=2;
	printf(" num_comunities is %i random is %f newsize is %f\n",num_comunities, num,newsize);
	return (unsigned int) newsize;
}

//  -2.254545 / ((log (0.2+ rand() % num_comunities)/(num_comunities+1))
