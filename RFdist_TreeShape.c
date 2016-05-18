///////////////////////////////////////////////
//
//  Created by Yujin Chung
//  Last modified date: July 14, 2012
//
//////////////////////////////////////////////


//////////////////////////////////////////////
// 
//
////////////////////////////////////////////



#include <sys/resource.h> 
#include <sys/time.h> 
#include <unistd.h> 
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> // malloc
#include "shape.h"
#include "misc.h"


int main()
{
  
  int i,j,k; // indexes for loops


  /*------------- Enter the number of taxa and the number of segments -------*/
  int no_taxa = 0;
  while(no_taxa <= 3)
  {
  	printf("Enter the number of taxa (>= 4): ");
	scanf("%d", &no_taxa);
  }  
  int flag_jointDist_givenTree = -1;
  while(flag_jointDist_givenTree != 0 &&flag_jointDist_givenTree != 1)
  {
  	printf("\nDo you want to compute the joint distribution for a given tree shape? (1:yes, 0: no): ");
  	scanf("%d", &flag_jointDist_givenTree);
  }
  int givenTreeShapeID = 0;
  int nshape = no_urshapes(no_taxa);
  if(flag_jointDist_givenTree == 1)
  {
     printf("\tThere are %d shapes.\n",nshape);
     while(givenTreeShapeID < 1 || givenTreeShapeID > nshape)
     {
     	printf("\tEnter the given tree shape ID (starting from 1): ");
     	scanf("%d", &givenTreeShapeID);
     }
  }

  int flag_fullJointDist = -1;
  while(flag_fullJointDist !=0 && flag_fullJointDist !=1)
  {
  	printf("\nDo you want to compute the joint distributions for any given tree shape? (1:yes, 0:no): ");
  	scanf("%d", &flag_fullJointDist);
  }

  int flag_partitionFt = -1;
  while(flag_partitionFt !=0 && flag_partitionFt !=1)
  {
  	printf("\nDo you want to compute the partition function of a Gibbs distribution? (1:yes, 0:no): ");
  	scanf("%d", &flag_partitionFt);
  }
  int nseg = 0;
  int flag_haveFullJointDist = -1;
  if(flag_partitionFt == 1)
  {
	while(nseg <=1)
 	{
  		printf("\tEnter the number of segments (>=2): ");
		scanf("%d",&nseg);
 	}
	
	if(flag_fullJointDist == 0)
	{
		while(flag_haveFullJointDist !=0 && flag_haveFullJointDist != 1)
		{
			printf("\tDo you have a file of the full joint distribution? (1:yes, 0:no): ");
			scanf("%d", &flag_haveFullJointDist);
		}
	}
  }


   printf("\n\n");
  
  /*------- computing the number of tree shapes------------*/
//  printf("\nComputing the number of tree shapes.....\n");
 // int nshape = no_urshapes(no_taxa);  
 // printf("\tthere are %d shapes. Done.\n\n",nshape);


  /*-------- calculating the joint distribution for a given tree shape (shape ID) ------------*/
  // int ID = 1; // shape ID
  long double** propDistance_2trees_each; 
  if(flag_jointDist_givenTree == 1)
  {
  	propDistance_2trees_each = malloc(nshape * sizeof(long double*));
  	for(j=0;j<nshape;j++)
    	{	    
      		propDistance_2trees_each[j] = malloc((no_taxa-2) * sizeof(long double));
    	}
  	printf("Computing the joint distribution of RF distance and a tree when a tree shape with ID=%d is given....\n",givenTreeShapeID);
  	eachfullJointDistribution(no_taxa,nshape,propDistance_2trees_each,1,givenTreeShapeID);
  	printf("Done.\n");
   }


  /*-------- calculating the full joint distribution ------------*/
  long double***  propDistance_2trees_full;
  if(flag_fullJointDist == 1 || flag_haveFullJointDist == 0)
  {
  	propDistance_2trees_full =  malloc(nshape * sizeof(long double**));
  	for(i=0;i<nshape;i++)
    	{
      		propDistance_2trees_full[i] = malloc(nshape * sizeof(long double*));
      		for(j=0;j<nshape;j++)
		{	     
	  		propDistance_2trees_full[i][j] = malloc((no_taxa-2) * sizeof(long double));  
		}
    	}
  	printf("Computing the full joint distribution of RF distance and a tree....\n");
  	print_fullJointDistribution(no_taxa,nshape,propDistance_2trees_full);
  	printf("Done.\n");
  }

  /*---------- Read joint distribution from a text file --------------*/
  char filename[10000];
  if(flag_fullJointDist == 0 && flag_haveFullJointDist == 1)
  {
  	printf("Reading jointDistribution....\n");
  	FILE *fopen(), *read_jointDist;
  	// char filename[10000];
  	sprintf(filename,"JointDistribution_%dtaxa.txt",no_taxa);
  	read_jointDist = fopen(filename,"r"); 
  	char buffer[100];
  	fgets(buffer, 100, read_jointDist);  fgets(buffer, 100, read_jointDist);  fgets(buffer, 100, read_jointDist); // skip the first two lines
  	char c;
  	// long double***  propDistance_2trees;
  	propDistance_2trees_full =  malloc(nshape * sizeof(long double**));
  	for(i=0;i<nshape;i++)
    	{
      		propDistance_2trees_full[i] = malloc(nshape * sizeof(long double*));
      		for(j=0;j<nshape;j++)
		{	    
	  		propDistance_2trees_full[i][j] = malloc((no_taxa-2) * sizeof(long double));
	  		// skip the first two columns
	  		do
	    			c = fgetc(read_jointDist);
	  		while (c != '\t');
	  			do
	    				c = fgetc(read_jointDist);
	  		while (c != '\t');
	  			for(k=no_taxa-3;k>=0;k--)
	    			{
	      				fscanf(read_jointDist, "%Le", &propDistance_2trees_full[i][j][k]);
	    			}					
		}
      		//printf("i=%d j=%d k=%d\n",i,j,k);
    	}
  	fclose(read_jointDist);
  	printf("Done.\n");
  }

  /*---------- Computing the normalizing constant and approximations --------------*/
  if(flag_partitionFt == 1)
  { 
  	long double var_normal, mean_L,var_L;
  	long double logZ, logZ_indep, logZ_Normalapprox;
  	FILE *out_file;
  	sprintf(filename,"partitionFnct_vs_beta_%dtaxa_%dsegments.txt",no_taxa,nseg);
  	out_file = fopen(filename,"w");
  	fputs("ntaxa\tnseg\tbeta\t",out_file);
  	//fputs("mean_L\tvar_L\t",out_file);
  	fputs("logZ\tlogZ_indep\tlogZ_NormalApprox\n",out_file);
  	fclose(out_file); 
  	printf("Computing the normalizing constant and approximations...");
  	int maxBeta = 15;
  	long double delta = 0.05;
  	int noBeta = (maxBeta - .001/maxBeta)/delta +1; 
  	long double Beta[noBeta],beta; 
  	int NoSeg[3]={nseg,50,100},n,b;
  	for(i=0;i<noBeta;i++)
    		Beta[i]=.00005/maxBeta+i*delta;
  	for(n=0;n<1;n++)
    	{
      		nseg=NoSeg[n]; 
      		b=0;
      		for(b=0;b<noBeta;b++)
		{
	  		beta=Beta[b];
	  		//printf("b=%Le. start computing the true logZ\n",beta);
	  		logZ= logZ_L(no_taxa, nshape, nseg, beta, propDistance_2trees_full);
	  		//print_cpu_time();
	  		//printf("b=%Le. start computing the indep approx logZ\n",beta);
	  		logZ_indep = log_normalizingConstant_indep(propDistance_2trees_full,no_taxa, nseg, nshape, beta);
	  		//print_cpu_time();
	  		//printf("b=%Le. start computing the normal approx logZ\n",beta);
	  		logZ_Normalapprox =log_normalizingConstant_Normalapprox(nseg,no_taxa, nshape, beta, propDistance_2trees_full);
	  		//print_cpu_time();
	  		mean_L  = (nseg - 1)*expected_dist(no_taxa, nshape,propDistance_2trees_full);
	  		var_L = (nseg-1) * (var_dist(no_taxa, nshape,propDistance_2trees_full)+cov_dist(no_taxa, nshape,propDistance_2trees_full));
	  
	  		out_file = fopen(filename,"a");
	  		fprintf(out_file,"%d\t%d\t%Le\t",no_taxa,nseg,beta);
	  		// fprintf(out_file,"%Le\t%Le\t",mean_L,var_L);
	  		fprintf(out_file,"%Le\t%Le\t%Le\n",logZ,logZ_indep,logZ_Normalapprox);
	  		fclose(out_file); 
		}		 
    	}
  	printf("DONE\n");
  	print_cpu_time();
  }


  return 0;
}

