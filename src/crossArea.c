#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>

#include "crossArea.h"
#include "getPdbStructure.h"
#include "usage.h"

#define MIN_GUESSES           100   /*this threshold stops the code from exiting 
				                    **before the measure of estimated error is stable.*/ 
#define PI                    3.1415926535897932  
#define _SCAN_EXHAUSTIVE

/* local data types */
typedef struct atom_tag ATOM_T;
struct atom_tag {
  ATOM_T *next;
  double  gridX;
  double  gridY;
  int     id;
}; /* (recursively) define this tag so as to make atoms listable (as well as addressable in an array) */

/* declare local functions */
void fillLookupGrid(          ATOM_T      **lookupGrid, 
			      int           lookupGridLength, 
			      double       *RS_proj, 
			      ATOM_T       *atoms, 
			      PDB_STRUCT_T *PDB_struct,
			      double        gasRadius);
int  lookupGridTestCollision( double        x, 
			      double        y, 
			      ATOM_T      **lookupGrid, 
			      int           lookupGridLength, 
			      PDB_STRUCT_T *PDB_struct,
			      double        gasRadius);
int  testGridSquare( ATOM_T *headAtom, PDB_STRUCT_T *PDB_struct, double x, double y, double gasRadius );
void randomRotation(PDB_STRUCT_T *struc);

/* interface section - read, typecast and verify arguments then call the calculation */
int main(int argc, char **argv)
{
  int    opt, optIndex, seed, nSteps, nOrients;
  double gasRadius;
  char  *fileName, *radiusFileName, *logFileName;
 
  static struct option longOptions[] = {
    {"gasradius",  1, 0, 'g'},
    {"infile",     1, 0, 'i'},
    {"logfile",    1, 0, 'l'},
    {"nsteps",     1, 0, 'n'},
    {"norients",   1, 0, 'o'},
    {"radlib",     1, 0, 'r'},
    {"seed",       1, 0, 's'},
    {"verbose",    0, 0, 'v'},
    {'\0', 0, 0, '\0'},
  };

  /* initialise arguments */
  gasRadius       = -1.0;
  fileName        = NULL;
  radiusFileName  = NULL;
  logFileName     = NULL;
  seed            =  0;
  nSteps          =  0;
  nOrients        = 10; //default number of initial orientations
  
  crossArea_verboseFlag = 0;
  crossArea_logFile      = stdout; /* write to stdout by default */
  crossArea_errorLogFile = stderr;

  /* loop over arguments passed in */
  for ( opt = getopt_long( argc, argv, "g:i:l:n:o:r:s:v",
        longOptions, &optIndex );
        opt != (char)-1;
        opt  = getopt_long( argc, argv, "g:i:l:n:o:r:s:v",
        longOptions, &optIndex ) ) {
    
    switch ( opt ) {
    case 'g':
      gasRadius = atof( optarg );
      break;
    case 'i':
      fileName = (char *)malloc( ( strlen( optarg ) + 1 ) * sizeof (char) );
      strcpy( fileName, optarg );
      break;
    case 'l':
      logFileName = (char *)malloc( ( strlen( optarg ) + 1 ) * sizeof (char) );
      strcpy( logFileName, optarg );
      crossArea_logFile = fopen( logFileName,"w" );
      if( crossArea_logFile == NULL )
      {
	fprintf( crossArea_errorLogFile, "Could not open logfile, name: \'%s\'", logFileName );
	return( EXIT_FAILURE );
      }
      crossArea_errorLogFile = crossArea_logFile; /* write both kinds of output to same file */
      break;
    case 'n':
      nSteps    = atoi( optarg );
      break;
    case 'o':
      nOrients  = atoi( optarg );
      break;
    case 'r':
      radiusFileName = (char *)malloc( ( strlen( optarg ) + 1 ) * sizeof (char) );
      strcpy( radiusFileName, optarg );
      break;
    case 's':
      /* seed the random number generator to ensure repeatability */
      seed = atoi( optarg );
      srand( seed );
      break;
    case 'v':
      fprintf( crossArea_logFile,"Verbose output engaged!\n");
      crossArea_verboseFlag = 1;
      break;
    default:
      fprintf( crossArea_errorLogFile, "Unrecognized option %c\n", opt );
      fprintf( crossArea_errorLogFile, usage );
      fprintf( stderr, "Unrecognized option %c\n", opt );
      fprintf( stderr, usage );
      
      return( EXIT_FAILURE );
      break;
    }
  }
 
  /* make sure everything is present */
  if( nSteps <= 0 || gasRadius < 0.0 || fileName == NULL || radiusFileName == NULL )
  {     
     fprintf( crossArea_errorLogFile, usage );
     fprintf( stderr, usage );
     return( EXIT_FAILURE );
  }

  crossArea(  nOrients, nSteps, gasRadius, fileName, radiusFileName );

  return( EXIT_SUCCESS );
}


double crossArea( int nOrients, int nSteps, double gasRadius, char *fileName, char *radiusFilename )
{
  PDB_STRUCT_T *PDB_struct;
  int           angleCount, fibStep, orientStep;
  double        theta, phi;
  double        Rx[3][3], Ry_andProj[3][3], *RStheta, *RS_proj;
  double        point[2], pointFromOrigin;
  double        usePointsInside;
  double        p_hit, areaEstimate, stdDevEstimate;
  double        meanArea, meanError, errorRatio;
  int           atomIndex, guessIndex, crdIndex;
  int           hitCount, lookupGridLength;
  ATOM_T      **lookupGrid, *atomTags;
  double        goldenRatio;

  fprintf( crossArea_logFile,  "Reading coordinates from file \"%s\" and radii from \"%s\"\n", fileName, radiusFilename);

  PDB_struct = getPdbStructure( fileName, radiusFilename, 1 );

  /*adjust the size of the box to allow a probe atom to fit in around the edges */
  PDB_struct->L = PDB_struct->L + 2.0 * ( gasRadius + PDB_struct->largestAtomicRadius );

  /*let the user know whats happening */
  fprintf( crossArea_logFile, "\nStarting Monte-Carlo area measurement.\n" );

  /* allocate space for rotated coordinates */
  RStheta = (double*)malloc( 3 * PDB_struct->nAtoms * sizeof(double) );
  RS_proj = (double*)malloc( 3 * PDB_struct->nAtoms * sizeof(double) );

  /* allocate and init space for atom lookup grid */
  atomTags   = (ATOM_T *)malloc( PDB_struct->nAtoms * sizeof(ATOM_T) );
  for( atomIndex = 0; atomIndex < PDB_struct->nAtoms; atomIndex++ )
  {
    atomTags[atomIndex].next = (ATOM_T *)NULL;
    atomTags[atomIndex].id   = atomIndex;
  }
  lookupGridLength = 1 + ((int)(PDB_struct->L / ( gasRadius + PDB_struct->largestAtomicRadius ) ) );
  lookupGrid = (ATOM_T **)malloc( lookupGridLength * lookupGridLength * sizeof(ATOM_T *) );

  /* unchanging parts of the rotation matrix about the x axis */
  Rx[0][0] = 1.0;
  Rx[0][1] = 0.0;
  Rx[0][2] = 0.0;
  Rx[1][0] = 0.0;
  Rx[2][0] = 0.0;
  
  /* finding mean tumbling area.. */
  angleCount   = 0;
  meanArea     = 0.0;
  meanError    = 0.0;

  /** The most irrational of all numbers */
  goldenRatio  = (1.0 + sqrt(5.0)) / 2.0; 
  
  /* loop over randomly oriented grids */
  for( orientStep = 0; orientStep < nOrients; orientStep++ ){

    /* shift the coordinates to a random orientation */
    randomRotation(PDB_struct);

  /** Sub-Loop over orientations on a Fibbonacci grid.
   **
   ** For pseudocode (and helpful diagram)
   ** see in-press article by Gonzalez in Mathematical 
   ** Geosciences.
   ** 
   ** ArXiv: http://dx.doi.org/10.1007/s11004-009-9257-x
   **
   */
  for( fibStep = -nSteps; fibStep <= nSteps; fibStep++ ){

    //asin: takes values on (-1,1) and returns values on [-pi/2....pi/2].
    double  r;
    r     =  (2*fibStep) / (double)(2*nSteps+1);
    theta = asin( r );
    theta = theta + 0.5*PI;

    phi   = 2.0 * PI * (fibStep) / goldenRatio;
   
    /*printf("##angle: %f %f\n", theta, phi);*/

    /* rotation matrix around x-axis 
    ** Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta) ]; */
    Rx[1][1] =        cos(theta);
    Rx[1][2] = -1.0 * sin(theta);
    Rx[2][1] =        sin(theta);
    Rx[2][2] =        cos(theta);

    /* rotate all the coordinates about the x-axis */
    for( crdIndex = 0; crdIndex < PDB_struct->nAtoms * 3; crdIndex += 3 ){ 
      /*RStheta = PDB_struct.crds(:,1:3) * Rx;*/
      RStheta[crdIndex]     = PDB_struct->crds[crdIndex]; /* rotation is about x axis */
      RStheta[crdIndex + 1] = Rx[1][1] * PDB_struct->crds[crdIndex + 1] +
                              Rx[1][2] * PDB_struct->crds[crdIndex + 2];

      RStheta[crdIndex + 2] = Rx[2][1] * PDB_struct->crds[crdIndex + 1] + 
                             Rx[2][2] * PDB_struct->crds[crdIndex + 2];
    }

    /* rotate around y-axis and project into x-y plane */
    /* Ry_andProj = [ cos(phi) 0 sin(phi); 0 1 0; 0 0 0]; */
    Ry_andProj[0][0] = cos(phi);
    Ry_andProj[0][2] = sin(phi);
    Ry_andProj[1][1] = 1.0;
      
    for( crdIndex = 0; crdIndex < PDB_struct->nAtoms * 3; crdIndex += 3 )
    { 
       /*RStheta = PDB_struct.crds(:,1:3) * Rx;*/
       RS_proj[crdIndex]     = Ry_andProj[0][0] * RStheta[crdIndex] + Ry_andProj[0][2] * RStheta[crdIndex + 2];
       RS_proj[crdIndex + 1] = RStheta[crdIndex + 1];
       RS_proj[crdIndex + 2] = 0.0;
    }

    /* allocate each atom centre to a grid square for quick lookup */
    fillLookupGrid( lookupGrid, lookupGridLength, RS_proj, atomTags, PDB_struct, gasRadius );

    if( crossArea_verboseFlag == 1 ){
      fprintf( crossArea_logFile,  "Projecting at angles: %g %g\n", theta, phi );
    }

    /* start guessing if points are in or out of the shape
    ** ...scale max number of guesses with size of shape
    ** but do not expect to reach this maximum */
    hitCount = 0;
    for( guessIndex = 1; guessIndex <= ( PDB_struct->L * PDB_struct->L ); guessIndex++ )     
    {
           /* random x-y coordinates in square of side L */
           point[0] = ((rand()/((double)RAND_MAX)) - 0.5) * PDB_struct->L;
           point[1] = ((rand()/((double)RAND_MAX)) - 0.5) * PDB_struct->L;
      
                      
	       /*Test that the point is not in the corners of the square*/
           pointFromOrigin  =  ( point[0] * point[0] + point[1] * point[1] );
           usePointsInside  =   PDB_struct->L; 
           usePointsInside *=   usePointsInside;
	  
           if( pointFromOrigin <= usePointsInside  )
	       {

/* compiler switch, can use this to swap out the lookup grid for testing */
#ifndef NO_LUT

	    /* test for collision using the grid */ 
	    hitCount += lookupGridTestCollision( point[0], 
						   point[1], 
						   lookupGrid, 
						   lookupGridLength, 
						   PDB_struct, 
						   gasRadius );
#else
        /* test for collision exhaustively over atoms*/
	    crdIndex = -3; 
	    for( atomIndex = 0; atomIndex < PDB_struct->nAtoms; atomIndex++ ) 
	    {
	      double clearance, distance, vectorDistance[2];
	      
	      crdIndex += 3;
	      vectorDistance[0] = point[0] - RS_proj[crdIndex];
	      vectorDistance[1] = point[1] - RS_proj[crdIndex + 1];

	      clearance         = PDB_struct->atomicRadii[atomIndex] + gasRadius;

	      /*Test for a collision
	      ** sqrt is numerically expensive */
	      distance = vectorDistance[0]*vectorDistance[0] + vectorDistance[1]*vectorDistance[1]; 
	      if( distance <= ( clearance * clearance ) )
	      { 
             hitCount += 1;
		     break;
	      }
	    } 
	       
#endif
	   }/* close if */


       if( guessIndex % 100 == 0 ) /* test for convergence every X guesses */
	   {
          if( guessIndex >= MIN_GUESSES )
	      {
            p_hit          = hitCount / (double)guessIndex;  

            /* estimated std deviation of a binomial distribution */
            stdDevEstimate = sqrt( p_hit * ( 1 - p_hit ) / (double)guessIndex ); 
            errorRatio     = stdDevEstimate / p_hit;
                 
		    if( guessIndex % 1000 == 0 && crossArea_verboseFlag )
		    {
                   fprintf( crossArea_logFile, "iteration: %i estimated error ratio: %f estimated area: %f\n", guessIndex, errorRatio, pow( PDB_struct->L, 2 ) * p_hit );
            }

            if( errorRatio < 0.001 && p_hit != 0.0 && p_hit != 1.0 )
            {
		      break; /*stop making guesses, we have a reliable estimate*/
            }
          }
       }
       }      /*close loop over guesses */

     /*binomial distribution: expected number of hits  = guessIndex * p_hit
     **unbiased estimate of p_hit                      = hitCount / guessIndex;
     **
     **estimated variance( hitCount )                 ~= p_hit * ( 1 - p_hit )/guessIndex
     **so for 2DP accuracy, require sqrt(variance)/p_hit << 0.001
     */       

       p_hit           = hitCount / (double)guessIndex;
       areaEstimate    = pow( PDB_struct->L, 2 );
       areaEstimate   *= p_hit;

       stdDevEstimate  = pow( PDB_struct->L, 2 );
       stdDevEstimate *= sqrt( p_hit * ( 1 - p_hit ) / (double)guessIndex );

       meanArea   = meanArea + areaEstimate;
       meanError  = meanError + stdDevEstimate;
       angleCount++;

       if( crossArea_verboseFlag == 1 || crossArea_verboseFlag == 0 )
       {
	fprintf( crossArea_logFile, "Calculation converged for this set of angles. Area: %g Estimated Error: %g\n", areaEstimate, stdDevEstimate);
       }
   }  //close loop over grid orientations

   if( crossArea_verboseFlag == 1 || crossArea_verboseFlag == 0 ){
     fprintf( crossArea_logFile, "Finished grid loop at given grid orientation, running estimates: %f pm %f\n",
              meanArea/(double)angleCount, meanError/(angleCount*sqrt(angleCount - 1.0)));
   }

  }   //close loop over randomly oriented grids


  /* work out mean cross sectional area */
  meanArea  = meanArea  / (double)angleCount; 
  meanError = meanError / (double)angleCount;
  fprintf( crossArea_logFile, "Mean Area Over All Projections: %g\nTotal ESE: %g\n", meanArea, meanError / sqrt(angleCount - 1.0) );

  return( meanArea );

}

void fillLookupGrid( ATOM_T      **lookupGrid, 
		     int           lookupGridLength, 
		     double       *RS_proj, 
		     ATOM_T       *atoms, 
		     PDB_STRUCT_T *PDB_struct,
		     double        gasRadius )
{
  int     crdIndex, i, j, atomIndex, gridIndex;
  int     nAtoms;
  double  x, y, l; 

  nAtoms = PDB_struct->nAtoms;

  /* refresh the grid */
  memset( lookupGrid, 0, lookupGridLength * lookupGridLength * sizeof(ATOM_T *) );
 
  /* half the grid size */
  l = PDB_struct->L * 0.5;  

  atomIndex = 0;
  for(crdIndex = 0; crdIndex < nAtoms * 3; crdIndex += 3 )
  {

    /* get coords and convert them to a grid square */
    x = RS_proj[crdIndex];
    y = RS_proj[crdIndex + 1];
    i = (int)( ( x + l )/ ( gasRadius + PDB_struct->largestAtomicRadius ) );
    j = (int)( ( y + l )/ ( gasRadius + PDB_struct->largestAtomicRadius ) );

    gridIndex = i * lookupGridLength + j;
   
    if( (void *)lookupGrid[gridIndex] == NULL )
    {
      /* add the first atom to this grid square */
      lookupGrid[gridIndex]  = &atoms[atomIndex];
      atoms[atomIndex].next  = NULL;
      atoms[atomIndex].gridX = x;
      atoms[atomIndex].gridY = y;
    }
    else
    {
      /* insert a subsequent atom at the head of the list of atoms in this square */
      atoms[atomIndex].next  =  lookupGrid[gridIndex];
      lookupGrid[gridIndex]  = &atoms[atomIndex];
      atoms[atomIndex].gridX = x;
      atoms[atomIndex].gridY = y;

    }


    atomIndex++;

  }

  return;

}
int lookupGridTestCollision( double        x, 
			     double        y, 
			     ATOM_T      **lookupGrid, 
			     int           lookupGridLength, 
			     PDB_STRUCT_T *PDB_struct,
			     double        gasRadius)
{
  int     i, j, gridIndex;
  double  l;
  ATOM_T *headAtom;

  l = PDB_struct->L * 0.5;

  i = (int)( ( l + x ) / ( gasRadius + PDB_struct->largestAtomicRadius ) );
  j = (int)( ( l + y ) / ( gasRadius + PDB_struct->largestAtomicRadius ) );

  /* check atoms in this grid square untill a collision, or reached last atom */
  gridIndex = i * lookupGridLength + j;
  headAtom  = lookupGrid[ gridIndex ];
  if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
  {
    return( 1 );
  }

  /* to be exhaustive, also have to check the eight neighbouring grid squares */ 
  if( i > 0 )
  {
    if( j > 0 )
    {
      gridIndex = ( i - 1 ) * lookupGridLength + j - 1;
      headAtom  = lookupGrid[ gridIndex ];
      if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
      {
        return( 1 );
      }
    }

    gridIndex = ( i - 1 ) * lookupGridLength + j;
    headAtom  = lookupGrid[ gridIndex ];

    if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
    {
      return( 1 );
    }

    if( j < lookupGridLength - 1 )
    {
      gridIndex = ( i - 1 ) * lookupGridLength + j + 1;
      headAtom  = lookupGrid[ gridIndex ];
      if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
      {
        return( 1 );
      }
    }
  }

  if( j > 0 )
  {
    gridIndex = i * lookupGridLength + j - 1;
    headAtom  = lookupGrid[ gridIndex ];
    if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
    {
      return( 1 );
    }
  }
  if( j < lookupGridLength - 1 )
  {
    gridIndex = i * lookupGridLength + j + 1;
    headAtom  = lookupGrid[ gridIndex ];
    if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
    {
      return( 1 );
    }
  }

  if( i < lookupGridLength - 1 )
  {

    if( j > 0 )
    {
      gridIndex = ( i + 1 ) * lookupGridLength + j - 1;
      headAtom  = lookupGrid[ gridIndex ];
      if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
      {
	return( 1 );
      }
    }

    gridIndex = ( i + 1 ) * lookupGridLength + j;
    headAtom  = lookupGrid[ gridIndex ];
    if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
    {
      return( 1 );
    }

    if( j < lookupGridLength - 1 )
    {
      gridIndex = ( i + 1 ) * lookupGridLength + j + 1;
      headAtom  = lookupGrid[ gridIndex ];
      if( testGridSquare( headAtom, PDB_struct, x, y, gasRadius ) == 1 )
      {
        return( 1 );
      }
    }
  }

  return( 0 );

}


int testGridSquare( ATOM_T *gridAtom, PDB_STRUCT_T *PDB_struct, double x, double y, double gasRadius )
{
  double xx, yy, clearance;

  while( gridAtom != NULL )
  {
    xx = x - gridAtom->gridX;
    yy = y - gridAtom->gridY;

    xx = xx * xx;
    yy = yy * yy;

    clearance = gasRadius + PDB_struct->atomicRadii[gridAtom->id];
    
    /* test for collision */
    if( xx + yy < clearance * clearance )
    {
      return( 1 );
    }

    gridAtom = gridAtom->next;
  }

  return( 0 );

}

/*Random rotation to apply before placing a 
 *Fibonacci grid over the molecule.*/
void randomRotation(PDB_STRUCT_T *struc){
  
  double theta, phi;
  double xx, yy, zz, cTheta, sTheta, cPhi, sPhi;
  int    i;

  /*random numbers on [0,1]*/
  theta = rand()/(double)RAND_MAX;
  phi   = rand()/(double)RAND_MAX;

  /* random angles on unit sphere */
  theta    = acos(2*theta-1.);
  phi     *= 2*PI;

  /*fprintf(stderr, "orientations: %f %f\n", theta, phi);*/

  cTheta = cos(theta);
  sTheta = sin(theta);
  cPhi   = cos(phi);
  sPhi   = sin(phi);
  
  /* reinitialise the coordinates */
  memcpy(struc->crds, struc->crds_store, 3*sizeof(double)*struc->nAtoms);

  /* reorient them */
  for(i = 0; i < struc->nAtoms; i++ ){
    
    /*apply longitude rotation*/
    xx = struc->crds[i*3+2]*sTheta + struc->crds[i*3]*cTheta;
    yy = struc->crds[i*3+1]; /* rot around y*/
    zz = struc->crds[i*3+2]*cTheta - struc->crds[i*3]*sTheta;

    /*apply latitude rotation*/
    struc->crds[i*3]   = yy*sPhi + xx*cTheta;
    struc->crds[i*3+1] = yy*cPhi - xx*sTheta;
    struc->crds[i*3+2] = zz; /* rot around z*/
  }

}
