
/*********************************************
 ** Copyright (C) 2008 University of Leeds.
**
** Authors: Josh Berryman and Tom Knapman <chm2twk@leeds.ac.uk>
**
*/

/* bring in a few standard utility functions */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

/* bring in interface definitions */
#include "getPdbStructure.h"
#include "crossArea.h"

/* define local functions */
int              countAtoms (FILE *file );
RECORD_NAME_ID_T getRecordType( char *inLine );
unsigned int     removeBlanks( char *in, char *out, unsigned int inFieldLength );
int              compareSubstr( char *shortString, char *longString );
RADIUS_LOOKUP_T *readRadii( FILE *radiusFile, int *count  );
REGEXP_RETURN_T  matches( char *mask, char *nameString );
double           assignRadius( RADIUS_LOOKUP_T *radii, char *atomName, char *resName, int numDefs );
void             centreCoords( PDB_STRUCT_T *PDB_struct );
double           setLargestRadius( PDB_STRUCT_T *PDB_struct  );

/****************************************************%
** This function reads in a pdb file and a file of atomic radius definitions
** and returns a structure containing the :
** coordinates, radii, number of atoms and longest axis length of the structure
**
*/

PDB_STRUCT_T *getPdbStructure( char *filename, char *radiusFilename, int modelNumber )
{


  RECORD_NAME_ID_T recordName;
  PDB_STRUCT_T    *PDB_struct;
  FILE            *radius_fid, *fid;
  double           radius, maxRadius, radiusFromCentre;
  int              numOfAtom, missedAtoms, readingModel, count, numDefs;
  int              crdIndex;
  char             tline[PDB_RECORD_LENGTH], atType[PDB_COLUMN_MAX_WIDTH], resType[PDB_COLUMN_MAX_WIDTH];
  char             crdStr[PDB_COLUMN_MAX_WIDTH];
  RADIUS_LOOKUP_T *radii;

  fid         = fopen( filename, "r" );
  radius_fid  = fopen( radiusFilename,"r" );


  if ( fid == NULL ) 
  {
   fprintf( crossArea_errorLogFile, "Unable to open specified pdb file : %s\n", filename );
    return( (void *)NULL );
  }
  if ( radius_fid == NULL ) 
  {
   fprintf( crossArea_errorLogFile, "Unable to open specified radiues file : %s\n", radiusFilename );
    return( (void *)NULL );
  }
			
	
  radii = readRadii( radius_fid, &numDefs );
  fprintf( crossArea_logFile, "read %i radius definitions from file %s\n", 
	   numDefs, radiusFilename );


  count = countAtoms( fid );

  /*Allocate and Initialise*/
  PDB_struct               = (PDB_STRUCT_T *)malloc( sizeof( PDB_STRUCT_T ) );
  PDB_struct->crds         = (double*)malloc( 3 * count * sizeof( double) );
  PDB_struct->atomicRadii  = (double*)malloc( count * sizeof( double) );
  PDB_struct->nAtoms       =  (unsigned int)count; 
  numOfAtom   = 0;
  crdIndex    = 0;
  maxRadius   = 0.0;
  missedAtoms = 0;

  if( crossArea_verboseFlag == 1 )
  {
    fprintf( crossArea_logFile, "Allocated memory for up to %i atoms\n", count);
    fprintf( crossArea_logFile, "Reading input pdb file and assigning radii to atoms by name...\n");
  }

  readingModel = 1;
  while ( readingModel <= modelNumber )
  {

    if( !fgets( tline, PDB_RECORD_LENGTH, fid ) )
    {
      break; /* end of file */
    }

    if ( strlen( tline ) > 6 )     /*skip lines too short to have a valid record name.*/
    {         
        recordName = getRecordType( tline ); /* record name is an enumerator type */
    		
	switch ( recordName )
	{
	    case RECORD_ATOM:  
	       if( readingModel == modelNumber )
	       {
		  /*read the atom type so we can then deduce the radius */
		  removeBlanks( &tline[PDB_ATOM_NAME_START], atType, PDB_ATOM_NAME_WIDTH );                     

		  /*radius is residue-type dependent also*/
		  removeBlanks( &tline[PDB_RESIDUE_NAME_START], resType, PDB_RESIDUE_NAME_WIDTH );

		  /* link the atom names and radii */
		  radius = assignRadius( radii, atType, resType, numDefs );	

                  /* radius 0.0 signifies none assigned because hydrogen */
		  if( radius < 0.0 )
		  {    
		      missedAtoms = missedAtoms + 1;
		      fprintf( crossArea_errorLogFile, "\nDid not assign radius for atom: %s %s \n", resType, atType );
                      fprintf( crossArea_errorLogFile, "status flag %g\n", radius);
		      return( 0 );
		  }
		  else if( radius > 0.0 )	
		  {

              /**record atomic radius*/                    
	   	      PDB_struct->atomicRadii[numOfAtom] = radius;
		    
              /**record the coords*/								   
              if( removeBlanks( &tline[PDB_X_START], crdStr, PDB_COORD_WIDTH  ) == 0 )
		      {
			        fprintf( crossArea_errorLogFile, "trouble reading coordinates, check formatting\n" );
			        fprintf( crossArea_errorLogFile, "line:\n%s\n", tline );
			        return( 0 );
              }
		      PDB_struct->crds[crdIndex] = atof( crdStr );


              if( removeBlanks( &tline[PDB_Y_START], crdStr, PDB_COORD_WIDTH ) == 0 )
		      {
			        fprintf( crossArea_errorLogFile, "trouble reading coordinates, check formatting\n" );
			        fprintf( crossArea_errorLogFile, "line:\n%s\n", tline );
			        return( 0 );
              }
		      PDB_struct->crds[crdIndex + 1] = atof( crdStr );		
                    

              if( removeBlanks( &tline[PDB_Z_START], crdStr, PDB_COORD_WIDTH ) == 0 )
		      {
			        fprintf( crossArea_errorLogFile, "trouble reading coordinates, check formatting\n" );
			        fprintf( crossArea_errorLogFile, "line:\n%s\n", tline );
			        return( 0 );
              }
              PDB_struct->crds[crdIndex + 2] = atof( crdStr );

                      /* record squared distance from centre */
		      radiusFromCentre  = PDB_struct->crds[crdIndex] * PDB_struct->crds[crdIndex];
		      radiusFromCentre += PDB_struct->crds[crdIndex + 1] * PDB_struct->crds[crdIndex + 1];
		      radiusFromCentre += PDB_struct->crds[crdIndex + 2] * PDB_struct->crds[crdIndex + 2];

		      if( radiusFromCentre >= maxRadius )
		      {
		        maxRadius = radiusFromCentre;
		      }

		      /**increment the atom and coordinate counters*/
		      numOfAtom++;
                      crdIndex += 3;

		      /**Display a minimal ptraj-style progress indicator**/
		      if( numOfAtom % 10  == 0 )
		      {   
			  fprintf( crossArea_logFile, " . " );
			  if( numOfAtom % 100 == 0 )
			  {
			    fprintf( crossArea_logFile, "\n" );
		     	  }
		       }
		    }
	         }
		 break;								       
	         case RECORD_ENDMDL:
		   readingModel = readingModel + 1; /** we have the next structure in an NMR pdb file. */
	         default:
		 break;
	} /* end of the SWITCH statement */ 
    } /* end of the IF statement checking the empty string */
    else {
      fprintf( crossArea_logFile, "Skipping line, too short: '%s'\n", tline );
      fprintf( stderr,            "Skipping line, too short: '%s'\n", tline );
    }

} /*end of the WHILE loop */

		    /*report missed atoms and offer advice */
if( missedAtoms > 0 )
{
  fprintf( crossArea_logFile, "\n%i atoms in the input did not have a radius defined in %s.\n", missedAtoms, radiusFilename );
  fprintf( crossArea_logFile, "Unless they were all hydrogens and united-atom radii are being used, this should probably be fixed.\n");
}
else
{
  fprintf( crossArea_logFile, "\nAll atoms were assigned radii successfully.\n" );
}

 PDB_struct->L              =  sqrt( maxRadius ); /**we have been working in squared-distances up to here. */
 PDB_struct->L              =  2 * PDB_struct->L; 
 PDB_struct->nAtoms         =  numOfAtom;         /* correcting the initial over-estimate */

 setLargestRadius( PDB_struct );

 centreCoords( PDB_struct );

 if( crossArea_verboseFlag == 1 ) 
 {
   printPDB_struct( PDB_struct, crossArea_logFile );
 }

 return( PDB_struct );

}

void centreCoords( PDB_STRUCT_T *PDB_struct )
{
  int    crdIndex;
  double centre[3];
  
  centre[0]=0.0;
  centre[1]=0.0;
  centre[2]=0.0;
  
  for( crdIndex = 0; crdIndex < PDB_struct->nAtoms * 3; crdIndex += 3 )
  {
    centre[0] += PDB_struct->crds[crdIndex];
    centre[1] += PDB_struct->crds[crdIndex + 1];
    centre[2] += PDB_struct->crds[crdIndex + 2];
  }

  centre[0] /= ((double)PDB_struct->nAtoms);
  centre[1] /= ((double)PDB_struct->nAtoms);
  centre[2] /= ((double)PDB_struct->nAtoms);

  for( crdIndex = 0; crdIndex < PDB_struct->nAtoms * 3; crdIndex += 3 )
  {
    PDB_struct->crds[crdIndex]     -= centre[0];
    PDB_struct->crds[crdIndex + 1] -= centre[1];
    PDB_struct->crds[crdIndex + 2] -= centre[2];
  }

  return;
}

/* function to print out pdb structure, for debug/test */
void  printPDB_struct( PDB_STRUCT_T *PDB_struct, FILE *outHandle   )
{
  int i, crdIndex;

  fprintf( outHandle, "L: %f\n", PDB_struct->L );
  fprintf( outHandle, "nAtoms: %i\n", PDB_struct->nAtoms );
  fprintf( outHandle, "largest atom radius: %f Angstrom\n",
	                         PDB_struct->largestAtomicRadius );
    
  crdIndex = 0;
  for( i = 0; i < PDB_struct->nAtoms; i++ )
  {
    fprintf( outHandle, "%i %3.3f %3.3f %3.3f  %3.3f\n", i, 
	     PDB_struct->crds[crdIndex],
	     PDB_struct->crds[crdIndex + 1],
	     PDB_struct->crds[crdIndex + 2], 
	     PDB_struct->atomicRadii[i] );

    crdIndex = crdIndex + 3;
  }


}

double setLargestRadius( PDB_STRUCT_T *PDB_struct  ){

  double r;
  int    i;

  r = DEFAULT_LARGEST_ATOMIC_RADIUS;

  for( i = 0; i < PDB_struct->nAtoms; i++ ) {
    if( PDB_struct->atomicRadii[i] > r ) {
      r = PDB_struct->atomicRadii[i];
    }
  }
  
  PDB_struct->largestAtomicRadius = r;

  return( r );

}




/**************locally scoped functions: */



/* REMOVEBLANKS removes spaces from either end of a string, writes to 'out'
** returns chars written not including terminator: ie buffer size if overrun */
unsigned int removeBlanks( char *in, char *out, unsigned int inFieldLength )
{
   unsigned int count = 0, writing = 0;

   while( in[count] != '\0' ) /*loop until string terminator*/
   {
     if( in[count] != ' ' )
     {
         out[writing++] = in[count];
     }
     else if ( writing > 0 )
     {
         out[writing] = '\0';
         return( writing ); /* return on first space after writing something */
     }
     count++;

     if( count == inFieldLength ) {
        out[writing] = '\0';
        return( writing ); /* return because have finished the input field */
     }
   }

   /* return on end of input string */
   out[writing] = '\0';
   return( writing );
}


int countAtoms (FILE *file )
{
  char  tline[PDB_RECORD_LENGTH];
  int   count;

  /* how many lines in the file? just get an upper bound, don't really expect memory to be a problem here.*/
  count = 0;
  while( NULL != fgets( tline, PDB_RECORD_LENGTH - 1, file ) )
  {
    
    count++;

  }
  fseek( file, 0, SEEK_SET ); /* back to beginning */

  return( count );
}

RADIUS_LOOKUP_T *readRadii( FILE *radiusFile, int *pt_count  )
{
  char             tline[88], radiusString[PDB_COLUMN_MAX_WIDTH];
  int              count; 
  RADIUS_LOOKUP_T *radii;
  
  /* how many lines in the file? */
  count = 0;
  while( NULL != fgets( tline, 88, radiusFile ) )
  {
    count++;
  }
  fseek( radiusFile, 0, SEEK_SET ); /* back to beginning */

  /* allocate an array of atom-res-radius triplets */
  radii = (RADIUS_LOOKUP_T *)malloc( count * sizeof(RADIUS_LOOKUP_T) );

  /* now read through properly */
  count = 0;
  while( NULL != fgets( tline, 88, radiusFile ) )
  {
       /*identify radius definition */
      removeBlanks( tline, radii[count].resId, RADIUS_COLUMN_MAX_WIDTH );
      removeBlanks(&tline[4], radii[count].atomId, RADIUS_COLUMN_MAX_WIDTH );
      removeBlanks(&tline[8], radiusString, RADIUS_COLUMN_MAX_WIDTH );
            
      radii[count].radius = atof( (const char*)radiusString );/* convert string to double*/
    
      fprintf(stderr, "%s %s\n%s\n%s\n", tline,radii[count].resId ,radii[count].atomId, radiusString );

      count++;
   }

  *pt_count = count;

  return( radii );

}


/**%function to link atomname and residue type to its radius*/
double assignRadius( RADIUS_LOOKUP_T *radii, char *atomName, char *resName, int numDefs )
{
  int    count, resMatch, atomMatch;
  double radius  = -1.0; /** 'fail' value by default*/

  for( count = 0; count < numDefs; count++ )
  { 

    resMatch  = matches( radii[count].resId, resName );
    if( resMatch == 0 ) /**no joy here, move along the loop */
    {
      continue;
    }

    atomMatch = matches( radii[count].atomId, atomName );
    if( atomMatch == 0 ) /**no joy here, move along the loop */
    {
      continue;
    }
   
    radius = radii[count].radius; /**have matched a pattern but could still achieve exact match */
    if( atomMatch + resMatch == 4 )   /**have achieved exact match - job done */
    {
       return( radius );
    }

  }

  return( radius );

}

   /**%function to compare strings*/
REGEXP_RETURN_T matches( char *mask, char *nameString )
{
  int numChars, count;
  REGEXP_RETURN_T matchStatus;

  numChars = strlen( nameString );

  if( numChars != strlen( mask ) )
  {
     return( RRT_NO_MATCH );   /*mask is longer or shorter than target string therefore no match possible*/
  }

  matchStatus = RRT_EXACT_MATCH;

  for( count = 0; count < numChars; count++ )
  {
    if( mask[count] != nameString[count] )
    {
        if( mask[count] == '*' ) 
	{
            matchStatus = RRT_PARTIAL_MATCH;
        }
	else
	{
            return( RRT_NO_MATCH );
	}
    }
  }

  return( matchStatus );

}

RECORD_NAME_ID_T getRecordType( char *inLine )
{ 
   if( compareSubstr( "ATOM", inLine ) == 1 )
   {
        return( RECORD_ATOM );
   }
   if( compareSubstr( "HETATM", inLine ) == 1 )
   {
        return( RECORD_ATOM );
   }
   if( compareSubstr( "ENDMDL", inLine ) == 1 )
   {
        return( RECORD_ENDMDL );
   }
   return( RECORD_OTHER );
}

int compareSubstr ( char *shortStr, char* longStr )
{
  int i = 0;
   
  while( shortStr[i] != '\0'  )
  {
    if( shortStr[i] != longStr[i] )
    {
      return( 0 );
    }

    i++;  
  }

  return( 1 );

}



