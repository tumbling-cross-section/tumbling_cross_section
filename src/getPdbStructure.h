
/* useful definitions */
#define PDB_RECORD_LENGTH 240 /*much longer than standard record but many pdb file variants exist */

#define PDB_COLUMN_MAX_WIDTH    8
#define RADIUS_COLUMN_MAX_WIDTH 16

#define PDB_ATOM_NAME_WIDTH     4
#define PDB_RESIDUE_NAME_WIDTH  3
#define PDB_COORD_WIDTH         8

#define PDB_ATOM_NAME_START     12
#define PDB_RESIDUE_NAME_START  17
#define PDB_X_START             30
#define PDB_Y_START             38
#define PDB_Z_START             46


/*have a default value just in case*/
#define DEFAULT_LARGEST_ATOMIC_RADIUS 1.0   

typedef struct pdb_struct_tag {
  double      *crds;
  double      *atomicRadii;
  double       L;
  double       largestAtomicRadius;
  unsigned int nAtoms;
} PDB_STRUCT_T;

typedef struct radius_lookup_tag {
  char         atomId[RADIUS_COLUMN_MAX_WIDTH + 1];
  char         resId[RADIUS_COLUMN_MAX_WIDTH + 1];
  double       radius;
} RADIUS_LOOKUP_T;

typedef enum record_name_id_tag {
  RECORD_OTHER = 0,
  RECORD_ATOM,
  RECORD_ENDMDL
} RECORD_NAME_ID_T;


typedef enum regexp_return_tag {
  RRT_NO_MATCH = 0,
  RRT_PARTIAL_MATCH,
  RRT_EXACT_MATCH
} REGEXP_RETURN_T;





PDB_STRUCT_T *getPdbStructure( char *filename, char *radiusFilename, int modelNumber );

void  printPDB_struct( PDB_STRUCT_T *PDB_struct, FILE *outHandle   );
