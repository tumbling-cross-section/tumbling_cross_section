
function PDB_struct = getPdbStructure( filename, radiusFilename, modelNumber )

%% Copyright (C) 2008 University of Leeds.
%%
%% This software is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or (at
%% your option) any later version.
%%
%% This software is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%% 
%% For the licence see
%% <http://www.gnu.org/licenses/>.

%% Authors: Josh Berryman and Tom Knapman <chm2twk@leeds.ac.uk>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function reads in a pdb file and a file of atomic radius definitions
%% and returns a structure containing the :
%% coordinates, radii, number of atoms and longest axis length of the structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radius_fid = -1;
fid = -1;

%% check arguments and open infile
switch( nargin )
     case 0
       disp( "getPdbStructure:\n requires pdb filename and radius filename,\n optionally accepts number of model in pdb file."  );
     case 1
       disp( "getPdbStructure:\n requires pdb filename and radius filename,\n optionally accepts number of model in pdb file."  );
     case 2
       fid         = fopen(filename,'r');
       radius_fid  = fopen( radiusFilename,'r');
       modelNumber = 1;
     case 3
       fid        = fopen(filename,'r');
       radius_fid = fopen( radiusFilename,'r');
     otherwise
       error( 'ERROR! getPdbStructure called with incorrect arguments' );
       return
end

if ( fid == -1 ) 
    error( strcat('Unable to open specified pdb file : ', filename ) );
    return
end
if ( radius_fid == -1 ) 
    error( strcat('Unable to open specified radius file : ', radiusFilename ) );
    return
end
				
radii = readRadii( radius_fid   );
fprintf( stdout, "read %i radius definitions from file %s\n", size(radii)(2), radiusFilename );


%%Initialise
PDB_struct = [];
NumOfATOM = 0;
maxRadius = 0;

fprintf(stdout, "Reading input pdb file and assigning radii to atoms by name...\n");
missedAtoms = 0;

readingModel = 1;
while ( readingModel <= modelNumber )

    tline = fgetl(fid);
    if ~ischar(tline)
        break; % For end of file recognition
    end

    sz = size(tline);
    if ( sz(1) > 0 )     %%Omit the empty lines to avoid error of invalid matrix index.
             
            Record_name = tline(1:6);
            Record_name = deblank(Record_name); % Assuming that the record name will be left-aligned 
    		
	    switch ( Record_name )
                case 'ATOM'
    		    
		    if( readingModel == modelNumber )									       

		      %%read the atom type so we can then deduce the radius
                      atType =  removeblanks( tline(12:16) );                     

		      %%radius is residue-type dependent also
                      resType = removeblanks( tline(18:20) );                      			

                      radius = assignRadius( radii, atType, resType );	

                      if( radius < 0.0 )
                           
                           missedAtoms = missedAtoms + 1;
			   fprintf( stdout, "\nCould not assign radius for atom: %s %s \n", resType, atType );
		      
                      elseif( radius > 0.0 )	

		         %%increment the count the number of atoms and pick up coordinates and radii
                         NumOfATOM = NumOfATOM + 1;
			
                         %%record the coords								   
                         crds(NumOfATOM,:) = [ str2num(tline(31:38)), str2num(tline(39:46)), str2num(tline(47:54))];

                         atomicRadius(NumOfATOM) = radius;
		         radiusFromCentre        = crds(NumOfATOM,1)**2 +crds(NumOfATOM,2)**2 + crds(NumOfATOM,3)**2;

		         if( radiusFromCentre >= maxRadius )
		            maxRadius = radiusFromCentre;
                         end

			 %%Display a minimal ptraj-style progress indicator
			 if( mod( NumOfATOM, 10 ) == 0 )
                            fprintf( stdout, " . " );
			    if( mod( NumOfATOM, 100 ) == 0 )
                               fprintf( stdout, "\n" );
                            end
                         end

                      end

		    end
										       
                 case 'ENDMDL'
		    readingModel = readingModel + 1; %% we have the next structure in an NMR pdb file.

            end % of the SWITCH statement 

     end % of the IF statement checking the empty string

end % of the WHILE loop

%%report missed atoms and offer advice
if( missedAtoms > 0 )
  fprintf( stdout, "\n%i atoms in the input did not have a radius defined in %s.\n", missedAtoms, radiusFilename );
  fprintf( stdout, "Unless they were all hydrogens and united-atom radii are being used, this should probably be fixed.\n");
else
  fprintf( stdout, "\nAll atoms were assigned radii successfully.\n" );
end



%%Organise output structure
PDB_struct.nAtoms         =  NumOfATOM;	 											       
PDB_struct.crds           =  crds;   	 	 											       PDB_struct.atomicRadii    =  atomicRadius;  
PDB_struct.L              =  sqrt( maxRadius ); %%we have been working in squared-distances up to here.
PDB_struct.L              =  PDB_struct.L * 2;


return;

%%%%%%%%%%%%%%End of main pdb-reading function



% REMOVEBLANKS removes spaces
function out = removeblanks( in )

out = deblank( in ); %%remove trailing blanks

%%remove leading blanks
for count = 1:length(out)
     if( out(count) != " " )
       out = out( count:end );     
       return;
     end
end



function radii = readRadii ( radiusFile  )

tline = fgetl( radiusFile );
count = 1;
while( ischar(tline) )
    %%identify radius definition
    resid  = removeblanks(tline(1:3));
    atomid = removeblanks(tline(5:9));
    radius = str2double( removeblanks(tline(11:end)) );
    radii(:,count) = {resid, atomid, radius}';
    
    tline = fgetl( radiusFile );
    count = count + 1;
end


%%%function to link atomname and residue type to its radius
function radius = assignRadius( radii, atomName, resName )

count   = 0;
numDefs = length(radii);
radius  = -1.0; %% 'fail' value by default

%%define these return codes for our 'regexp parser'
%NO_MATCH      = 0;
%GENERAL_MATCH = 1;
%EXACT_MATCH   = 2;

%%crudest possible search over defined atom names
%%this should still be fast relative to the real calculation though, so not worried
while ( count < numDefs )

   count = count + 1;

   resMatch  = matches( radii{1,count}, resName );
   if( resMatch == 0 ) %%no joy here, move along the loop
     continue;
   end

   atomMatch = matches( radii{2,count}, atomName );
   if( atomMatch == 0 ) %%no joy here, move along the loop
     continue;
   end
   
   radius = radii{3, count}; %%have matched a pattern but could still achieve exact match

   if( atomMatch + resMatch == 4 )   %%have achieved exact match - job done
      return;
   end

end


%%%function to compare strings
function matchStatus = matches( mask, nameString )

%%define these return codes for our 'regexp parser'
%NO_MATCH      = 0;
%GENERAL_MATCH = 1;
%EXACT_MATCH   = 2;

numChars    = length( nameString );

if( numChars != length( mask ) )
   matchStatus = 0;   %%mask is longer or shorter than target string therefore no match possible
   return;
end

matchStatus = 2;
count       = 1;

while ( count <= numChars )
   if( mask(count) != nameString(count) )
        if( mask(count) == "*") 
            matchStatus = 1;
        else
            matchStatus = 0;
            return;
        end
   end
   count = count + 1;
end




