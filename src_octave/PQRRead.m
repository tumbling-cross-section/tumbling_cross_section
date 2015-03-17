function PDB_struct = PQRRead(filename, modelNumber)

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
%% This function reads in a pqr file (for our purposes a pdb wihch also contains atomic radii)
%% and returns a structure containing the :
%% coordinates, radii, number of atoms and longest axis length of the structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check arguments and open infile
switch( nargin )
     case 0
       disp( "PQRRead:\n requires filename \n optionally accepts number of model in file."  );
     case 1
       fid         = fopen(filename,'r');
       modelNumber = 1;
     case 2
       fid        = fopen(filename,'r');
     otherwise
       error( 'ERROR! PQRRead called with incorrect arguments' );
       return
end

if fid == -1,
    error('Unable to open specified file');
    return
else
    PDB_struct = [];
    NumOfATOM = 0;
    maxRadius = 0;

end

disp('Reading PQR file............');

readingModel = 1;
while ( readingModel <= modelNumber )
    
    tline = fgetl(fid);
    
    if ~ischar(tline)
        break; % For end of file recognition
    end
    
    if size(tline)>0 % Omit the empty lines to avoid error of invalid matrix index.
             
            sz = size(tline);
            tline = [tline blanks(80-sz(2))];   % Trim input line to 80 chars. 
            Record_name = upper(tline(1:6));    % Assuming that the record name will be left-alligned
            Record_name = deblank(Record_name); 
    
            switch Record_name
                case 'ATOM'

		    if( readingModel == modelNumber )
    
			%%increment the count the number of atoms and pick up coordinates and radii								       
			NumOfATOM = NumOfATOM+1;
			crds(NumOfATOM,:) = [ str2num(tline(31:38)), str2num(tline(39:46)), str2num(tline(47:54))];
			atomicRadius(NumOfATOM) = str2num(tline(63:70));

    			radiusFromCentre         = crds(NumOfATOM,1)**2 +crds(NumOfATOM,2)**2 + crds(NumOfATOM,3)**2;			       
			if( radiusFromCentre >= maxRadius )
			   maxRadius = radiusFromCentre;
			end

                    end			
               		        
                case 'ENDMDL'
		    readingModel = readingModel + 1; %% we have the next structure in an NMR pdb file.

                otherwise  %%If its not an atom we aren't interested
                  %%  disp('skipping non-atom record');
                  %%  disp( Record_name );												       
												       
            end % of the SWITCH statement 
     end % of the IF statement checking the empty string
end % of the WHILE loop

PDB_struct.nAtoms         =  NumOfATOM;	 	
PDB_struct.crds           =  crds;   	 
PDB_struct.atomicRadii    =  atomicRadius;  
PDB_struct.L              =  sqrt( maxRadius );
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

