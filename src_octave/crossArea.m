function area = crossArea( nThetaSteps, nPhiStepsMax, gasRadius, fileName, radiusFilename )

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

%%seed the random number generator for repeatability
%%rand("seed", 42);



largestRadius = 3.0; %%don''t expect to see many atoms bigger than this in biology
minGuesses    = 100; %%this threshold stops the code from exiting before the measure of estimated error is stable.

if( nargin == 4 )
   fprintf( stdout,  "No radius file provided, so reading coordinates and radii from file \"%s\"\n", fileName);
   PDB_struct = PQRRead( fileName, 1 );
else
   fprintf( stdout,  "Reading coordinates from file \"%s\" and radii from \"%s\"\n", fileName, radiusFilename);
   PDB_struct = getPdbStructure( fileName, radiusFilename, 1 );
end

%%adjust the size of the box to allow a probe atom to fit in around the edges
PDB_struct.L = PDB_struct.L + 2 * largestRadius + 2 * gasRadius;

%%let the user know whats happening
fprintf( stdout, "\nStarting Monte-Carlo area measurement.\n" );

angleCount = 1;

for thetaStep = 1:nThetaSteps
   theta = pi * thetaStep / nThetaSteps;

   %% rotation matrix around x-axis
   Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta) ];
   RStheta = PDB_struct.crds(:,1:3) * Rx;

   nPhiSteps = floor( nPhiStepsMax * sin( theta ) ); %%sample evenly over the unit sphere
   if( nPhiSteps == 0 )
     nPhiSteps = 1; %%zero steps is a bit funny
   end

   for phiStep = 1:nPhiSteps
       phi        = 2 * pi * phiStep / nPhiSteps;
       Ry_andProj = [ cos(phi) 0 sin(phi); 0 1 0; 0 0 0];
       RS_proj    = Ry_andProj * RStheta'; %%' gives matrix transpose

       hitCount = 0;

       fprintf( stdout,  "Projecting at angles: %g %g\n", theta, phi );

       for guessIndex = 1:floor( 1 * PDB_struct.L**2 )     

           point = unifrnd ( -0.5 * PDB_struct.L, 0.5 * PDB_struct.L, 2, 1);
           
           %%Test that the point is not in the corners of the square
	   pointFromOrigin = ( point( 1 )**2 + point( 2 )**2 );
	   if( pointFromOrigin <=  ( largestRadius + gasRadius + 0.5 * PDB_struct.L  )**2 )
	
             for atomIndex = 1:PDB_struct.nAtoms 

	         vectorDistance  = point - RS_proj( 1:2 , atomIndex );
	         clearance       = PDB_struct.atomicRadii( atomIndex ) + gasRadius;

		 %%Test for a collision
		 %% sqrt is numerically expensive
		 sqDistance = vectorDistance( 1 )**2 + vectorDistance( 2 )**2;
 
                 if( sqDistance <= ( clearance**2 ) )
		      hitCount += 1;
		      break;
		 end %%if

              end          %%for
           end             %%if


           if( rem( guessIndex, 100) == 0 ) %% test for convergence every ten guesses


              if( guessIndex >= minGuesses )

                  p_hit          = hitCount / guessIndex;  
                  stdDevEstimate = sqrt( p_hit * ( 1 - p_hit ) / guessIndex ); %%estimated std deviation of a binomial distribution
                  errorRatio     = stdDevEstimate / p_hit;
                 
                  fprintf( stdout, "iteration: %i estimated error ratio: %f estimated area: %f\n", guessIndex, errorRatio, p_hit * (PDB_struct.L**2) );

                  if( errorRatio < 0.01 && p_hit != 0.0 && p_hit != 1.0 )
                  
                     break; %%stop making guesses, we have a reliable estimate

                  end
              end
           end
       end      %%big for

       %%%binomial distribution: expected number of hits  = guessIndex * p_hit
       %%%unbiased estimate of p_hit                      = hitCount / guessIndex;

       %%%estimated variance( hitCount )                 ~= p_hit * ( 1 - p_hit )/guessIndex
       %%%so for 2DP accuracy, require sqrt(variance)/p_hit << 0.01
       
       p_hit = hitCount / guessIndex;
       stdDevEstimate = sqrt( p_hit * ( 1 - p_hit ) / guessIndex );
     
     
       areaEstimate = (PDB_struct.L**2) * p_hit;
       stdDevEstimate =  (PDB_struct.L**2) * stdDevEstimate;

       A( angleCount, 1:4) = [ theta, phi, areaEstimate, stdDevEstimate ];
       angleCount = angleCount + 1;

       fprintf( stdout, "Calculation converged for this set of angles. Area: %g Estimated Error: %g\n", areaEstimate, stdDevEstimate);

   end
end

%% work out mean cross sectional area
meanArea  = mean( A(:, 3) ); 
meanError = mean( A(:, 4) );
fprintf( stdout, "Mean Area Over All Projections: %g\nTotal ESE: %g\n", meanArea, meanError / sqrt(length(A) - 1.0) );
