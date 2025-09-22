function vecPosi = Position(obj, seleAng)
%{
This function computes the torus embedding position at "seleAng".

Model: (R > r)

[outAng, inAng] -> [(R + r*cos(inAng))*cos(outAng), (R + r*cos(inAng))*sin(outAng), r*sin(inAng)].

NOTE:
1.  This function is a "method" of the class "TorusFitting".
2.  Since the torus model is not polynomial-based, we separate the "position (0 order)" and the 1st order
      (tangent & normal vector) functions.
-----------------------------------------------------------------------------------------------------------

Inputs:
--------------- Essential
seleAng:        a 2D array [numSamp, 2]. Each row is a 2D angular coordinate [outAng, inAng] (in radian).

Outputs:
--------------- Essential
vecPosi:        a 2D array [numSamp, obj.dim]. Each row is the 3D position of corresponding "seleAng".
%}

%% Sanity check.
% Extract variables.
outR   = obj.outR;
inR    = obj.inR;
dimObj = obj.dim;

% Constant.
[numSamp, dimAng] = size(seleAng);
% Rename variables for simplification.
outAng = seleAng(:,1);
inAng  = seleAng(:,2);

% Assert: the torus intrinsic coordinate must be 2D.
assert( dimAng == 2, 'The row dimension of "seleAng" (%d) must be 2.', dimAng ); 

%% Compute the position.
vecPosi = zeros(numSamp, dimObj);
vecPosi(:,1) = (outR + inR * cos(inAng)) .* cos(outAng);
vecPosi(:,2) = (outR + inR * cos(inAng)) .* sin(outAng);
vecPosi(:,3) = inR * sin(inAng);

end

