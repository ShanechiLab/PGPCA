function [vecArray, resColl] = LocalSpace(obj, spaceType, seleAng, varargin)
%{
This function computes the various local spaces at "seleAng". "spaceType" is for multi-functionality.

(A) Model: (R > r)

    [outAng, inAng] -> [(R + r*cos(inAng))*cos(outAng), (R + r*cos(inAng))*sin(outAng), r*sin(inAng)].

(B) spaceType = 

    1.  "tangent":      vecArray = [2, obj.dim, numSamp]. The tangent vector w.r.t "outAng & inAng".
    2.  "norm-D":       vecArray = [1, obj.dim, numSamp]. The norm-D vector orthogonal to "tangent".
    3.  "complete":     vecArray = [3, obj.dim, numSamp]. It combines [tangent; norm-D] together.

NOTE:
1.  This function is a "method" of the class "TorusFitting".
2.  "spaceType" integrates various local coordinate sets in one function.
3.  The format of "vecArray" is always 3D, even when dim(norm-D) = 1. This is for handling consistency and
      future generalization (e.g., if this torus is embedded in R^10, then dim(norm-D) = 10 - 2 = 8).
-----------------------------------------------------------------------------------------------------------

Inputs:
--------------- Essential
spaceType:      a string. It decides the type of local space vectors. Options are above.
seleAng:        a 2D array [numSamp, 2]. Each row is a 2D angular coordinate [outAng, inAng] (in radian).
--------------- Optional ('tag' +　'value' with (number))
"normSW":       (1) (normSW) a boolean. (def: true)
                    "true" means all vectors "||vecArray(i,:,j)|| = 1" are normalized.

Outputs:
--------------- Essential
vecArray:       a 3D array [numVec, obj.dim, numSamp]. Each [:,:,j] is a local space at "seleAng(j,:)".
resColl:        a struct. It includes all other results in this function as following:
    spaceType:      a string. It's equal to input "spaceType". Just for record.
    normSW:         a boolean. It's equal to option "normSW". Just for record.
    vecScale:       a 2D array [numSamp, numVec]. Each [i,j] is the scaling factor of vecArray(j,:,i).
%}

%% Assign inputs.
% Extract variables.
outR   = obj.outR;
inR    = obj.inR;
dimObj = obj.dim;

% Constant.
[numSamp, dimAng] = size(seleAng);
% Rename variables for simplification.
outAng = seleAng(:,1);
inAng  = seleAng(:,2);

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('normSW',  true,   @islogical);

% Parse the inputs.
optionP.parse(varargin{:});
normSW = optionP.Results.normSW;

% Assert: the torus intrinsic coordinate must be 2D.
assert( dimAng == 2, 'The row dimension of "seleAng" (%d) must be 2.', dimAng ); 
% The norm-D dimension.
dimnormD = dimObj - dimAng;

%% Decide the included types of vectors.
switch spaceType
    case "tangent"
        % CASE: tangent.
        tangSW = true;
        normDSW = false;
        
    case "norm-D"
        % CASE: norm-D.
        tangSW = false;
        normDSW = true;
        
    case "complete"
        % CASE: [tangent; norm-D].
        tangSW = true;
        normDSW = true;
        
    otherwise
        error('The given "spaceType" (%s) is invalid.', spaceType);
end

%% Compute: tangent vectors.
if( tangSW )
    % The scaling factor.
    tangScale = [outR + inR * cos(inAng), inR * ones(numSamp, 1)];
    % Assert: all scales are positive.
    assert( all(tangScale > 0, 'all'), ...
            '"tangScale" (min: %2.2e) must > 0.', min(tangScale, [], 'all') );
    
    % The normalized tangent space.
    tangSpace = cat(3, [-sin(outAng)            ,  cos(outAng)            , zeros(numSamp, 1)], ...
                       [-sin(inAng).*cos(outAng), -sin(inAng).*sin(outAng), cos(inAng)       ]);
    
    % If no normalization, multiply the scale back.
    if( ~normSW )
        for(j = 1:1:dimAng)
            tangSpace(:,:,j) = tangSpace(:,:,j) .* tangScale(:,j);
        end
    end
    
    % Organize the 3D array.
    tangSpace = permute(tangSpace, [3,2,1]);
    
else
    % No computation needed.
    tangScale = [];
    tangSpace = [];
    
end

%% Compute: norm-D vectors.
if( normDSW )
    % The scaling factor.
    normDScale = (outR + inR * cos(inAng)) * inR;
    % Assert: all scales are positive.
    assert( all(normDScale > 0, 'all'), ...
            '"normDScale" (min: %2.2e) must > 0.', min(normDScale, [], 'all') );
    
    % The normalized norm-D space.
    normDSpace = [cos(inAng).*cos(outAng), cos(inAng).*sin(outAng), sin(inAng)];
    
    % If no normalization, multiply the scale back.
    if( ~normSW )
        for(j = 1:1:dimnormD)
            normDSpace(:,:,j) = normDSpace(:,:,j) .* normDScale(:,j);
        end
    end

    % Organize the 3D array.
    normDSpace = permute(normDSpace, [3,2,1]);
    
else
    % No computation needed.
    normDScale = [];
    normDSpace = [];

end

%% Prepare outputs.
% The final local space.
vecArray = cat(1, tangSpace, normDSpace);

% Collect information.
resColl = struct('spaceType',   spaceType, ...
                 'normSW',      normSW, ...
                 'vecScale',    []);

resColl.vecScale = cat(2, tangScale, normDScale);

end

