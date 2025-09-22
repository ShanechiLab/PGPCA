function probZ = PGPCATorusFunProbZ( torusObj, sampZ, probZType, varargin )
%{
This function computes the probability p(z) on torus given different probability settings.

Input "probZType" can be:

(1) "uniTorus":     z is uniformly distributed on the "embedded torus surface".
(2) "uniAng":       z is uniformly distributed on the "2D angular covering space" (NOT torus surface!).

NOTE:
1.  This function should keep simple to make it flexible. It's just a basic/auxiliary function.
2.  We add option "normSW" for sanity check.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
torusObj:           a "TorusFitting" object. It's a user-defined torus class.
sampZ:              a 2D array [numZ, 2]. Each [j,:] = [outAng, inAng] is a 2D torus coordinate.
probZType:          a string "uniTorus / uniAng". The type of p(z) on the torus. Options are above.
------------------- Optional ('tag' + 'value' with (number))
"normSW":           (1) (normSW) a boolean. (def: true)
                        "true" means the "probZ" is normalized so sum(probZ) = 1.

Outputs:
------------------- Essential
probZ:              a column [numZ, 1]. The p(z) of "sampZ".
%}

%% Assign inputs.
% Constant.
numZ = size(sampZ, 1);

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('normSW',  true,   @islogical);

% Parse the inputs.
optionP.parse(varargin{:});
normSW = optionP.Results.normSW;

%% Compute "probZ".
switch probZType
    case "uniAng"
        % CASE: uniform "augular" distribution.
        probZ = ones(numZ, 1);
        
    case "uniTorus"
        % CASE: uniform ""
        probZ = torusObj.LocalArea(sampZ, "compMd", "coarea");
        
    otherwise
        error('The given "probZType" (%s) is invalid.', probZType);
end

% (if asked) normalize p(z).
if( normSW )
    probZ = probZ / sum(probZ);
end

end

