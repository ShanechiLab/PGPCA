function areaScale = LocalArea(obj, seleAng, varargin)
%{
This function computes the local area scaling factor at "seleAng".

Algorithm (coarea formula):  areaScale = sqrt( det( matJ'*matJ ) );

NOTE:
1.  This function is a "method" of the class "TorusFitting".
2.  When the torus is standard, the area factor = (R + r*cos(inAng))*r. But considering generalizabiliy,
      we use "coarea formula" to handle all types of tangent spaces w.r.t different torus.
3.  Since dim(torus) = 2 no matter how high the "obj.dim" is, "matJ'*matJ" is always a 2*2 matrix, so the
      computation is fast.
-----------------------------------------------------------------------------------------------------------

Inputs:
--------------- Essential
seleAng:        a 2D array [numSamp, 2]. Each row is a 2D angular coordinate [outAng, inAng] (in radian).
--------------- Optional ('tag' +　'value' with (number))
"compMd":       (1) (compMd) a string "coarea / standFix". (def: "coarea")
                    The computing method of the area scaling factor. They are:
    "coarea":           areaScale = sqrt( det( matJ'*matJ ) ); The coarea formula
    "standFix":         areaScale = (R + r*cos(inAng))*r; The fixed formula for standard torus.

Outputs:
--------------- Essential
areaScale:      a column [numSamp, 1]. The area scaling factor of corresponding "seleAng(i,:)".
%}

%% Assign inputs.
% Extract variables.
outR   = obj.outR;
inR    = obj.inR;
dimObj = obj.dim;

% Constant.
[numSamp, dimAng] = size(seleAng);
% Rename variables for simplification.
inAng  = seleAng(:,2);

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('compMd',  "coarea",  @(x) any(x == ["coarea", "standFix"]));

% Parse the inputs.
optionP.parse(varargin{:});
compMd = optionP.Results.compMd;

% Assert: the torus intrinsic coordinate must be 2D.
assert( dimAng == 2, 'The row dimension of "seleAng" (%d) must be 2.', dimAng ); 

%% Compute area scaling factor based on selected method.
switch compMd
    case "coarea"
        % CASE: coarea formula.
        tangArray = obj.LocalSpace("tangent", seleAng, "normSW", false);
        
        % The symmetric matrix "matJ' * matJ".
        tangInner = pagemtimes(tangArray, 'none', tangArray, 'transpose');
        tangInner = (tangInner + permute(tangInner, [2,1,3]))/2;
        
        % The determinant with squared root.
        areaScale = zeros(numSamp, 1);
        for(j = 1:1:numSamp)
            areaScale(j) = sqrt(det(tangInner(:,:,j)));
        end

    case "standFix"
        % CASE: the fixed formula for a standard torus.
        % Assert: the embedding dimension = 3.
        assert( dimObj == 3, 'The "obj.dim" (%d) must be 3 to use option "standFix".', dimObj );
        
        % Compute the area factor.
        areaScale = (outR + inR * cos(inAng)) * inR;
        
    otherwise
        error('The given "compMd" (%s) is invalid.', compMd);
end

% Assert: areaScale is real & >= 0.
assert( all( isreal(areaScale) ), 'All "areaScale" must be real values.' );
assert( all( areaScale >= 0 ), 'All "areaScale" (min: %2.2e) must be >= 0.', min(areaScale) );

end

