function [sampY, sampZ, resColl] = PGPCATorusSample( numSamp, probZType, torusObj, funStru )
%{
This function generates random samples given a "torus-based" PGPCA model. There are two steps:

1.  Sample z:       "z" is a 2D torus state following the distribution p(z).
2.  Sample y|z:     "y|z" follows a Gaussian distribution p(y|z).

The sampling methods of "z" on the torus are:

1.  "uniTorus":     "z" is uniformly distributed on the "embedded torus surface".
2.  "uniAng":       "z" is uniformly distributed on the "2D angular covering space" (NOT torus surface!).

The I/O of function handles in inputs:

1.  funCOV:         [dimY, dimY, numZ] = funCOV(z);         Each [:,:,j] is the COV at z(j,:).
2.  funProbY:       [numY, 1]          = funProbY(y);       The output is p(y).

NOTE:
1.  We sample "z" by function "torusObj.GeneSample" not function handle "funProbZ" from "PGPCATorusModel"
      because "funProbZ" is a discretized approximation of the continuous z distribution on the torus, so
      we sample "z" by a real continuous function.
2.  Field "funStru.funProbY" is optional so this function is more general.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
numSamp:            a +integer. The #sample.
probZType:          a string "uniTorus / uniAng". The type of p(z). Details are above.
torusObj:           a "TorusFitting" object. The torus object.
funStru:            a struct. It includes all functions for generating samples & probability below:
    funCOV:             a func. handle. The Gaussian COV of p(y|z) at z. I/O is above.
op  funProbY:           a func. handle. The p(y) of each function. I/O is above.

Outputs:
------------------- Essential
sampY:              a 2D array [numSamp, dimY]. Each [j,:] is a "PGPCA + torus" sample y(j).
sampZ:              a 2D array [numSamp, 2]. Each [j,:] is a "torus" sample z(j).
resColl:            a struct. It includes all other information in following fields:
    probZType:          a string. It's equal to input "probZType" just for recording.
    sampZEmbed:         a 2D array [numSamp, dimY]. The embedded torus coordinate w.r.t "sampZ".
    sampCov:            a 3D array [dimY, dimY, numSamp]. Each [:,:,j] is the COV at "sampZ(j,:)".
    sampYlZ:            a 2D array [numSamp, dimY]. Each [j,:] is sampled from a local Gaussian prob.
    probY:              a column [numSamp, 1]. The probability of "sampY".
%}

%% Compute: Y = Z + Y|Z
% Basic variables.
dimY = torusObj.dim;

% Compute: Z.
[sampZEmbed, sampZ] = torusObj.GeneSample( probZType, numSamp );

% Gaussian Cov at Z.
sampCov = funStru.funCOV( sampZ );

% Compute: Y|Z.
sampYlZ = mvnrnd( zeros(numSamp, dimY), sampCov );

% Compute: Y = Z + Y|Z.
sampY = sampZEmbed + sampYlZ;

% (if given) compute p(y).
if( isfield(funStru, 'funProbY') )
    probY = funStru.funProbY( sampY );
else
    probY = [];
end

%% Prepare outputs.
resColl = struct('probZType',       probZType, ...
                 'sampZEmbed',      sampZEmbed, ...
                 'sampCov',         sampCov, ...
                 'sampYlZ',         sampYlZ, ...
                 'probY',           probY);

end

