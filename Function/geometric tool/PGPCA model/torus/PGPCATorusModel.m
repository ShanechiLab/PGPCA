function [torusObj, inteZPoint, funStru, resColl] = PGPCATorusModel( torusSet, covSet )
%{
This function creates a torus-based PGPCA stochastic model. It includes all probability functions.

PGPCA model:  y = g(z) + K(z)*C*x + r (COV: R = (sigma^2)*I_n)

      where:  z ~ [outAng, inAng], the torus 2D coordinate following a probability distribution on torus.
              x ~ N(0, I_m), a multivariate normal distribution.
-----------------------------------------------------------------------------------------------------------

True model:   y = g(z) + w(z)

     where:   z ~ [outAnd, inAng], the torus 2D coordinate as above.
              w ~ N(0, W(z)). W(z) can be a location-dependent/independent COV (e.g., GeCOV / EuCOV).

Types of p(z):

1.  "uniTorus":     z is uniformly distributed on the "embedded torus surface".
2.  "uniAng":       z is uniformly distributed on the "2D angular covering space" (NOT torus surface!).

Type of p(y|z):

1.  "EuCOV":        p(y|z) COV follows the extrinsic coordinate, the Euclidean one R^n.
2.  "GeCOV":        p(y|z) COV follows the instinsic coordinate, the geometric one.
-----------------------------------------------------------------------------------------------------------

The output func. handles of this torus model:

(A) z-dependent function.
    (1) funCOV:         [dimY, dimY, numZ] = funCOV(z);         Each [:,:,j] is the COV at z(j,:).

(B) For PGPCA inputs.
    (1) funProbZ:       [numZ, 1]          = funProbZ(z);       The output is p(z).
    (2) funMean:        [numZ, dimY]       = funMean(z);        Each [j,:]' is g(z(j,:)).
    (3) funMatK:        [dimY, dimL, numZ] = funMatK(z);        Each [:,:,j] is K(z(j,:)).

(C) Other prob. functions.
    (1) funProbYlZ:     [numY, numZ]       = funProbYlZ(y, z);  Each [i,j] is p(y(i,:)|z(j,:)).
    (2) funProbY:       [numY, 1]          = funProbY(y);       The output is p(y).
-----------------------------------------------------------------------------------------------------------

Since "probY = sum_z(probYlZ * probZ)", some torus points "z" must be selected smartly to approximate the
  integral efficiently. We use "numZInte" to sample a "meshgrid" on torus:

1.  numInteOut:     uniformly sample along outer radius [0, outR].
2.  numInteIn:      uniformly sample along inner radius [0, inR].

So the #total sampling "numInteAll = numInteOut * numInteIn". Make sure that it's not too large!

NOTE:
1.  We separate the "modeling & sampling" into two functions since their mechanisms are different.
2.  The PGPCA model is a low-D version of the true model, like low-D PCA vs. true samples.
3.  The integrating z (centInte) is smartly selected on the torus uniformly, and weights are determined
      by p(z). This way models different continuous p(z) on the torus efficiently.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
------------------- (ABOUT TORUS)
torusSet:           a struct. It determines the torus part by following fields:
    outR:               a +scalar. The outer radius of the torus.
    inR:                a +scalar. The inner radius of the torus. (note: outR > inR)
    probZType:          a string "uniTorus / uniAng". The type of p(z) on the torus. Options are above.
    numZInte:           a row [1, 2] = [numInteOut, numInteIn]. The integrating Z points for "probY".
------------------- (ABOUT NORMAL COV)
covSet:             a struct. It determines the attached normal distribution by following fields:
    covBasic:           a 3D array [3, 3]. The basic COV matrix (w/o local coordinate transformation).
    probYlZType:        a string "EuCOV / GeCOV". It decides the COV coordinate type. Options are above.

Outputs:
------------------- Essential
torusObj:           a "TorusFitting" object. It's a user-defined torus class.
inteZPoint:         a 2D array [numInteAll, 2]. Each [j,:] is a torus integrating point for "probY".
funStru:            a struct. It includes all func. handles in following fields (details are above):
    funProbZ:           a func. handle. The probability p(z).
    funMean:            a func. handle. The manifold embedding coordinate g(z).
    funMatK:            a func. handle. The local COV coordinate at z.
    funCOV:             a func. handle. The local COV at z.
    funProbYlZ:         a func. handle. The conditional probability p(y|z).
    funProbY:           a func. handle. The probability p(y).
resColl:            a struct. It records other information in this function in following fields:
    inteZProb:          a column [numInteAll, 1]. The probability p(inteZPoint).
    inteZMean:          a 2D array [numInteAll, dimY]. Each [j,:] is a torus position at "inteZPoint".
    inteZCov:           a 3D array [dimY, dimY, numInteAll]. Each [:,:,j] is a COV at "inteZPoint".
%}

%% Extract variables for simplicity.
% About the torus.
outR      = torusSet.outR;
inR       = torusSet.inR;
probZType = torusSet.probZType;
numZInte  = torusSet.numZInte;

% About the normal COV.
covBasic    = covSet.covBasic;
probYlZType = covSet.probYlZType;

%% The torus model & related functions.
% The torus object.
torusObj = TorusFitting(outR, inR);

% Functions of the torus.
funMean = @(z) torusObj.Position( z );
funMatK = @(z) PGPCATorusFunMatK( torusObj, z );

% p(z) on the torus.
funProbZ = @(z) PGPCATorusFunProbZ( torusObj, z, probZType, "normSW", true );

%% Add normal distribution on the torus.
% Function: z-dependent (maybe) COV.
funCOV = @(z) PGPCAFunNormalCOV( covBasic, z, probYlZType, ...
                                 @(x) torusObj.LocalSpace("complete", x, "normSW", true) );

% Conditional probability p(y|z).
funProbYlZ = @(y, z) MVNPDFBatch(y, funMean(z)', funCOV(z));

%% Special: marginal probability p(y).
% The integrating points along the outer/inner radius.
inteZOut = [0:1:numZInte(1)-1] / numZInte(1) * (2*pi);
inteZIn  = [0:1:numZInte(2)-1] / numZInte(2) * (2*pi);

% Create "inteZPoint".
[meshOut, meshIn] = meshgrid(inteZOut, inteZIn);
inteZPoint = [meshOut(:), meshIn(:)];

% We compute all "inteZPoint"-related values first to save the running time of "funProbY".
inteZProb = funProbZ(inteZPoint);
inteZMean = funMean(inteZPoint);
inteZCov  = funCOV(inteZPoint);

% Marginal probability p(y).
funProbY = @(y) MVNPDFBatch(y, inteZMean', inteZCov) * inteZProb;

%% Prepare outputs.
% Collect all functions.
funStru = struct('funProbZ',        funProbZ, ...
                 'funMean',         funMean, ...
                 'funMatK',         funMatK, ...
                 'funCOV',          funCOV, ...
                 'funProbYlZ',      funProbYlZ, ...
                 'funProbY',        funProbY);

% Collect other information.
resColl = struct('inteZProb',       inteZProb, ...
                 'inteZMean',       inteZMean, ...
                 'inteZCov',        inteZCov);

end

