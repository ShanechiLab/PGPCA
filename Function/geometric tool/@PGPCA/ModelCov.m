function [covVal, resColl] = ModelCov(covMode, sampZ, paraStru, funStru)
%{
This function computes all types of covariances in PGPCA model. It's good for coding integrity.

Model:  y_t = phi(z_t) + K(z_t)*C*x_t + r_t [COV: R -> Ku(z_t) * diag(mainVar, sideVar) * Ku(z_t)']

where   (A) Ku(z_t) = [K(z_t), Kc(z_t)]. (u: unified basis, c: complement basis).
        (B) Basis "K(z_t) & Kc(z_t)" correspond to variance "mainVar & sideVar", respectively.
-----------------------------------------------------------------------------------------------------------

We compute 3 types of covariance:

1.  "covCPart":     covVal = the noise covariance of R(sampZ).   
2.  "covRPart":     covVal = the covariance "K(z_t)*C*x_t" -> K(sampZ) * CC' * K(sampZ)'.
3.  "covY":         covVal = the covariance of "y - phi(z)" = covCPart + covRPart.

NOTE:
1.  No matter what "covMode" is selected, both "covCPart & covRPart" are computed. The output "covVal" is
      just good for anonymous function application.
2.  The I/O interface follows other functions, so it's generalizable.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
sampZ:              a 2D array [numZ, dimZ]. Every row is a sample of z_t.
paraStru:           a struct. It contains all (learned) system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. It collects all function handles in following fields:
    funMatK:            a func. handle. The "orthonormal" coordinate matrix "K(z)" for each "z" sample.

Outputs:
------------------- Essential
covVal:             a 3D array [dimY, dimY, numZ]. The selected COV at each "sampZ" sample.
resColl:            a struct. It records all computing information in following fields:
    covCPart:           a 3D array [dimY, dimY, numZ]. The COV of C part.
    covRPart:           a 3D array [dimY, dimY, numZ]. The COV of R part.
%}

%% Prepare basic/shared parameters.

% Extract values for convenience.
matC    = paraStru.matC;
mainVar = paraStru.mainVar;
sideVar = paraStru.sideVar;

% Generate "matK" for "sampZ" and extract constants.
matKAll = funStru.funMatK(sampZ);
[dimY, ~, numZ] = size(matKAll);

%% Compute: C part.
if( isempty(matC) )
    % CASE: "matC" is not included.
    covCPart = zeros(dimY, dimY, numZ);
else
    % CASE: "matC" is included.
    covCPart = pagemtimes(matKAll, matC);
    covCPart = pagemtimes(covCPart, 'none', covCPart, 'transpose');
    covCPart = (covCPart + permute(covCPart, [2,1,3]))/2;
end

%% Compute: R part.
if( isempty(sideVar) )
    % CASE: All noise is modeled by "mainVar" only.
    covRPart = repmat( eye(dimY)*mainVar, [1,1,numZ] );
else
    % CASE: "sideVar" captures the unmodeled space.
    matKCov = pagemtimes(matKAll, 'none', matKAll, 'transpose');
    matKCov = (matKCov + permute(matKCov, [2,1,3]))/2;
    covRPart = mainVar * matKCov + sideVar * (eye(dimY) - matKCov);
end

%% Prepare outputs.
% Select the COV type.
switch covMode
    case "covCPart"
        covVal = covCPart;

    case "covRPart"
        covVal = covRPart;

    case "covY"
        covVal = covCPart + covRPart;
        
    otherwise
        error('The given "covMode" (%s) is invalid.', covMode);
end

% Collect other info.
resColl = struct('covCPart', covCPart, ...
                 'covRPart', covRPart);

end

