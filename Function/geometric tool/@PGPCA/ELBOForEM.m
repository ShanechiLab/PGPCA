function logLLProbQ = ELBOForEM(typeELBO, sampY, inteZStru, paraStru, funStru)
%{
11-02-2023: (A) Include function "CrossEntropy" to handle the extreme case 0*log(0) = 0 in cost function.
            (B) Add "zeroProbLimit" to handle the numerical problem between "zProb & probYZ".
-----------------------------------------------------------------------------------------------------------

This auxiliary function computes the ELBO of PGPCA EM. Note that this bound depends on EM's type, so we 
  need this function for coding integrity. The steps are:

1.  All [j, k] pairs of p(y_j, z_k).
2.  For each j, compute its ELBO: sum_k( q_j(z_k) * ln(p(y_j, z_k)) ) - sum_k( q_j(z_k) * ln(q_j(z_k)) );

NOTE:
1.  For some types (e.g., "delta_func"), the last term is not necessary.
2.  This is the most general form of ELBO, so the I/O interface is also generalizable!
3.  Function "funStru.funProbZ" is normalized for discrete integration already. No further normalization.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
typeELBO:           a string "posterior / delta_func". The types of ELBO depending on the learning Algo.
sampY:              a 2D array [numY, dimY]. Each [j, :] is a sample y_j.
inteZStru:          a struct. It includes variables for integrating probability w.r.t state z.
    zGp:                a 2D array [numZGp, dimZ]. Each [k, :] is a sample z_k.
    zIndSet:            a 2D array [numY, numZInd]. Each [j, k] is p(y_j, z_inteZIndSet(j,k)).
    zProb:              a 2D array [numY, numZInd]. Each [j, k] is q_j(z_inteZIndSet(j,k)).
paraStru:           a struct. It contains all system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. It collects all function handles in following fields:
    funProbZ:           a func. handle. The probability p(z).
    funProbYlZ:         a func. handle. The probability p(y|z).

Outputs:
------------------- Essential
logLLProbQ:         a column [numY, 1]. Each [j, 1] is the ELBO of log( p(y_j) ).
%}

%% Prepare variables.
% Constant.
zeroProbLimit = PGPCA.zeroProbLimit;

% Extract variables.
zGp     = inteZStru.zGp;
zIndSet = inteZStru.zIndSet;
zProb   = inteZStru.zProb;

% Compute constants.
numY = size(sampY, 1);

%% Compute all pairs of p(y_j, z_k) (= "probYZ").
switch typeELBO
    case "posterior"
        % CASE: q_j(z_k) = posterior distribution.
        probZ   = funStru.funProbZ(zGp);
        probYlZ = funStru.funProbYlZ(sampY, zGp, paraStru);
        probYZ  = probYlZ .* probZ';

    case "delta_func"
        % CASE: q_j(z) = delta function \delta(z - z_j^*) from IVT.
        funProbYlZPair = @(y,z) funStru.funProbYlZ(y, z, paraStru);
        probYZ = PGPCA.PairCondiProb(sampY, zGp(zIndSet, :), funProbYlZPair);

    otherwise
        error('The "typeELBO" (%s) is invalid.', typeELBO);
end

%% Calibrate "zProb" due to numerical error.
%{
This is a special session for the numerical problem: precision limitation. Mathematically, 1e-324 ~= 0.
  But 1e-324 -> 0 in MATLAB! Therefore, when some "probYZ(i,j) = 0", we need to make their corresponding 
  "zProb(i) = 0" to have "0*log(0) = 0". However, "zProb(i)" must be small to justify this manipulation.
%}

if( typeELBO == "posterior" )
    % Due to "zIndSet(j,:)" operation, we need to clean "zProb" row by row.
    for(j = 1:1:numY)
        zeroIndYZ = (probYZ(j, zIndSet(j,:)) == 0);
        
        % Assert: corresponding zProb(j,k) are small enough.
        maxRowZProb = max(zProb(j, zeroIndYZ));
        assert( all(maxRowZProb <= zeroProbLimit), ...
                'The selected "zProb" (max: %2.3e) must be small enough (limit: %2.3e).', ...
                maxRowZProb, zeroProbLimit );
        
        % Put them into zeros.
        zProb(j, zeroIndYZ) = 0;
    end
end

%% Compute ELBO for each y_j.
switch typeELBO
    case "posterior"
        % CASE: posterior distribution.        
        % Compute the 1st term.
        logLL1stTerm = zeros(numY, 1);
        for(j = 1:1:numY)
            logLL1stTerm(j) = CrossEntropy("complete", zProb(j, :), probYZ(j, zIndSet(j,:)), ...
                                           'negaSW',     true, ...
                                           'sumCheckSW', [true, false]);
        end
        
        % Compute the 2nd term.
        logLL2ndTerm = CrossEntropy("complete", zProb, zProb, 'negaSW', true, 'sumCheckSW', true(1,2));

        % Final log-likelihood.
        logLLProbQ = logLL1stTerm - logLL2ndTerm;

    case "delta_func"
        % CASE: delta function \delta(z - z_j^*) from IVT.
        logLLProbQ = log( probYZ );
        assert( all( IsNumericProper(logLLProbQ) ), ...
                '"logLLProbQ" in type "delta_func" must be properly numerical!' );

    otherwise
        error('The "typeELBO" (%s) is invalid.', typeELBO);
end

end

