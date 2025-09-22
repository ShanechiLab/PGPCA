%{
This is a tutorial for running PGPCA + torus in ICLR 2025 submission. Basically, it includes 3 parts:

1.  Create a true torus model and generate samples from it.
2.  Learn a PGPCA torus model from the training samples.
3.  Plot 3D slices for the probability distributions of the PGPCA torus model.

NOTE:
1.  Add all folder "Function" and its subfolders in the path first.
%}

scriptSW.TorusModelAndSample = 0;
scriptSW.PGPCALearning = 0;
scriptSW.PlotPGPCADistribution = 1;

%% Part 1: Create the torus model + generate random samples.
%{
Keyword can be changed:

1.  probZType   = "uniAng" or "uniTorus".       The p(z) distribution.
2.  probYlZType = "EuCOV" or "GeCOV".           The distribution coordinate K(z).
%}
if(scriptSW.TorusModelAndSample)
    
% Model: the torus & p(z) on top of it.
torusSet = struct('outR',           3, ...
                  'inR',            1, ...
                  'probZType',      "uniAng", ...
                  'numZInte',       [50,20]);
              
% Model: the normal distribution p(y|z).
covSet   = struct('covBasic',       diag([0.1, 0.3, 0.5]), ...
                  'probYlZType',    "EuCOV");
              
% The #sample.
numSamp = 50000;

% -------------------------------------- Create the model. -------------------------------------------
% Generate the model.
[torusObj, inteZPoint, torusFunStru, resColl] = PGPCATorusModel( torusSet, covSet );

% -------------------------------------- Generate samples. -------------------------------------------    
[sampY, sampZ, sampInfo] = PGPCATorusSample( numSamp, torusSet.probZType, torusObj, torusFunStru );

end

%% Part 2: Learn the PGPCA model.
if(scriptSW.PGPCALearning)
    
% Setting: learning p(z) or not.
probZFitSW = true;
% Setting: the type of covariance (or distribution coordinate K(z)) = "EuCOV"/"GeCOV".
covType = "EuCOV";
% Setting: the PGPCA model dimension 0 <= m <= 3.
dimC = 2;

% -------------------------------------- Prepare inputs for PGPCA. -----------------------------------    
% Only use the "funMean & funProbZ" in the true model. 
funStru = struct('funProbZ',    torusFunStru.funProbZ, ...
                 'funMean',     torusFunStru.funMean);

% Select "funMatK" based on the COV type.
dimY = size(sampY, 2);
switch covType
    case "EuCOV"
        funStru.funMatK = @(x) repmat(eye(dimY), [1,1,size(x,1)]);
    case "GeCOV"
        funStru.funMatK = torusFunStru.funMatK;
    otherwise
        error('The "covType" (%s) is invalid.', covType);
end

% Parameter: system variables.
paraStru = struct('matC', [], 'mainVar', [], 'sideVar', []);

% -------------------------------------- Run PGPCA. -------------------------------------------
[dimCNew, paraStruNew, funStruNew, resEStep, resMStep, resLogLL] = ...
    PGPCA.EMAlgo(sampY, inteZPoint, dimC, paraStru, funStru, ...
                 'iterLimit', 40, "saveType", "core", "fitFunProbZSW", probZFitSW);

end

%% Part 3: Plot PGPCA 3D probability distribution with slices.
%{
Since the PGPCA training takes long, we provide two ways to see the result.

1.  From the true model above.
2.  From the matfile we provided.
%}
if(scriptSW.PlotPGPCADistribution)
    
% Setting: the source of plotting variable. ("script"/"matfile")
dataSource = "matfile";

% Setting: the plotting data index (only for "PGPCA_Simu_T2_3D_ICLR")
matInd = 8;

% The XYZ grid for plotting the distribution.
gridCell = {[-6:0.1:6], [-6:0.1:6], [-6:0.1:6]}; 

% The limit of colorbar to normalize the probability of all figures.
cbLimit = [0, 0.012];

% -------------------------------------- Prepare plotting variables. -----------------------------------
switch dataSource
    case "script"
        plotMd       = "given";
        plotZSamp    = inteZPoint;
        plotParaStru = [];
        plotFunStru  = torusFunStru;
        
    case "matfile"
        load("PGPCA_Simu_T2_3D_ICLR.mat", 'matColl');
        plotMd       = "parameter";
        plotZSamp    = matColl(matInd).sampZ;
        plotParaStru = matColl(matInd).paraStruNew;
        plotFunStru  = matColl(matInd).funStruNew;

    otherwise
        error('The given "dataSource" (%s) is invalid.', dataSource);
end

% -------------------------------------- Plot the result. -----------------------------------
[axisH, objHColl, probInfo] = ...
        PGPCA.PlotProbDensity(gridCell, plotMd, plotZSamp, plotParaStru, plotFunStru, ...
                              "tickGap", 1, "cbLimit", cbLimit);

% Axis labels.
xlabel(axisH, 'Y_1');
ylabel(axisH, 'Y_2');
zlabel(axisH, 'Y_3');

end

