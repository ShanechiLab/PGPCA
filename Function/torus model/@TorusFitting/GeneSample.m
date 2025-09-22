function [sampY, sampZ, resColl] = GeneSample(obj, sampType, sampPara)
%{
This function generates various types of samples on the torus. Type "sampType" can be:

1.  "lattice":      sampPara = a row [1, 2] = [numSampOut, numSampIn]. Lattice along outR/inR.
2.  "given":        sampPara = a 2D array [numSamp, 2]. The given 2D torus angular coordinate.
3.  "uniTorus":     sampPara = +integer. The #sample distributed uniformly "on the torus".
4.  "uniAng":       sampPara = +integer. The #sample distributed uniformly "on the 2D angular plane".

NOTE:
1.  This function handle both static/random samples on the torus for convenience & coding integrity.
2.  The type "given" is for unifying the user-given case in many applications.
3.  For each case, I rename "sampPara" to make the code readable, or the "sampPara" meaning is unclear.
-----------------------------------------------------------------------------------------------------------

Inputs:
--------------- Essential
sampType:       a string "lattice / given / uniTorus / uniAng". The sampling types. Details are above.
sampPara:       a parameter. It depends on "sampType". Details are above.

Outputs:
--------------- Essential
sampY:          a 3D array [numSamp, 3]. Each [j,:] is an embedded sample based on "sampZ(j,:)".
sampZ:          a 2D array [numSamp, 2]. Each [j,:] is a sample on the 2D angular space.
resColl:        a struct. It includes all other information in following fields:
    numSamp:        a +integer. The #sample in total.
    lattZOut:       a column [numSampOut, 1]. The lattice Z along "outR" (only "sampType = lattice").
    lattZIn:        a column [numSampIn, 1]. The lattice Z along "inR" (only "sampType = lattice").
%}

%% Prepare variables.
% Extract variables.
outR   = obj.outR;
inR    = obj.inR;

% Constant: scaling factor for sampling points.
rejSampFactor = 4;

% Initialize output parameters.
resColl = struct('numSamp', [], 'lattZOut', [], 'lattZIn', []);

%% Generate samples.
switch sampType
    case "lattice"
        %% CASE: Build lattice grid along outR & inR.
        % Rename variable.
        numZLatt = sampPara;
        clear('sampPara');
        
        % The lattice points along the outer/inner radius.
        lattZOut = [0:1:numZLatt(1)-1]' / numZLatt(1) * (2*pi);
        lattZIn  = [0:1:numZLatt(2)-1]' / numZLatt(2) * (2*pi);

        % Compute "sampZ".
        [meshOut, meshIn] = meshgrid(lattZOut, lattZIn);
        sampZ = [meshOut(:), meshIn(:)];
        
        % Prepare specific outputs.
        resColl.lattZOut = lattZOut;
        resColl.lattZIn  = lattZIn;
        
    case "given"
        %% CASE: The "sampZ" is given.
        % Rename variable.
        sampZ = sampPara;
        clear('sampPara');
        
    case "uniTorus"
        %% CASE: Random samples on the "torus" uniformly.
        % Rename variable.
        numSamp = sampPara;
        clear('sampPara');
        
        % Compute "sampZOut".
        sampZOut = ShiftRand([0,2*pi], [numSamp,1]);
        
        % "sampZIn" is distorted in the manifold mapping, so it must be balanced by "sampling rejection".
        numOverSamp = numSamp * rejSampFactor;
        % The candidate of "sampZIn" and its gated number
        sampZInCandi = ShiftRand([0,2*pi], [numOverSamp,1]);
        sampZInGate  = (outR + inR * cos(sampZInCandi)) / (outR + inR);
        % The random rejecting samples to decide which "sampZInCandi" is kept!
        rejRandSamp = rand([numOverSamp, 1]);
        sampZIn = sampZInCandi( rejRandSamp <= sampZInGate );
        
        % Assert: the #sampZIn is enough.
        numSampZIn = length(sampZIn);
        if( numSampZIn >= numSamp )
            % #sampZIn is enough. We just cut it into the correct length.
            sampZIn = sampZIn(1:numSamp, 1);
        else
            % #sampZIn is not enough! We show an error not resampling to avoid infinite iteration.
            error('The #accepted "sampZIn" (%d over %d) is too little (need %d).', ...
                  numSampZIn, numOverSamp, numSamp);
        end
        
        % Bind the result.
        sampZ = [sampZOut, sampZIn];
        
    case "uniAng"
        %% CASE: Random samples on the "2D angular plane" uniformly.
        % Rename variable.
        numSamp = sampPara;
        clear('sampPara');
        
        % Compute "sampZ" by shifted MATLAB function "rand".
        sampZ = ShiftRand([0,2*pi], [numSamp,2]);
        
    otherwise
        error('The given "sampType" (%s) is invalid.', sampType);
end

%% Prepare "shared" outputs.
% Compute "sampY".
sampY = obj.Position(sampZ);

% The #sample.
resColl.numSamp = size(sampZ, 1);

end

