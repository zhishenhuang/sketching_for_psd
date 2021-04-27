function [fcn,S] = sketch( m, M, typeOfSketch, performTest, varargin )
% fcn = sketch( m, M, 'type' )
%   returns a function that implements a sketch of the requested type
%   so that Y = fcn(A) is a sketched version of A.
%   Every time you call this function, you get a new random sketch
%   Y has m rows, A has M rows. The number of columns of A is arbitrary
%   (Y will have the same number of columns)
%
% Valid types of sketches:
%   gaussian, haar, count, fjlt, hadamard, sparse, subsample, fd (frequent
%   direction), precondhd (preconditioning + uniform sampling)
%
% [fcn,S]   = sketch( ... )
%   also returns the explicit matrix representation of the sketch
%   e.g., fcn(A) is the same as S*A.  S is of size m x M
%
% sketch( ..., performTest )
%   if performTest = true, then this checks that the sketch
%   really has the property E[ S'S ] = I_M
%   errorHistory = sketch( ..., performTest )
%       will return the full history of the errrors
%
% sketch( ..., parameterName, parameterValue, ... )
%   allows optional parameter, such as:
%       'sparsity' (for sparse sketches)
%       'weights'  (for subsample, if you want it non-uniform)
%       'nReps'    (for how many repetitions to use when testing)
%
% Stephen Becker, Feb 2019
% modified by Leo Huang to add frequent direction sketching scheme

S   = [];
if nargin < 4 || isempty(performTest)
    performTest = false;
end


prs = inputParser;
addParameter(prs,'sparsity',.01);
addParameter(prs,'weights',[]);
addParameter(prs,'nReps',100);
addParameter(prs,'tensormode',false);


parse(prs,varargin{:});
sparsity     = prs.Results.sparsity;
weights      = prs.Results.weights;
nReps        = prs.Results.nReps;
tensormode   = prs.Results.tensormode;


if performTest
    sumS    = zeros(M);
    if nargout > 0
        errHist     = zeros(nReps,1);
    end
    fprintf('\nRunning test to see of E[S''S] = I (for sketch of type %s)\n', typeOfSketch);
    printEvery  = round( nReps / 10 );
    for rep = 1:nReps
        % Call this own function recursively
        [~,S]   = sketch( m, M, typeOfSketch );
        sumS    = sumS + S'*S;
        if nargout > 0
            errHist(rep) = norm( sumS/rep - eye(M), 'fro' )/M;
        end
        if ~mod(rep,printEvery)
            err = norm( sumS/rep - eye(M), 'fro' )/M;
            fprintf('%3d trials, error || sampleMean(S''S)-I ||_F is %4.1e', ...
                rep, err );
            if rep > printEvery
                fprintf(', %.2f change', err/errOld);
            end
            fprintf('\n');
            errOld = err;
        end
    end
    sumS    = sumS/nReps;
    fprintf('The first 5 x 5 block of sampleMean is: \n');
    disp( sumS(1:5,1:5) );
    fprintf('Average diagonal entry is %.3f, should be 1\n', mean(diag(sumS)) );
    if nargout > 0
        fcn     = errHist;
    end
    return;
end

switch lower(typeOfSketch)
    
    case 'fd'
        if mod(m,2)~=0
            m = m+1;
        end
        fcn = @(A) fastFD(A,m);
        S   = [];
    
    case 'gaussian'
        S       = randn(m,M)/sqrt(m);
        if tensormode
            fcn     = @(A,MODE) ttm(A, S ,MODE);
        else
            fcn     = @(A) S*A;
        end

    case 'haar'  % see http://arxiv.org/abs/math-ph/0609050  by Mezzadri

        [Q,R]   = qr( randn(M,m), 0 );
        d       = sign(diag(R));
        Q       = Q*spdiags(d,0,m,m);
        S       = sqrt(M/m)*Q';
        if tensormode
            fcn     = @(A,MODE) ttm(A, S, MODE);
        else
            fcn     = @(A) S*A;
        end
        
    case {'count', 'countsketch'}
        if tensormode
            % use slowCountSketch algo
            d       = sign(randn(M,1));
            fcn     = @(A,MODE) ttm(A,spdiags(d,0,M,M),MODE);
            S       = randi(m,M,1); 
        else
            d       = sign(randn(M,1));
            D       = spdiags(d,0,M,M); % bsxfun() is another efficient way to do this
            useTranspose    = true;
            indx_map        = int64(randi(m,M,1)); % don't do this in C!
            if exist( 'countSketch_BLAS', 'file' )
                fcn     = @(A) countSketch_BLAS(A'*D,indx_map,m,useTranspose)';
            elseif exist( 'countSketch', 'file' )
                fcn     = @(A) countSketch(A'*D,indx_map,m,useTranspose)';
            else
                fcn     = @(A) slowCountSketch( D*A, double(indx_map) );
            end
        end

    case 'fjlt'
        if tensormode
            d       = sign(randn(M,1));
            HD      = dctmtx(M) * spdiags(d,0,M,M); % bsxfun() is another efficient way to construct D
            S       = randperm( M, m ); % indices permuation for the target dimension      
            fcn     = @(A,MODE) sqrt(M/m) * ttm( A ,HD ,MODE);
        else
            d       = sign(randn(M,1));
            D       = spdiags(d,0,M,M); % bsxfun() is another efficient way to do this
            ind     = randperm( M, m );
            subsample   = @(X) X(ind,:);
            fcn     = @(A) sqrt(M/m)*subsample( dct( D*A ) ); % FIXME
        end
        
    case {'fjlt_hadamard','hadamard'} % Hadamard version of FJLT
        M2  = 2^nextpow2(M);
        if M ~= M2
            % need to zero pad
            upsample = @(X) [X; zeros(M2 - M, size(X,2) ) ];
        else
            upsample = @(X) X; % do nothing
        end
        
        d       = sign(randn(M2,1));
        D       = spdiags(d,0,M2,M2); % bsxfun() is another efficient way to do this
        m       = round(m * M2/M); % adjust m to maintain compress ratio, LH, Jan 20
        ind     = randperm( M2, m );
        subsample   = @(X) X(ind,:);
        
        if exist('hadamard_pthreads','file')==3 
            fcn     = @(x) 1/sqrt(m)*subsample( hadamard_pthreads( D*upsample(full(x))) );
        elseif exist('hadamard','file')==3
            fcn     = @(x) 1/sqrt(m)*subsample( hadamard( D*upsample(full(x))) );
        elseif exist('Hadamard_teaching_code','file')==2
            % It turns out our naive Matlab implementation is better than
            %   the fwht function!
            fcn     = @(x) 1/sqrt(m)*subsample( Hadamard_teaching_code( ...
                D*upsample(full(x)) ) );
        else
            % This is slow!
            fcn     = @(x) (M2*sqrt(1/m))*subsample( fwht( D*upsample(full(x)), [], 'hadamard') );
        end
        
    case 'sparse'
        S   = sign(sprandn(m,M,sparsity));
        sparsity_actual = nnz(S)/(m*M); % often slightly under sparse
        S   = sqrt(1/(m*sparsity_actual))*S; % we may not have had exactly sparsity*m*n entries
        if tensormode
            fcn = @(A,MODE) ttm(A,S,MODE);
        else
            fcn = @(A) S*A;
        end
        
    case 'subsample'
        if isempty(weights)
            ind     = randperm( M, m );
            subsample   = @(X) X(ind,:);
            fcn     = @(A) sqrt(M/m)*subsample(A);
        else
            ind     = randsample( M, m, true, weights );
            subsample   = @(X) X(ind,:);
            fcn     = @(A) spdiags(sqrt(weights./m),0,m,m)*subsample(A);
        end
        
    case {'canonical','precondhd','uniformhd'}
        fcn = @(A) precondHD(A, m, true);
    case {'canonicalcompressed','precondhdcompressed','uniformhdcompressed'}
        fcn = @(A) precondHD(A, m, false);
        
    otherwise
        error('Invalid type of sketch requested');
end

if isempty(S) && nargout >= 2
    S   = fcn(eye(M));
end

end % end of main function

function Y = slowCountSketch( DX, targetRows )
% slow version of count sketch
    m   = length( targetRows );
    Y   = zeros(m, size(DX,2) );
    for j = 1:size(DX,1)
        i   = targetRows(j);
        Y(i,:) = Y(i,:) + DX(j,:);
    end
end

function Y = precondHD(X, m, undoF)
    if nargin < 3, undoF = true; end
% The preconditioning is done to the FIRST dimension of the input data.

    [T,N] = size(X);
    Y     = zeros(T,N);
    D     = spdiags( sign(randn(T,1)), 0, T, T );
    F     = @(x) dct( D*x );
    Ft    = @(x) D*idct(x);
    X     = F(X);
    if m < T
        for ColInd = 1:N
            randomIndices = randperm(T,m);
            Y(randomIndices,ColInd) = X(randomIndices,ColInd);
        end
        if undoF
            % This is what Leo was doing
            Y     = Ft(Y) * T/m; % unbiased estimator
        else
            % Keep it compressed
            Y   = sparse(Y)*T/m;
        end
    else
        Y   = X;
        if ~undoF
            error('If there is no sub-sampling, then not using "undoF" makes no sense');
        end
    end
end

function varargout = fastFD( X, l )
% fast frequent direction
% here l = 2*l_paper, l must be even
% Ghashami et al. 2015
[n,d] = size(X);
B = zeros(l, d);
if mod(l,2) ~= 0
    fprintf('\t input dimension augmented to be even! \n');
    l = l + 1;
end
if l > n
    fprintf('\t Subsample dimension given for Fast FD is large enough for no compression! \n')
    B = X;
    varargout{1} = B;
    return;
else
    pointer = 1;
    for ind = 1:n
        if pointer == l+1
            [~,Sigma,V] = svd(B,'econ');
            delta = Sigma(l/2,l/2)^2;
            core  = diag(Sigma(1:l/2,:)).^2 - delta * ones(l/2,1); 
            arr   = sqrt([core; zeros(l/2,1)]);
            matr = spdiags(arr,0,l,l);
%             B(find(abs(B)<=1e-16)) = 0;
            B = matr * V';
            pointer = l/2 + 1;
        else
            B(pointer,:) = X(ind,:);
            pointer = pointer + 1;
        end
    end
    varargout{1} = B;
    [~,~, varargout{2}] = svd(B);
    
end

end