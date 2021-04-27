function [C,lags] = findAutocorr_TimeBaseLine(data,varargin)
% [C,lags] = findAutocorr_TimeBaseLine(data,varargin)
%   Zhishen Huang 2019
%   Modified by Stephen Becker Dec 16 2019

prs = inputParser;
addParameter( prs, 'RowsMarker', []);
addParameter( prs, 'posOnly', true );
addParameter( prs, 'SCALEOPT', 'unbiased');
addParameter( prs, 'removeMean', true);
addParameter( prs, 'scale', 1);
addParameter( prs, 'interptype','spline');
addParameter( prs, 'iscanonical',false);    % false means time compression
addParameter( prs, 'm', []);
addParameter( prs, 'isCompressed', []);     % for time compression

parse(prs,varargin{:});
RowsMarker = prs.Results.RowsMarker;
posOnly    = prs.Results.posOnly;
removeMean = prs.Results.removeMean;
SCALEOPT   = prs.Results.SCALEOPT;
scale      = prs.Results.scale;
interptype = prs.Results.interptype;
iscanonical= prs.Results.iscanonical;
m          = prs.Results.m;
isCompressed    = prs.Results.isCompressed;

if ~iscanonical && ~isempty(RowsMarker)
    % Sub-sampling only in the time domain
    %     check RowsMarker are just 0 and 1
    if ~all(ismember(RowsMarker,[0;1]))
        error('The given Rows Markers are not in the correct format: {0,1} ~');
    end
    if issparse(RowsMarker)
        RowsMarker  = full(RowsMarker);
    end
    if isCompressed
        % This means that "data" is not really the full size
        nz      = data;
        data    = zeros( length(RowsMarker), size(nz,2) );
        data(RowsMarker,:) = nz;
    end
end

[T,N] = size( data );
if removeMean
    data = data - mean(data);
end

if issparse(data) || issparse( RowsMarker )
    f   = @(col) full(col);
else
    f   = @(col) col;
end

if ~iscanonical && ~isempty(RowsMarker)
    % We just sub-sampled the time domain
    
    [normalisers,lags] = xcorr(f(RowsMarker),round(T/2),SCALEOPT);
    normalisers        = normalisers(lags>=0);
    normalisers(normalisers<1/(T-1)) = inf;
    
    j = 1;
    [C,lags] = xcorr( f(data(:,j)), round(T/2), SCALEOPT );
    % Cmat     = zeros(N,length(C));
    % Cmat(j,:)= C;
    for j = 2:N
        [CC,lags] = xcorr( f(data(:,j)), round(T/2), SCALEOPT );
        C = C + CC;
        %     Cmat(j,:)= CC;
    end
    
    if nargin < 2 || isempty( posOnly ) || posOnly
        ind  = find( lags >= 0 ); % it is symmetric, so find just
        lags = lags(ind);
        C    = C(ind) ./ normalisers;
    end

elseif iscanonical && ~isempty(RowsMarker)
%     check RowsMarker are just 0 and 1
    if ~all(ismember(RowsMarker,[0;1]))
        error('The given Rows Markers are not in the correct format: {0,1} ~');
    end
    
    [normalisers,lags_normalisers] = xcorr(f(RowsMarker(:,1)),round(T/2),SCALEOPT);
    normalisers = normalisers(lags_normalisers>=0);
    normalisers(normalisers<1/(T-1)) = inf;    
    j = 1;
    [C,lags] = xcorr( f(data(:,j)), round(T/2), SCALEOPT );
    C = C(lags>=0)./ normalisers;
    
    for j = 2:N
        [CC,lags]          = xcorr( f(data(:,j)), round(T/2), SCALEOPT );
        [normalisers,lags_normalisers] = xcorr(f(RowsMarker(:,j)),round(T/2),SCALEOPT);
        normalisers        = normalisers(lags_normalisers>=0);
        normalisers(normalisers<1/(T-1)) = inf;
        CC = CC(lags>=0) ./ normalisers;
        C = C + CC;
    end
    
else
    % uniform sampling, process each particle in for-loop, inside which
    % another for-loop process each lag
    % this approach is 100% correct
    
    maxLag  = round(T/2);
    C       = zeros(maxLag+1,1);
    for col = 1:N
        C_tmp     = zeros(maxLag+1,1);
        counter   = zeros(maxLag+1,1);
        omegaCol  = find( data(:,col) );
        ynz       = nonzeros( data(:,col) );
        yy        = nan(T,1);
        yy(omegaCol) = ynz;
        for i = omegaCol'
            for k = 0:min(T-i,maxLag)
                if ~isnan(yy(k+i))
                    C_tmp(k+1)   = C_tmp(k+1) + yy(i)*yy(k+i);
                    counter(k+1) = counter(k+1) + 1;
                end
            end
        end
        counter(counter==0) = inf;
        C = C + C_tmp ./ counter;
    end
    
end

interpLoc = find(C==0);
filledLoc = setdiff(1:length(C),interpLoc);
if ~isempty(interpLoc)
    C(interpLoc) = interp1(filledLoc,C(filledLoc),interpLoc,interptype);
end

if iscanonical && ~isempty(m)
    C(1) = C(1) - (N-m)/(N*T)* norm(data,'fro')^2;
end

% C = C/N; % normalise particle dimension by compressed dimension size
C = C * scale; % normalise particle dimension by original dimension


end

