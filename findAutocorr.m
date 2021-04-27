function [C,lags] = findAutocorr( data, varargin)
% [C,lags] = findAutocorr( data, [posOnly], [removeMean], [ScaleOpt] )
%   if posOnly = true (default), only shows positive lags
%   Stephen Becker, Nov 5 2018, based on old 2016, 2017 scripts
%
%   "data" should be of size T x N, where T is the time-axis (for
%   correlation) and N is # of observations (eg., spatial) which will
%   be summed over.
%
%   2019: TODO, can I avoid for loop? No, xcorr( matrix ) does wrong
%       thing.
%   12/16/2019, allows for data to be sparse

prs = inputParser;
addParameter( prs, 'posOnly', true );
addParameter( prs, 'SCALEOPT', 'unbiased');
addParameter( prs, 'removeMean', true);
addParameter( prs, 'scale', 1);
addParameter( prs, 'iscanonical',false);
addParameter( prs, 'm', []);

parse(prs,varargin{:});
posOnly    = prs.Results.posOnly;
removeMean = prs.Results.removeMean;
SCALEOPT   = prs.Results.SCALEOPT;
scale      = prs.Results.scale;
iscanonical= prs.Results.iscanonical;
m          = prs.Results.m;

% if nargin < 4 || isempty( SCALEOPT )
%     SCALEOPT = 'coeff';
%     %SCALEOPT = 'unbiased';
% end
% 
% if nargin < 3 || isempty( removeMean )
%     removeMean  = true;
% end

[T,N] = size( data );
if removeMean
    data = data - mean(data); % helps, so that correlation --> 0 as t --> 0
    % mean(data) is 1 x Phi, so we are removing the same spatial
    %   average from all the data points
end
   
j = 1;
[C,lags] = xcorr( full(data(:,j)), round(T/2), SCALEOPT );
% Cmat     = zeros(N,length(C));
% Cmat(j,:)= C;
for j = 2:N
    [CC,lags] = xcorr( full(data(:,j)), round(T/2), SCALEOPT );
    C = C + CC;
%     Cmat(j,:)= CC;
end

% C = C/N; % normalise particle dimension by compressed dimension size
if nargin < 2 || isempty( posOnly ) || posOnly
    ind = find( lags >= 0 ); % it is symmetric, so find just
    lags     = lags(ind);
    C        = C(ind);
end

if iscanonical && ~isempty(m)
   C(1) = C(1) - (N-m)/(N*T)* norm(data,'fro')^2; 
end

C = C * scale; % normalise particle dimension by original dimension

end