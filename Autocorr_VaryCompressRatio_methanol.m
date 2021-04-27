%{ Autocorrelation Evaluation
%
%}

%% Load data
clear all; close all; clc;
addpath ~/Desktop/research/sketching/MDSpectralAnalysis-master/Data/Methanol_Velocity_Data/
% addpath /home/zhhu7269/autocorr/data/Methanol_Velocity_Data/
rng(2343);
[t,Vx,Vy,Vz] = loadMetOHVel();

[T,N] = size(Vx');
trueX     = Vx(:,:)'; % trueX is TIME BY SPACE

trueX = trueX(ceil(T/2):T,:); % burnout for equilibrium, June 9, 2020
[T,~] = size(trueX);

% removeMean  = true;
removeMean  = false;
X    = trueX;
if removeMean
    X    = X - mean(X);
end
mem0    = whos('X').bytes;

posOnly     = true;
% SCALEOPT    = 'coeff';
SCALEOPT    = 'unbiased';
AUTOCORR    = @(X,varargin) findAutocorr(X,'posOnly',posOnly,...
    'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,varargin{:} );
pos_frequencies     = @( f ) f( 1:round(length(f)/2)+1 );
psd_from_autocorr   = @( cc ) pos_frequencies( abs(fft(cc)) );

[cc_orig,lags] = AUTOCORR(X);
psd_orig = psd_from_autocorr( cc_orig );

%
importantFrequencies = 1:5e2;
% lowLags  = 1:33;
% midLags  = 34:67;
% highLags = 68:101;
lowLags  = 1:200;
midLags  = 200:400;
highLags = 400:600;
% importantFrequencies = round(length(psd_orig)*2/3):length(psd_orig);

%
importantLags = 1:101;
figure(1); clf;
subplot(2,1,1)
plot(lags(importantLags),cc_orig(importantLags));     
title('Autocorr ground truth')
set(gca,'FontSize',23)
% axis tight
% xlim([0,10]);
% set(gca,'yscale','log')

subplot(2,1,2)
plot( psd_orig  )
title('PSD ground truth'); axis tight
set(gca,'FontSize',23)

% and plot the low frequency region
hold all
plot( importantFrequencies, psd_orig(importantFrequencies),'r-'  )

%% Data Preparation: Remove Mean; Ground Truth Autocorrelation visualization

compressRatio_Grid = 10.^(linspace(-2,-1,10));
% Tol    = 5e-2;
repeat = 1e3;
% ALGONAMES = {'Haar','Gaussian','FJLT-Hadamard','UniformHD',...
%     'Particle Baseline', 'Time Baseline','Naive Uniform'};
% ALGOLIST = [1:7];

ALGONAMES = {'Haar','Gaussian','FJLT-Hadamard',...
    'Particle Baseline', 'Time Baseline','Naive Uniform'};
ALGOLIST = [1:6];

[err_psd_l2, err_psd_l1, err_psd_lInf,err_psd_lInfElem, err_psd_W, err_psd_zoom,...
    err_cc_l1,err_cc_l2,err_cc_lInf,err_cc_lInfElem,err_cc_low,err_cc_mid,err_cc_high,...
    storage ] =...
    deal( zeros(length(compressRatio_Grid), length(ALGONAMES), repeat) ) ;
errInfElem = @(x,x0) max( abs(x(x0~=0)- x0(x0~=0))./abs(x0(x0~=0)) );

%% Computation, and plot compression ratio v.s. accuracy
for CompInd = 1:length(compressRatio_Grid)
% for ind = 1:length(compressRatio_Grid)
    fprintf('   Current Compress Ratio Grid Point %u out of %u\n',CompInd,length(compressRatio_Grid));
    compressRatio = compressRatio_Grid(CompInd);
    m             = round(compressRatio*N);
    rep = 1;
    while rep<=repeat
        fprintf('   ~~~Repeat %u at Compress Ratio Grid Point %u~~~\n',rep,CompInd);
        
        for ALGO = ALGOLIST
            NAME = ALGONAMES{ALGO};
            switch NAME
                case 'Haar'
                    HAAR    = sketch( m, N, 'haar' );
                    XX      = HAAR(X')';
                    cc      = AUTOCORR(XX);
                case 'Gaussian'
                    GAUSSIAN= sketch( m, N, 'gaussian' );
                    XX      = GAUSSIAN(X')';
                    cc      = AUTOCORR(XX);
                case 'FJLT-Hadamard'
                    FJLT    = sketch( m, N, 'fjlt_hadamard' );
                    XX      = FJLT(X')';
                    cc      = AUTOCORR(XX);
                case 'UniformHD'
%                     PRECONDHD    = sketch( m, N, 'canonical' );
                    PRECONDHD    = sketch( m, N, 'UniformHDcompressed' );
                    XX  = PRECONDHD(X')';
                    cc  = AUTOCORR( XX, 'iscanonical', true,'m',m );
                case 'Particle Baseline'
                    XX  = sqrt(N/m)*baseline(X, m, '2');
                    cc  = AUTOCORR( XX );
%                     XX      = baseline(X, m, '2');
%                     cc      = findAutocorr(XX,'posOnly',posOnly,...
%                         'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/m);
                case 'Naive Uniform'
                    [XX,RowsMarker_uni] = baseline(X,compressRatio, 'all');
                    cc       = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_uni,...
                        'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
                        'iscanonical',true,'m',m);
            %         X_UNI        = baseline(X',m, 'all-old')';
            %         cc_uni_old       = findAutocorr_TimeBaseLine(X_UNI,'posOnly',posOnly,...
            %             'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
            %             'iscanonical',true,'m',m);
                case 'Time Baseline'
%                     [XX,RowsMarker_sub1]  = baseline(X, round(compressRatio*T),'1-alternative');
%             %         cc_sub1_old     = findAutocorr_TimeBaseLine(X_SUB1,...
%             %             'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
%                     cc     = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_sub1,...
%                         'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
                    [XX,RowsMarker_sub1]  = baseline(X, round(compressRatio*T),'1-alternative');
                    RowsMarker_sub1       = logical( RowsMarker_sub1 );
                    XX  = XX( RowsMarker_sub1, :); % "isCompressed" is now True
                    RowsMarker_sub1     = sparse( RowsMarker_sub1 );
                    cc = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_sub1,'isCompressed',true,...
                        'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
            end
            mem    = whos('XX').bytes;
            if ALGO==6, mem = mem + whos('RowsMarker_sub1').bytes; end
            fprintf('Method is %18s, X is compressed by %g\n', NAME, mem0/mem );
            if issparse(X), fprintf('\tthough nnz is reduced by %g\n', numel(X)/nnz(XX) ); end
            
            storage( CompInd, ALGO, rep) = mem/mem0;
            
            psd = psd_from_autocorr( cc );
            
            err_cc_l1( CompInd, ALGO, rep )        = norm( cc - cc_orig,1 )/norm(cc_orig,1);
            err_cc_l2( CompInd, ALGO, rep )        = norm( cc - cc_orig )/norm(cc_orig);
            err_cc_low( CompInd, ALGO, rep )       = norm( cc(lowLags) - cc_orig(lowLags) )/norm(cc_orig(lowLags));
            err_cc_mid( CompInd, ALGO, rep )       = norm( cc(midLags) - cc_orig(midLags) )/norm(cc_orig(midLags));
            err_cc_high( CompInd, ALGO, rep )      = norm( cc(highLags)- cc_orig(highLags))/norm(cc_orig(highLags));
            err_cc_lInf( CompInd, ALGO, rep )      = norm( cc - cc_orig, Inf );
            err_cc_lInfElem( CompInd, ALGO, rep )  = errInfElem( cc,cc_orig );
            
            err_psd_l1( CompInd, ALGO, rep )       = norm( psd - psd_orig ,1  )/norm(psd_orig ,1);
            err_psd_l2( CompInd, ALGO, rep )       = norm( psd - psd_orig  )/norm(psd_orig);
            err_psd_lInf( CompInd, ALGO, rep )     = norm( psd - psd_orig, Inf );
            err_psd_lInfElem( CompInd, ALGO, rep ) = errInfElem(psd,psd_orig);
%             err_psd_W( CompInd, ALGO, rep )        = wasserstein_distance( psd, psd_orig, 'p', 1 );
            errZoomFcn  = @(psd) norm( psd(importantFrequencies) - psd_orig(importantFrequencies) )/norm( psd_orig(importantFrequencies) );
            err_psd_zoom( CompInd, ALGO, rep )     = errZoomFcn( psd );
            
        end
     
        rep = rep + 1;
    end
end

%%
clear t Vx Vy Vz trueX X XX
save('MD_Res_Jun.mat');
