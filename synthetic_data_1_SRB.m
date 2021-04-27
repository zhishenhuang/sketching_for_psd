%{
Stephen's version of the synthetic data
Dec 16 2019

%}
clear all; close all; clc;
cd ~/Desktop/research/sketching/code_becker/
rng(0);
omega_max = 1/200 * pi;
period    = 2*pi/omega_max;

% MODE = 'time';
% MODE = 'particle';
MODE = 'comp';

Time = round(100*period);
N    = 1e4;

% Amplitude   = 5e1/N;  % wlog, only scale sine, not noise
Amplitude   = 100*1.5/sqrt(N); % 100* makes it easy, and we have an easy story
omega       = .5*omega_max;
% omega       = 5*omega_max;

freqBdd = 0.1;
dt      = freqBdd * pi/omega_max;
time    = (1:dt:Time)';
T       = length(time);
X       = Amplitude*sin( omega*time + 2*pi*rand(1,N) ) + randn(T,N); % random phase sine, and white noise
% If we want particle compression to not work well...
%   suppose there are a few particles that have a special mode at omega2
omega2              = 5*omega;
switch N
    case 1e3
        Amplitude2          = 20*Amplitude;
    case 1e4
        Amplitude2          = 80*Amplitude;
end

switch MODE
    case 'particle'
        specialParticles    = randperm( N, 2 );
        X(:,specialParticles) = X(:,specialParticles) + ...
            Amplitude2*sin( omega2*time + 2*pi*rand(1,length(specialParticles)) );

    case 'time'
        pulseAmp    = 100*Amp;
        PulseCenter = randi([100,T-100],2);
        pulse = @(t,t0,A,L) A*sin(pi/L*(t-t0)); 
        halfPulseWidth = 25;
        Pulse1 = [PulseCenter(1)-halfPulseWidth:PulseCenter(1)+halfPulseWidth]';
        Pulse2 = [PulseCenter(2)-halfPulseWidth:PulseCenter(2)+halfPulseWidth]';
        X(Pulse1,:) = X(Pulse1,:) +...
            pulse(Pulse1,PulseCenter(1)-halfPulseWidth,pulseAmp,2*halfPulseWidth);
        X(Pulse2,:) = X(Pulse2,:) +...
            pulse(Pulse2,PulseCenter(2)-halfPulseWidth,pulseAmp,2*halfPulseWidth);
        figure(1); clf;
        plot(X(:,2));
        
    case 'comp'
        specialParticles    = randperm( N, 2 );
        X(:,specialParticles) = X(:,specialParticles) + ...
            Amplitude2*sin( omega2*time + 2*pi*rand(1,length(specialParticles)) );
        pulseAmp    = 10*Amplitude;
        NumPulse    = 2;
        PulseCenter = randi([100,T-100],NumPulse,1);
        pulse  = @(t,t0,A,L) A*sin(pi/L*(t-t0)); % cosine can also work
        halfPulseWidth = 25;
        fprintf('dt = %.5f\n',dt);
        fprintf('T = 2*pi/omega = %.5f\n omega = %.5f\n grid_per_period = %.5f\n',...
            2*pi/omega,omega,2*pi/omega/dt);
        Pulses = zeros(length(time),1);
        for pulseInd = 1:NumPulse
            PulsesLoc = [PulseCenter(pulseInd)-halfPulseWidth:PulseCenter(pulseInd)+halfPulseWidth]';
            Pulses(PulsesLoc) = Pulses(PulsesLoc) +...
                pulse(PulsesLoc,PulseCenter(pulseInd)-halfPulseWidth,pulseAmp,2*halfPulseWidth);
        end
        fprintf('grid_per_halfpulse = %.5f\n',halfPulseWidth);
        fprintf('T_pulse width = %.5f\n', 2*halfPulseWidth*dt );
        fprintf('pulse omega = %.5f\n',  pi/(halfPulseWidth*dt));
        X = X + Pulses;    
        figure('units','normalized','outerposition',[0 0 1 1/2]); clf;
        subplot(1,2,1);
        plot(X(:,2));
%         xlim([900 1100])
        title('common particle')
        subplot(1,2,2);
        plot(X(:,specialParticles(1)));
%         xlim([900 1100])
        title('special particle')
end

% Data Preparation: Remove Mean; Ground Truth Autocorrelation visualization
% removeMean  = true;
removeMean  = false;
if removeMean
    X    = X - mean(X);
end
posOnly     = true; 
% SCALEOPT    = 'coeff';
SCALEOPT    = 'unbiased';
% AUTOCORR    = @(X,varargin) findAutocorr(X,'posOnly',posOnly,...
%     'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,varargin{:} );
AUTOCORR    = @(X,varargin) findAutocorr(X,'posOnly',posOnly,...
    'removeMean',removeMean,'SCALEOPT',SCALEOPT,varargin{:} );
pos_frequencies     = @( f ) f( 1:round(length(f)/2)+1 );
psd_from_autocorr   = @( cc ) pos_frequencies( abs(fft(cc)) );

% First, data for the true original signal
[cc_orig,lags]  = AUTOCORR(X,'scale',1/N);
psd_orig        = psd_from_autocorr(cc_orig);

figure(2); clf;
subplot(2,1,1)
plot(lags,cc_orig);     
title('Autocorrelation ground truth')
set(gca,'FontSize',23)
% axis tight
% xlim([0,10]);

subplot(2,1,2)
plot( psd_orig  )
title('Power Density ground truth'); axis tight
set(gca,'FontSize',23)

% and plot the special region with the small spike
if strcmpi(MODE,'particle')|| strcmpi(MODE,'comp')
    importantFrequencies    = 120:140;
    hold all
    plot( importantFrequencies, psd_orig(importantFrequencies),'r-'  )
end

%% plot compression ratio v.s. accuracy
% sketch.m not same as ~/Repos/randomized-algorithm-class/Code

% NGrid   = logspace(2,log10(N), 3 );
NGrid   = N;
compressRatio_Grid = [.003,.005,.007,.01,0.03,.05,.07,.1];
repeat  = 1e2;

mem0    = whos('X').bytes;

% ALGONAMES = {'Haar','Gaussian','FJLT_Hadamard','UniformHD',...
%     'Particle Baseline', 'Time Baseline','Naive Uniform'};
ALGONAMES = {'Haar','Gaussian','FJLT_Hadamard',...
    'Particle Baseline', 'Time Baseline','Naive Uniform'};
ALGOLIST = [1:6];

% ALGONAMES   = {'Haar','Gaussian','FJLT','UniformHD','UniformHD v2',...
%     'UniformHD v3','Particle Baseline', 'Time Baseline','Naive uniform',...
%     'Naive uniform iid'};
% ALGOLIST    = [1:3,5,7,8,9];
% [err_psd_l2, err_psd_lInf, err_psd_W, err_cc_l2, err_psd_zoom, storage ] = deal( ...
%     zeros(length(NGrid), length(compressRatio_Grid), length(ALGONAMES), repeat ) );

[err_psd_l1,err_psd_l2, err_psd_lInf,err_psd_lInfElem, err_psd_W, err_psd_zoom,...
    err_cc_l1, err_cc_l2, err_cc_lInf, err_cc_lInfElem, storage ] = deal( ...
zeros(length(NGrid), length(compressRatio_Grid), length(ALGONAMES), repeat ) );

errInfElem = @(x,x0) max( abs(x(x0~=0)- x0(x0~=0))./abs(x0(x0~=0)) );
errL1      = @(x,xTrue) norm(x-xTrue,1)/norm(xTrue,1);

%% Computation
for NInd = 1:length(NGrid)
    fprintf('~~~ NGrid Point %u out of %u\n',NInd,length(NGrid));
    NN = round(NGrid(NInd));
    XSlice = X(:,1:NN);
    for CompInd = 1:length(compressRatio_Grid)
        compressRatio = compressRatio_Grid(CompInd);
        m = round(compressRatio*NN);  
        fprintf('~~~~ Compression ratio %d of %d (ratio is %g)\n', CompInd, length(compressRatio_Grid), compressRatio );
        
        for rep = 1:repeat
            fprintf('~~~~~ Rep %u\n',rep);
            for ALGO    = ALGOLIST
                NAME    = lower(ALGONAMES{ALGO});
%                 if rep == 5 && strcmpi(NAME, 'time baseline'), keyboard; end
                switch ALGO
                    case {1,2,3}
                        sk  = sketch( m, NN, NAME );
                        XX  = sk(XSlice')';
                        cc  = AUTOCORR( XX ,'scale',1/NN);
%                     case 4  
%                         % Leo's version of unformHD
%                         sk  = sketch( m, NN, 'uniformHD');
%                         XX  = sk(XSlice')';
%                         cc  = AUTOCORR( XX, 'iscanonical', true,'m',m );
%                     case 5
%                         % Stephen's version 1 of unformHD
%                         sk  = sketch( m, NN, 'UniformHDcompressed');
%                         XX  = sk(XSlice')';
%                         cc  = AUTOCORR( XX,'scale',1/NN ,'iscanonical', true,'m',m );
%                     case 6
%                         % Stephen's version 2 of unformHD
%                         sk  = sketch( NN, NN, 'uniformHD'); % only mixing, no subsampling
%                         XX  = sk(XSlice')';
%                         % Now... 
%                         % Naive uniform, Leo's version
%                         [XX,RowsMarker_uni] = baseline(XX,compressRatio, 'all');
%                         cc = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_uni,...
%                             'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
%                             'iscanonical',true,'m',m);
                    case 4
                        % Particle (baseline), Leo's version
                        XX  = sqrt(NN/m)*baseline(XSlice, m, '2');
                        cc  = AUTOCORR( XX ,'scale',1/NN );
                    case 5
                        % Time (baseline), Leo's version
                        try
                            [XX,RowsMarker_sub1]  = baseline(XSlice, round(compressRatio*T),'1-alternative');
                            RowsMarker_sub1       = logical( RowsMarker_sub1 );
                            XX  = XX( RowsMarker_sub1, :); % "isCompressed" is now True
                            RowsMarker_sub1     = sparse( RowsMarker_sub1 );
                            cc = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_sub1,'isCompressed',true,...
                                'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/NN);
                        catch ME
                            [storage( NInd, CompInd, ALGO, rep),...
                            err_cc_l1(  NInd, CompInd, ALGO, rep ),...
                            err_cc_l2(  NInd, CompInd, ALGO, rep ),...
                            err_cc_lInfElem( NInd, CompInd, ALGO, rep ),...
                            err_psd_l1( NInd, CompInd, ALGO, rep ),...
                            err_psd_l2( NInd, CompInd, ALGO, rep ),...
                            err_psd_lInfElem( NInd, CompInd, ALGO, rep ) ]= deal(NaN);                            
%                             err_cc_lInf( NInd, CompInd, ALGO, rep ),...
%                             err_psd_lInf( NInd, CompInd, ALGO, rep ),...
%                             err_psd_W(  NInd, CompInd, ALGO, rep ),...
%                             err_psd_zoom( NInd, CompInd, ALGO, rep )] = deal(NaN);
                            fprintf('Method is %18s, Rep %u failed!!! \n', NAME, rep );
                            continue;
                        end
                    case 6
                        % Naive uniform, Leo's version
                        [XX,RowsMarker_uni] = baseline(XSlice,compressRatio, 'all');
                        cc = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_uni,...
                            'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/NN,...
                            'iscanonical',true,'m',m);
%                     case 10
%                         % Naive uniform, iid, so only control expected # of
%                         % nonzeros
%                         [XX,RowsMarker_uni] = baseline(XSlice,compressRatio, 'all-iid');
%                         cc = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_uni,...
%                             'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
%                             'iscanonical',true,'m',m);
                end
                mem    = whos('XX').bytes;
                if ALGO==6, mem = mem + whos('RowsMarker_sub1').bytes; end
                fprintf('Method is %18s, X is compressed by %g\n', NAME, mem0/mem );
                if issparse(X), fprintf('\tthough nnz is reduced by %g\n', numel(X)/nnz(XX) ); end
                
                storage( NInd, CompInd, ALGO, rep) = mem/mem0;
                
                psd = psd_from_autocorr( cc );
                err_cc_l1(  NInd, CompInd, ALGO, rep )   = errL1(cc, cc_orig);
                err_cc_l2(  NInd, CompInd, ALGO, rep )   = norm(cc - cc_orig)/norm(cc_orig);
%                 err_cc_lInf( NInd, CompInd, ALGO, rep )  = norm(cc - cc_orig, Inf);
                err_cc_lInfElem( NInd, CompInd, ALGO, rep ) = errInfElem(cc,cc_orig);

                err_psd_l1( NInd, CompInd, ALGO, rep )   = errL1(psd, psd_orig);
                err_psd_l2( NInd, CompInd, ALGO, rep )   = norm(psd - psd_orig)/norm(psd_orig);
%                 err_psd_lInf( NInd, CompInd, ALGO, rep ) = norm(psd - psd_orig, Inf);
%                 err_psd_W(  NInd, CompInd, ALGO, rep )   = wasserstein_distance(psd, psd_orig, 'p', 1);
                err_psd_lInfElem( NInd, CompInd, ALGO, rep ) = errInfElem(psd,psd_orig);
                
                if strcmpi(MODE,'particle') || strcmpi(MODE,'comp')
                    errZoomFcn  = @(psd) norm( psd(importantFrequencies) - psd_orig(importantFrequencies) )/norm( psd_orig(importantFrequencies) );
                    err_psd_zoom( NInd, CompInd, ALGO, rep ) = errZoomFcn( psd );
                end
                
%                 % Optional: plot
%                 if length(compressRatio_Grid)==1 && length(NGrid)==1 && repeat==1
%                     figure(1); clf;
%                     plot( psd_orig , 'DisplayName','Reference','linewidth',2 )
%                     % and plot the special region with the small spike
%                     importantFrequencies    = 120:140;
%                     hold all
%                     plot( importantFrequencies, psd_orig(importantFrequencies),'r-'  )
%                     plot( psd, '--','linewidth',2,'DisplayName', NAME );
%                     legend
%                     pause
%                 end
                
            end
        end
    end
end

%
clear X XSlice XX
% save('Pulses_Res.mat');
save('Comp_Res_July.mat');

%% Plot error

% == Choice 1: which error metric to show (pick one)  
% errTensor   = err_psd_l2; str = 'relative l2 error for PSD';
% errTensor   = err_psd_lInf; str = 'l_{inf} error';
% errTensor   = err_psd_zoom;   str = 'relative l2 error in mid-frequencies';
errTensor   = err_cc_l2; str = 'relative l2 error for autocorrelation';
% errTensor   = err_psd_W; str = 'Wasserstein p=1 distance';

% == Choice 2: take mean/median/min ?
% Indices are N x compression x ALGO x rep, so look at 4th index
% err         = mean(   errTensor, 4 ); type = 'mean';
err         = median( errTensor, 4 ); type = 'median';
% err         = min(    errTensor, [], 4 ); type = 'min';
compression = mean(   storage, 4 );       % average compression

% figure(1); clf;     % at a fixed compression
% CompInd     = 1;
% compressRatio = compressRatio_Grid(CompInd);
% holdMarker( 'reset' );
% for ALGO    = 1:length(ALGONAMES)
%     plot( NGrid, err(:,CompInd,ALGO), 'linewidth',2,'DisplayName',ALGONAMES{ALGO} );
%     hold all
% end
% holdMarker( 0 );
% legend
% xlabel('NGrid');
% title(sprintf('Error vs N, for approximate compression factor %g',compressRatio) );

%%
figure(1); clf;     % at a fixed N
holdMarker( 'reset' );
NInd    = 1;
NN = round(NGrid(NInd));

for ALGO    = [1:3,5,7,8,9]
    % Wrong way (uses wrong compression ratio for x points)
%     plot( compressRatio_Grid,err(NInd,:,ALGO), 'o--', ...
%         'linewidth',2,'DisplayName',ALGONAMES{ALGO},'markersize',16  );
    
    % Right way:
    semilogx( compression(NInd,:,ALGO),err(NInd,:,ALGO), 'o--',...
        'linewidth',2,'DisplayName',ALGONAMES{ALGO},'markersize',16 );
    hold all
end
holdMarker( 0 );
set(gca,'FontSize',20);
ylim([0,1]);
legend
xlabel('Compression');
ylabel(sprintf('Error (%s, %s over %d trials)',str,type,repeat));
title(sprintf('Error vs compression, for N =  %g',NN) );

%% Save
% saveFigure is old, use savefig

% savefig(gcf,'N1e4_twoSignals_Dec16_relError','compact')
% export_fig N1e4_twoSignals_Dec16_relError -pdf -transparent
