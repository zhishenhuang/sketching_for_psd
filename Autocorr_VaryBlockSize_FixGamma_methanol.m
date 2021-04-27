%{
Stephen's version of the synthetic data
Dec 16 2019

%}
clear all; close all; clc;
rng(0);
addpath ~/Desktop/research/sketching/MDSpectralAnalysis-master/Data/Methanol_Velocity_Data/

% DATASOURCE = 'synthetic';
DATASOURCE = 'MD';

%%

switch DATASOURCE
    case 'synthetic'
        omega_max = 1/200 * pi;
        period = 2*pi/omega_max;
        Time= round(10^(3)*period);
        N   = 1e3;
        
        % Amplitude   = 5e1/N;  % wlog, only scale sine, not noise
        Amplitude   = 1e2*1.5/sqrt(N); % 100* makes it easy, and we have an easy story
        omega       = .5*omega_max;
        freqBdd = 0.1;
        dt      = freqBdd * pi/omega_max;
        time    = (1:dt:Time)';
        T       = length(time);
        % If we want particle compression to not work well...
        %   suppose there are a few particles that have a special mode at omega2
        omega2  = 5*omega;        
        Amplitude2          = 3*Amplitude;
        Vx = zeros(T,N);
        species1Perc = .9;
        Species1 = randperm(N,round(N*species1Perc));
        Species2 = setdiff(1:1:N,Species1);
        NoiseAmp1     = 0.1 * Amplitude;
        NoiseAmp2     = 0.1 * Amplitude2;
        Vx(:,Species1) = Amplitude *sin(omega *time+2*pi*rand(1,length(Species1))) + NoiseAmp1*randn(T,length(Species1));
        Vx(:,Species2) = Amplitude2*sin(omega2*time+2*pi*rand(1,length(Species2))) + NoiseAmp2*randn(T,length(Species2));
        
        Grid = 10.^(linspace(2.5,log10(T),6));
        compressRatio = 0.1;
        repeat  = 500;
        
%         WINDOW = false;
        WINDOW = 'bartlett';
        MaxLag = 40;
        
%% plot compression ratio v.s. accuracy
% sketch.m not same as ~/Repos/randomized-algorithm-class/Code
    case 'MD'
        [t,VxFull,VyFull,VzFull] = loadMetOHVel();
        [T,N] = size(VxFull');
        Vx = VxFull(:,ceil(T/2):T); % burnout for equilibrium, June 9, 2020
        [~,T] = size(Vx);
        Grid = 10.^(linspace(2,4,6));
%         Grid = 10.^(linspace(1.8,4,6));
        compressRatio = .1;
        repeat  = 20;
        
        WINDOW = 'bartlett';
        MaxLag = 15;
%         MaxLag = 3e4;
end


ALGONAMES = {'Haar','Gaussian','FJLT_Hadamard',...
    'Particle Baseline', 'Time Baseline','Naive Uniform'};
% ALGOLIST = [2];
ALGOLIST = [1:6];

[err_psd_l2, err_psd_l1,err_psd_lInfElem,...
err_cc_l1, err_cc_l2, err_cc_lInf, err_cc_lInfElem, storage ] = deal( ...
    zeros(length(Grid), length(ALGONAMES), repeat ) );

errInfElem = @(x,x0) max( abs(x(x0~=0)- x0(x0~=0))./abs(x0(x0~=0)) );
% plotFlag = true;
plotFlag = false;
skrep = 1; % repeat of sketching within one window

%% Main Body:
% for gridInd = [4:6]
for gridInd = 1:length(Grid)
%     gridInd = 5;
    m = round(compressRatio*N);
    TMax = round(Grid(gridInd));
    fprintf('~~~ Grid Point %u out of %u\n',gridInd,length(Grid));
%% Computation    
    for rep = 1:repeat          
          switch DATASOURCE
              case 'synthetic'
                  spurious = T - TMax;
                  if repeat >= spurious && rep <= spurious
                      shift = rep-1;
                  else
                      shift = randi([0,spurious]);
                  end
                  X = Vx(1+shift:TMax+shift,:);
              case 'MD'
                  spurious = T - TMax;
                  if repeat >= spurious && rep <= spurious
                      shift = rep;
                  else
                      shift = randi([0,spurious]);
                  end
                  X = Vx(:,1+shift:TMax+shift)';
          end
          mem0    = whos('X').bytes;
          
          %% Data Preparation: Remove Mean; Ground Truth Autocorrelation visualization
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
%           psd_orig        = psd_from_autocorr(cc_orig);
          
          if plotFlag
              %%
              %         [cc_orig,lags]  = AUTOCORR(Vx,'scale',1/N);
              %         psd_orig        = psd_from_autocorr(cc_orig);
              figure(1); clf;
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
              
          end
    
        fprintf('~~~~~ Rep %u\n',rep);
        for ALGO    = ALGOLIST
            NAME    = lower(ALGONAMES{ALGO});
            switch ALGO
                case {1,2,3}
                    switch WINDOW
                        case false
                            sk  = sketch( m, N, NAME );
                            XX  = sk(X')';
                            cc  = AUTOCORR( XX ,'scale',1/N);
                            for skrepInd = 2:skrep
                                sk  = sketch( m, N, NAME );
                                XX  = sk(X')';
                                cc_tmp  = AUTOCORR( XX ,'scale',1/N);
                                cc = cc + cc_tmp;
                            end
                            cc = cc/skrep;
                        case 'bartlett'
%                             batchSize = min([max(round(sqrt(T)/1.8),2*MaxLag),TMax]);
                            batchSize = min([max(round(sqrt(T)),2*MaxLag),TMax]);
                            batchNum  = floor(TMax/batchSize);
                            sk  = sketch( m, N, NAME );
                            XSlice = X(1:batchSize,:);
                            XX  = sk(XSlice')';
                            cc  = AUTOCORR( XX ,'scale',1/N);
                            for skrepInd = 2:skrep
                                sk  = sketch( m, N, NAME );
                                XX  = sk(XSlice')';
                                cc_tmp = AUTOCORR( XX ,'scale',1/N);
                                cc  = cc + cc_tmp;
                            end
                            cc = cc/skrep;
%                             for batchInd = 2:batchNum
%                                 sk = sketch(m,N,NAME);
%                                 XX = sk(X((batchInd-1)*batchSize+1:min([batchInd*batchSize,TMax]),:)')';
%                                 cc_tmp = AUTOCORR(XX,'scale',1/N);
%                                 cc(1:length(cc_tmp)) = (cc(1:length(cc_tmp))*(batchInd-1) + cc_tmp)/batchInd;
%                             end
                            for batchInd = 2:batchNum
                                sk = sketch(m,N,NAME);
                                XSlice = X((batchInd-1)*batchSize+1:min([batchInd*batchSize,TMax]),:);
                                XX = sk(XSlice')';
                                cc_batch = AUTOCORR(XX,'scale',1/N);
                                for skrepInd = 2:skrep
                                    sk  = sketch( m, N, NAME );
                                    XX  = sk(XSlice')';
                                    cc_tmp = AUTOCORR( XX ,'scale',1/N);
                                    cc_batch  = cc_batch + cc_tmp;
                                end
                                cc_batch = cc_batch/skrep;
                                cc = cc + cc_batch;
                            end
                            cc = cc/batchNum;
                            
%                             psd = psd_from_autocorr(cc);
                        case 'welch'
                            batchSize = min([max(round(sqrt(T)),MaxLag),TMax]);
                            Shift  = max(round(batchSize/50),1);
                            assert(Shift >= 1);
                            sk  = sketch( m, N, NAME );
                            XX  = sk(X(1:batchSize,:)')';
                            cc  = AUTOCORR( XX ,'scale',1/N);
                            batchEnd = batchSize+Shift;
                            batchBegin = 1+Shift;
                            batchNum = 1;
                            while batchEnd < TMax
                                sk = sketch(m,N,NAME);
                                batchComp = batchBegin:batchEnd;
                                XX = sk(X(batchComp,:)')';
                                cc = cc + AUTOCORR(XX,'scale',1/N);
                                batchBegin = batchBegin + Shift;
                                batchEnd   = batchEnd + Shift;
                                batchNum   = batchNum + 1;
                            end
                            cc = cc/batchNum;
%                             while batchEnd < TMax
%                                 sk = sketch(m,N,NAME);
%                                 batchComp = batchBegin:batchEnd;
%                                 XX = sk(X(batchComp,:)')';
%                                 cc_tmp = AUTOCORR(XX,'scale',1/N);
%                                 cc(1:length(cc_tmp)) = (cc(1:length(cc_tmp))*(batchInd-1) + cc_tmp)/batchInd;
%                                 batchBegin = batchBegin + Shift;
%                                 batchEnd = min(batchEnd + Shift,TMax);
%                                 batchInd = batchInd + 1;
%                             end
                    end
                case 4
                    % Stephen's version 1 of unformHD
                    sk  = sketch( m, N, 'UniformHDcompressed');
                    XX  = sk(X')';
                    cc  = AUTOCORR( XX,'scale',1/N ,'iscanonical', true,'m',m );
                case 5
                    % Particle (baseline), Leo's version
                    XX  = sqrt(N/m)*baseline(X, m, '2');
                    cc  = AUTOCORR( XX ,'scale',1/N );
                    for skrepInd = 2:skrep
                        XX  = sqrt(N/m)*baseline(X, m, '2');
                        cc_tmp  = AUTOCORR( XX ,'scale',1/N );
                        cc = cc + cc_tmp;
                    end
                    cc = cc/skrep;
                case 6
                    % Time (baseline), Leo's version
%                     try
                        [XX,RowsMarker_sub1]  = baseline(X, round(compressRatio*TMax),'1-alternative');
                        RowsMarker_sub1       = logical( RowsMarker_sub1 );
                        XX  = XX( RowsMarker_sub1, :); % "isCompressed" is now True
                        RowsMarker_sub1     = sparse( RowsMarker_sub1 );
                        cc = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_sub1,'isCompressed',true,...
                            'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
                        
                        for skrepInd = 2:skrep
                            [XX,RowsMarker_sub1]  = baseline(X, round(compressRatio*TMax),'1-alternative');
                            RowsMarker_sub1       = logical( RowsMarker_sub1 );
                            XX  = XX( RowsMarker_sub1, :); % "isCompressed" is now True
                            RowsMarker_sub1     = sparse( RowsMarker_sub1 );
                            cc_tmp = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_sub1,'isCompressed',true,...
                                'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
                            cc = cc+cc_tmp;
                        end
                        cc = cc/skrep;
                        
%                     catch ME
%                         [storage( gridInd, ALGO, rep),...
%                             err_cc_l2(  gridInd, ALGO, rep ),...
%                             err_cc_lInf( gridInd, ALGO, rep ),...
%                             err_cc_lInfElem( gridInd, ALGO, rep ),...
%                             err_cc_l1( gridInd, ALGO, rep )] = deal(NaN);
% %                             err_psd_l2( gridInd, ALGO, rep ),...
% %                             err_psd_l1( gridInd, ALGO, rep ),...
% %                             err_psd_lInfElem( gridInd, ALGO, rep ),...
%                             
%                         fprintf('Method is %18s, Rep %u failed!!! \n', NAME, rep );
%                         continue;
%                     end
                case 7
                    % Naive uniform, Leo's version
                    [XX,RowsMarker_uni] = baseline(X,compressRatio, 'all');
                    cc = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_uni,...
                        'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
                        'iscanonical',true,'m',m);
                    for skrepInd = 2:skrep
                        [XX,RowsMarker_uni] = baseline(X,compressRatio, 'all');
                        cc_tmp = findAutocorr_TimeBaseLine(XX,'RowsMarker',RowsMarker_uni,...
                            'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
                            'iscanonical',true,'m',m);
                        cc = cc+cc_tmp;
                    end
                    cc = cc/skrep;
            end
            mem    = whos('XX').bytes;
            if ALGO==6, mem = mem + whos('RowsMarker_sub1').bytes; end
            fprintf('Method is %18s, X is compressed by %g\n', NAME, mem0/mem );
            if issparse(X), fprintf('\tthough nnz is reduced by %g\n', numel(X)/nnz(XX) ); end
            
            storage( gridInd, ALGO, rep) = mem/mem0;
            
%             psd = psd_from_autocorr( cc );
%             FocusLags = 1:30;
            FocusLags = 1:MaxLag;
            err_cc_l2(  gridInd, ALGO, rep )   = norm(cc(FocusLags) - cc_orig(FocusLags))/norm(cc_orig(FocusLags));
            err_cc_lInf( gridInd, ALGO, rep )  = norm(cc(FocusLags) - cc_orig(FocusLags), Inf);
            err_cc_lInfElem( gridInd, ALGO, rep ) = errInfElem(cc(FocusLags),cc_orig(FocusLags));
            err_cc_l1(  gridInd, ALGO, rep )   = norm(cc(FocusLags)-cc_orig(FocusLags), 1)/norm(cc_orig(FocusLags),1);
            
%             err_psd_l2( gridInd, ALGO, rep )   = norm(psd - psd_orig)/norm(psd_orig);
%             err_psd_l1( gridInd, ALGO, rep ) = norm(psd - psd_orig, Inf);
%             err_psd_lInfElem( gridInd, ALGO, rep ) = errInfElem(psd,psd_orig);
            
        end
    end
    
end

%%
clear X XSlice XX
% save('Pulses_Res.mat');
save('BlockSimu_Jun.mat');

%% Plot
% 
% ALGONAMES = {'Haar','Gaussian','FJLT-Hadamard','UniformHD',...
%     'Particle Baseline', 'Time Baseline','Naive Uniform'};

figure('units','normalized','outerposition',[0 0 1 1]); clf;
errTensor_l2   = err_cc_l2; str_l2 = '\textbf{Relative $\ell_2$ error}';
% errTensor_lInf = err_cc_lInf; str_lInf = '\textbf{$\ell_{\infty}$ error}';
errTensor_lInfElem = err_cc_lInfElem; str_lInfElem = '\textbf{Relative $\ell_{\infty}$ error}';
% errTensor   = err_psd_zoom;   str = 'relative $\ell_2$ error in low frequencies';
% errTensor   = err_cc_l2; str = 'relative $\ell_2$ error for autocorrelation';
errTensor_l1   = err_cc_l1; str_l1 = '\textbf{Relative $\ell_1$ error}';

% compression  = mean(   storage, 3 );       % average compression
% err         = mean(   errTensor, 3 ); type = 'median';
err_l2       = mean( errTensor_l2, 3 ); % typel2 = 'median';
% err_lInf     = mean( errTensor_lInf, 3 ); % typelInf = 'median';
err_lInfElem = mean( errTensor_lInfElem, 3 ); % typelInf = 'median';
err_l1       = mean( errTensor_l1, 3 ); % type = 'min';


xlabstr = '\textbf{Total length of time series}';
% ALGOList= [1,2,3];
ALGOList= [1:6];

subplot(1,3,3)
for ALGO    = ALGOList
% for ALGO    = [1:3,5,7,8,9]
    % Right way:
    if ALGO == 1 || ALGO == 2 || ALGO == 3
        plot( Grid,err_lInfElem(:,ALGO), 'o--',...
            'linewidth',3.5,'DisplayName',ALGONAMES{ALGO},'markersize',16 );
        hold all
    else
        plot( Grid,err_lInfElem(:,ALGO), 'o--',...
            'linewidth',2,'DisplayName',ALGONAMES{ALGO},'markersize',16 );
        hold all
    end
end
holdMarker( 0 );
set(gca,'FontSize',24);
% ylim([0,1]);
% legend
xlabel(xlabstr,'interpreter','latex');
% ylabel(sprintf('Error (%s, %s over %d trials)',str,type,repeat),'interpreter','latex');
% ylabel(sprintf('%s, %s over %d trials',str_l2,type,repeat),'interpreter','latex');
ylabel(sprintf('%s',str_lInfElem),'interpreter','latex');
% axis tight;
set(gca,'xscale','log')
set(gca,'yscale','log')
% set(gca,'FontSize',22)


subplot(1,3,1)
% for ALGO    = [1:3,5,7,8,9]
for ALGO    = ALGOList
    % Wrong way (uses wrong compression ratio for x points)
%     plot( compressRatio_Grid,err(NInd,:,ALGO), 'o--', ...
%         'linewidth',2,'DisplayName',ALGONAMES{ALGO},'markersize',16  );
    
    % Right way:
    if ALGO == 1 || ALGO == 2 || ALGO == 3
        plot( Grid,err_l1(:,ALGO), 'o--',...
            'linewidth',3.5,'DisplayName',ALGONAMES{ALGO},'markersize',16 );
        hold all
    else
        plot( Grid,err_l1(:,ALGO), 'o--',...
            'linewidth',2,'DisplayName',ALGONAMES{ALGO},'markersize',16 );
        hold all
    end
end
holdMarker( 0 );
set(gca,'FontSize',24);
axis tight
% ylim([0,1]);
% xlim([1e-2 10^(-0.5)])
% axis tight 
xlabel(xlabstr,'interpreter','latex');
% ylabel(sprintf('Error (%s, %s over %d trials)',str,type,repeat),'interpreter','latex');
% ylabel(sprintf('%s, %s over %d trials',str_lInf,type,repeat),'interpreter','latex');
ylabel(sprintf('%s',str_l1),'interpreter','latex');
set(gca,'xscale','log')
set(gca,'yscale','log')


subplot(1,3,2)
for ALGO    = ALGOList
% for ALGO    = [1:3,5,7,8,9]
    % Right way:
    if ALGO == 1 || ALGO == 2 || ALGO == 3
        plot( Grid,err_l2(:,ALGO), 'o--',...
            'linewidth',3.5,'DisplayName',ALGONAMES{ALGO},'markersize',16 );
        hold all
    else
        plot( Grid,err_l2(:,ALGO), 'o--',...
            'linewidth',2,'DisplayName',ALGONAMES{ALGO},'markersize',16 );
        hold all
    end

end
holdMarker( 0 );
set(gca,'FontSize',24);
ylim([10^(-1.35),1]);
axis tight
xlabel(xlabstr,'interpreter','latex');
% ylabel(sprintf('Error (%s, %s over %d trials)',str,type,repeat),'interpreter','latex');
% ylabel(sprintf('%s, %s over %d trials',str_W,type,repeat),'interpreter','latex');
ylabel(sprintf('%s',str_l2),'interpreter','latex');
set(gca,'xscale','log')
set(gca,'yscale','log')
% set(gca,'FontName','Times')
legend
legend('location','best','interpreter','latex','NumColumns',7);


%%
export_fig('~/Desktop/research/sketching/figures/TvsError.pdf', '-pdf','-transparent');
