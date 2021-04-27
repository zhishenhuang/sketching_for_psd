%{ Autocorrelation Evaluation
%
%}

%% Load data
clear all; close all; clc;
% mode = 'synthetic';
mode = 'real';
rng(0)
addpath ~/Desktop/research/sketching/code_becker/
addpath ~/Desktop/research/sketching/MDSpectralAnalysis-master/Data/Methanol_Velocity_Data/
rng(2343);
[t,Vx,Vy,Vz] = loadMetOHVel();

[T,N] = size(Vx');
trueX     = Vx(:,ceil(T/2):T)'; % trueX is TIME BY SPACE
[T,~] = size(trueX);
%%
% Data Preparation: Remove Mean; Ground Truth Autocorrelation visualization
% removeMean  = true;
removeMean  = false;
if strcmpi(mode, 'real')
    X    = trueX;
end
if removeMean
    X    = X - mean(X);
end

posOnly     = true; 
% SCALEOPT    = 'coeff';
SCALEOPT    = 'unbiased';
% First, data for the true original signal
[cc_orig,lags] = findAutocorr(X,'posOnly',posOnly,...
    'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N );
% pwd_orig = periodogram(mean(X,2));
% pwd_orig = periodogram( cc_orig ); pwd_orig = sqrt(pwd_orig);
pwd_orig = abs( fft(cc_orig) ); 

if strcmpi(mode,'synthetic')
    pwd_orig = pwd_orig( round(1:end/2) + 1 );
else
    pwd_orig = pwd_orig( 1:end/2 + 1 );
end

% [pks_orig,loc_orig] = findpeaks(smoothdata(pwd_orig,'loess'),...
%     'MinPeakProminence',5e-6);
% [loc_orig,order] = sort(loc_orig);
% pks_orig = pks_orig(order);

figure(1); clf;
subplot(2,1,1)
plot(lags,cc_orig,'LineWidth',.5);
title('Autocorrelation ground truth')
set(gca,'FontSize',23)
axis tight
subplot(2,1,2)
freqBdd = 0.2;
lowfreq_cutoff = floor(length(pwd_orig)*freqBdd);
plot( pwd_orig  );
hold on;
plot(pwd_orig(1:lowfreq_cutoff),'color','r','LineWidth',2)
title('Power Density ground truth')
axis tight
set(gca,'FontSize',23)

%% Plot Test
figure(2);  %clf;
% periodogram( cc_orig )
% h = gca;
% h.Children.YData = .5*h.Children.YData;
hold all
semilogy( pwd_orig  )
%
figure(3); clf; 
periodogram( mean(X,2) )
figure(4); clf;
semilogy( smoothdata(pwd_orig));

%%
figure(1); clf;
findpeaks(smoothdata(pwd_orig,'loess'),'MinPeakProminence',5e-6);
hold on;
plot(pwd_orig,'r--','LineWidth',0.5);
ylim([0 7e-4])
%%

% compressRatio_Grid = [0.01];
compressRatio_Grid = 10.^(linspace(-2,-1,10));
% Tol    = 5e-2;
repeat = 5;
% [error_canonical,error_fjlt,error_haar,error_gaussian,error_uni,error_sub1,error_sub2] = ...
%     deal(zeros(length(compressRatio_Grid),repeat));
% [pwderror_canonical,pwderror_fjlt,pwderror_haar,pwderror_gaussian,pwderror_uni,...
%     pwderror_sub1,pwderror_sub2] = deal(zeros(length(compressRatio_Grid),repeat));
% [w2_canonical,w2_fjlt,w2_haar,w2_gaussian,w2_uni,...
%     w2_sub1,w2_sub2] = deal(zeros(length(compressRatio_Grid),repeat));

[da_canonical,da_fjlt,da_haar,da_gaussian,...
    da_uni,da_sub1,da_sub2] = deal(zeros(length(compressRatio_Grid),1));


% [pwdlowfreqerror_canonical,pwdlowfreqerror_fjlt,pwdlowfreqerror_haar,pwdlowfreqerror_gaussian,...
%     pwdlowfreqerror_uni,pwdlowfreqerror_sub1,pwdlowfreqerror_sub2] = deal(zeros(length(compressRatio_Grid),1));
%% plot compression ratio v.s. accuracy
for ind = 1:length(compressRatio_Grid)
    fprintf('   Current Compress Ratio Grid Point %u out of %u\n',ind,length(compressRatio_Grid));
    compressRatio = compressRatio_Grid(ind);
    m             = round(compressRatio*N);
    count = 1;
    while count<=repeat
        fprintf('   ~~~Repeat %u at Compress Ratio Grid Point~~~\n',count);
        
        HAAR         = sketch( m, N, 'haar' );
        X_HAAR       = HAAR(X')';
        cc_haar      = findAutocorr(X_HAAR,'posOnly',posOnly,...
            'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
        
        GAUSSIAN     = sketch( m, N, 'gaussian' );
        X_GAUSSIAN  = GAUSSIAN(X')';
        cc_gaussian  = findAutocorr(X_GAUSSIAN,'posOnly',posOnly,...
            'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
        
        FJLT         = sketch( m, N, 'fjlt' );
        X_FJLT       = FJLT(X')';
        cc_fjlt      = findAutocorr(X_FJLT,'posOnly',posOnly,...
            'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
        
        PRECONDHD    = sketch( m, N, 'canonical' );
        X_PRECONDHD  = PRECONDHD(X')';
        cc_canonical = findAutocorr(X_PRECONDHD,'posOnly',posOnly,...
            'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
            'iscanonical',true,'m',m);

        X_SUB2       = baseline(X, m, '2');
        cc_sub2      = findAutocorr(X_SUB2,'posOnly',posOnly,...
            'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/m);
        
        [X_UNI,RowsMarker_uni] = baseline(X,compressRatio, 'all');
        cc_uni       = findAutocorr_TimeBaseLine(X_UNI,'RowsMarker',RowsMarker_uni,...
            'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
            'iscanonical',true,'m',m);
%         X_UNI        = baseline(X',m, 'all-old')';
%         cc_uni_old       = findAutocorr_TimeBaseLine(X_UNI,'posOnly',posOnly,...
%             'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N,...
%             'iscanonical',true,'m',m);

        [X_SUB1,RowsMarker_sub1]  = baseline(X, round(compressRatio*T),'1-alternative');
%         cc_sub1_old     = findAutocorr_TimeBaseLine(X_SUB1,...
%             'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
        cc_sub1     = findAutocorr_TimeBaseLine(X_SUB1,'RowsMarker',RowsMarker_sub1,...
            'posOnly',posOnly,'removeMean',removeMean,'SCALEOPT',SCALEOPT,'scale',1/N);
        
        if strcmpi(mode,'synthetic')
            pwd_canonical= abs( fft(cc_canonical) ); pwd_canonical = pwd_canonical(  round(1:end/2) + 1 );
            pwd_gaussian = abs( fft(cc_gaussian) );  pwd_gaussian  = pwd_gaussian( round(1:end/2) + 1 );
            pwd_fjlt     = abs( fft(cc_fjlt) );      pwd_fjlt      = pwd_fjlt( round(1:end/2) + 1 );
            pwd_haar     = abs( fft(cc_haar) );      pwd_haar      = pwd_haar( round(1:end/2) + 1 );
            pwd_sub1     = abs( fft(cc_sub1) );      pwd_sub1      = pwd_sub1( round(1:end/2) + 1 );
            pwd_sub2     = abs( fft(cc_sub2) );      pwd_sub2      = pwd_sub2( round(1:end/2) + 1 );
            pwd_uni      = abs( fft(cc_uni) );       pwd_uni       = pwd_uni( round(1:end/2) + 1 );
        else
            pwd_canonical= abs( fft(cc_canonical) ); pwd_canonical = pwd_canonical( 1:end/2 + 1 );
            pwd_gaussian = abs( fft(cc_gaussian) );  pwd_gaussian  = pwd_gaussian( 1:end/2 + 1 );
            pwd_fjlt     = abs( fft(cc_fjlt) );      pwd_fjlt      = pwd_fjlt( 1:end/2 + 1 );
            pwd_haar     = abs( fft(cc_haar) );      pwd_haar      = pwd_haar( 1:end/2 + 1 );
            pwd_sub1     = abs( fft(cc_sub1) );      pwd_sub1      = pwd_sub1( 1:end/2 + 1 );
            pwd_sub2     = abs( fft(cc_sub2) );      pwd_sub2      = pwd_sub2( 1:end/2 + 1 );
            pwd_uni      = abs( fft(cc_uni) );       pwd_uni       = pwd_uni( 1:end/2 + 1 );
        end
        
        da_canonical(ind,count)  = dalpha(pwd_canonical, pwd_orig);
        da_gaussian(ind,count)   = dalpha(pwd_gaussian, pwd_orig);
        da_fjlt(ind,count)       = dalpha(pwd_fjlt, pwd_orig);
        da_haar(ind,count)       = dalpha(pwd_haar, pwd_orig);
        da_uni(ind,count)        = dalpha(pwd_uni, pwd_orig);
        da_sub1(ind,count)       = dalpha(pwd_sub1, pwd_orig(1:length(pwd_sub1)));
        da_sub2(ind,count)       = dalpha(pwd_sub2, pwd_orig);
        
%         w2_canonical(ind,count) = WassersteinDistance(pwd_canonical,pwd_orig);
%         w2_gaussian(ind,count) = WassersteinDistance(pwd_gaussian,pwd_orig);
%         w2_fjlt(ind,count) = WassersteinDistance(pwd_fjlt,pwd_orig);
%         w2_haar(ind,count) = WassersteinDistance(pwd_haar,pwd_orig);
%         w2_sub1(ind,count) = WassersteinDistance(pwd_sub1,pwd_orig);
%         w2_sub2(ind,count) = WassersteinDistance(pwd_sub2,pwd_orig);
%         w2_uni(ind,count) = WassersteinDistance(pwd_uni,pwd_orig);
        
%         [pks_canonical, loc_canonical] = FindPeaks(pwd_canonical);
%         [pks_gaussian,  loc_gaussian]  = FindPeaks(pwd_gaussian);
%         [pks_fjlt,      loc_fjlt]      = FindPeaks(pwd_fjlt);
%         [pks_haar,      loc_haar]      = FindPeaks(pwd_haar);
%         [pks_sub1,      loc_sub1]      = FindPeaks(pwd_sub1);
%         [pks_sub2,      loc_sub2]      = FindPeaks(pwd_sub2);
%         [pks_uni,       loc_uni]       = FindPeaks(pwd_uni);
        
%         pkCount_canonical(ind,count)   = length(ComparePeaks(pks_canonical,pks_orig,loc_canonical,loc_orig,length(pwd_orig),'Tol',Tol));
%         pkCount_gaussian(ind,count)    = length(ComparePeaks(pks_gaussian,pks_orig,loc_gaussian,loc_orig,length(pwd_orig),'Tol',Tol));
%         pkCount_fjlt(ind,count)        = length(ComparePeaks(pks_fjlt,pks_orig,loc_fjlt,loc_orig,length(pwd_orig),'Tol',Tol));
%         pkCount_haar(ind,count)        = length(ComparePeaks(pks_haar,pks_orig,loc_haar,loc_orig,length(pwd_orig),'Tol',Tol));
%         pkCount_sub1(ind,count)        = length(ComparePeaks(pks_sub1,pks_orig,loc_sub1,loc_orig,length(pwd_orig),'Tol',Tol));
%         pkCount_sub2(ind,count)        = length(ComparePeaks(pks_sub2,pks_orig,loc_sub2,loc_orig,length(pwd_orig),'Tol',Tol));
%         pkCount_uni(ind,count)         = length(ComparePeaks(pks_uni,pks_orig,loc_uni,loc_orig,length(pwd_orig),'Tol',Tol));
        
%         error_canonical(ind,count)     = norm(cc_canonical - cc_orig)/norm(cc_orig);
%         error_gaussian(ind,count)      = norm(cc_gaussian - cc_orig)/norm(cc_orig);
%         error_fjlt(ind,count)          = norm(cc_fjlt - cc_orig)/norm(cc_orig);
%         error_haar(ind,count)          = norm(cc_haar - cc_orig)/norm(cc_orig);
%         error_uni(ind,count)           = norm(cc_uni - cc_orig)/norm(cc_orig);
%         error_sub1(ind,count)          = norm(cc_sub1 - cc_orig)/norm(cc_orig);
%         error_sub2(ind,count)          = norm(cc_sub2 - cc_orig)/norm(cc_orig);
        
%         pwderror_canonical(ind,count)  = norm(pwd_canonical - pwd_orig)/norm(pwd_orig);
%         pwderror_gaussian(ind,count)   = norm(pwd_gaussian - pwd_orig)/norm(pwd_orig);
%         pwderror_fjlt(ind,count)       = norm(pwd_fjlt - pwd_orig)/norm(pwd_orig);
%         pwderror_haar(ind,count)       = norm(pwd_haar - pwd_orig)/norm(pwd_orig);
%         pwderror_uni(ind,count)        = norm(pwd_uni - pwd_orig)/norm(pwd_orig);
%         pwderror_sub1(ind,count)       = norm(pwd_sub1 - pwd_orig(1:length(pwd_sub1)))/norm(pwd_orig(1:length(pwd_sub1)));
%         pwderror_sub2(ind,count)       = norm(pwd_sub2 - pwd_orig)/norm(pwd_orig);
          
%         freqBdd = 0.15;
%         lowfreq_cutoff = floor(length(pwd_orig)*freqBdd);
%         pwdlowfreqerror_canonical(ind)  = pwdlowfreqerror_canonical(ind) + norm(pwd_canonical(1:lowfreq_cutoff) - pwd_orig(1:lowfreq_cutoff))/norm(pwd_orig(1:lowfreq_cutoff));
%         pwdlowfreqerror_gaussian(ind)   = pwdlowfreqerror_gaussian(ind)  + norm(pwd_gaussian(1:lowfreq_cutoff) - pwd_orig(1:lowfreq_cutoff))/norm(pwd_orig(1:lowfreq_cutoff));
%         pwdlowfreqerror_fjlt(ind)       = pwdlowfreqerror_fjlt(ind)      + norm(pwd_fjlt(1:lowfreq_cutoff) - pwd_orig(1:lowfreq_cutoff))/norm(pwd_orig(1:lowfreq_cutoff));
%         pwdlowfreqerror_haar(ind)       = pwdlowfreqerror_haar(ind)      + norm(pwd_haar(1:lowfreq_cutoff) - pwd_orig(1:lowfreq_cutoff))/norm(pwd_orig(1:lowfreq_cutoff));
%         pwdlowfreqerror_uni(ind)        = pwdlowfreqerror_uni(ind)       + norm(pwd_uni(1:lowfreq_cutoff) - pwd_orig(1:lowfreq_cutoff))/norm(pwd_orig(1:lowfreq_cutoff));
%         pwdlowfreqerror_sub1(ind)       = pwdlowfreqerror_sub1(ind)      + norm(pwd_sub1(1:lowfreq_cutoff) - pwd_orig(1:lowfreq_cutoff))/norm(pwd_orig(1:lowfreq_cutoff));
%         pwdlowfreqerror_sub2(ind)       = pwdlowfreqerror_sub2(ind)      + norm(pwd_sub2(1:lowfreq_cutoff) - pwd_orig(1:lowfreq_cutoff))/norm(pwd_orig(1:lowfreq_cutoff));

        count = count + 1;
    end
end

%%
figure(1); clf;
plot(compressRatio_Grid,mean(da_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(da_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(da_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(da_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(da_uni,2),'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid,mean(da_sub1,2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid,mean(da_sub2,2),'--','DisplayName','particle compression','LineWidth',2);
axis tight
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('$\ell_{\alpha\inf}$\textbf{ Error}','interpreter','latex')
legend('location','best','interpreter','latex','NumColumns',2);
% title('low freq')
set(gca,'FontSize',32)
set(gca,'FontName','Times')
% set(gca,'yscale','log')

%%
figure(1); clf;
subplot(1,2,1)
plot(compressRatio_Grid,pwdlowfreqerror_canonical,'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,pwdlowfreqerror_gaussian,'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,pwdlowfreqerror_fjlt,'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,pwdlowfreqerror_haar,'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,pwdlowfreqerror_uni,'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid,pwdlowfreqerror_sub1,'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid,pwdlowfreqerror_sub2,'--','DisplayName','particle compression','LineWidth',2);

xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('\textbf{Relative Error}','interpreter','latex')
legend('location','best','interpreter','latex');
title('low freq')
set(gca,'FontSize',32)
set(gca,'FontName','Times')
set(gca,'yscale','log')

subplot(1,2,2)
plot(compressRatio_Grid,pwdhighfreqerror_canonical,'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,pwdhighfreqerror_gaussian,'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,pwdhighfreqerror_fjlt,'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,pwdhighfreqerror_haar,'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,pwdhighfreqerror_uni,'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid,pwdhighfreqerror_sub1,'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid,pwdhighfreqerror_sub2,'--','DisplayName','particle compression','LineWidth',2);

xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('\textbf{Relative Error}','interpreter','latex')
legend('location','best','interpreter','latex');
title('high freq')
set(gca,'FontSize',32)
set(gca,'FontName','Times')
set(gca,'yscale','log')

%%
addpath ~/Downloads/altmany-export_fig-5b3965b/


%% Low sampling rate: autocorr accuracy
SmoothMethod = 'rlowess';
% figure(1); clf;
figure('units','normalized','outerposition',[0 0 1 1]); clf;
% plot(compressRatio_Grid,mean(error_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
plot(compressRatio_Grid,smooth(mean(error_canonical,2),SmoothMethod),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,smooth(mean(error_gaussian,2),SmoothMethod),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,smooth(mean(error_fjlt,2),SmoothMethod),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,smooth(mean(error_haar,2),SmoothMethod),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(error_uni,2),'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid,mean(error_sub1,2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid,mean(error_sub2,2),'--','DisplayName','particle compression','LineWidth',2);
axis tight
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('\textbf{Relative Error}','interpreter','latex')
legend('location','best','interpreter','latex');
set(gca,'FontSize',32)
ylim([0 3])
yticks(0:.25:3)
set(gca,'FontName','Times')

% export_fig('autocorr_accuracy_low_samp.pdf', '-pdf','-transparent');

%% high sampling rate
figure('units','normalized','outerposition',[0 0 1 1]); clf;
plot(compressRatio_Grid(1:end-1),mean(error_canonical(1:end-1,:),2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid(1:end-1),mean(error_gaussian(1:end-1,:),2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid(1:end-1),mean(error_fjlt(1:end-1,:),2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid(1:end-1),mean(error_haar(1:end-1,:),2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid(1:end-1),mean(error_uni(1:end-1,:),2),'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid(1:end-1),mean(error_sub1(1:end-1,:),2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid(1:end-1),mean(error_sub2(1:end-1,:),2),'--','DisplayName','particle compression','LineWidth',2);
axis tight
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('\textbf{Relative Error}','interpreter','latex')
legend('location','best','interpreter','latex');
set(gca,'FontSize',32)
ylim([10^(-2.5) 1e0])
set(gca,'FontName','Times')
set(gca,'yscale','log')

export_fig('autocorr_accuracy_high_samp.pdf', '-pdf','-transparent');

%% pwd error low samp
figure(2);clf;
plot(compressRatio_Grid,mean(pwderror_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(pwderror_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(pwderror_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(pwderror_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(pwderror_uni,2),'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid,mean(pwderror_sub1,2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid,mean(pwderror_sub2,2),'--','DisplayName','particle compression','LineWidth',2);
axis tight
% ylim([0 2.5])
% yticks(0:.25:2.5)
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('\textbf{Relative Error}','interpreter','latex')
legend('location','best','interpreter','latex');
set(gca,'FontSize',32)
set(gca,'FontName','Times')



%% pwd error lowfreq
figure('units','normalized','outerposition',[0 0 1 1]); clf;

subplot(1,2,1)
lowfreq_cutoff = floor(length(pwd_orig)/10);
f = (1:1:round(5002/10))/50;
% plot(f,pwd_haar,'LineWidth',1,'DisplayName','Haar');
% hold all;
% plot(f,pwd_gaussian,'LineWidth',1,'DisplayName','Gaussian');
% plot(f,pwd_canonical,'LineWidth',1,'DisplayName','Uniform-HD');
% plot(f,pwd_fjlt,'LineWidth',1,'DisplayName','fjlt');
plot(f,smoothdata(pwd_orig(1:length(f)),'loess'),'r','LineWidth',2,'DisplayName','Ground Truth');
% plot(f,pwd_sub1,'DisplayName','Time Compression','LineWidth',1);
% plot(f,pwd_sub2,'DisplayName','Particle Compression','LineWidth',1);
% plot(f,pwd_uni,'DisplayName','Naive Uniform','LineWidth',1);
% axis tight
legend('location','best','interpreter','latex');
xlabel('\textbf{Frequency} (THz)','interpreter','latex')
ylabel('\textbf{Power Spectral Density} (a.u.)','interpreter','latex')
% ylim([-2 10])
% yticks([-2:2:10])
set(gca,'FontSize',28)
set(gca,'FontName','Times')
set(gca,'YTick',[])


subplot(1,2,2)
plot(compressRatio_Grid,mean(pwdlowfreqerror_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(pwdlowfreqerror_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(pwdlowfreqerror_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(pwdlowfreqerror_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(pwdlowfreqerror_uni,2),'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid,mean(pwdlowfreqerror_sub1,2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid,mean(pwdlowfreqerror_sub2,2),'--','DisplayName','particle compression','LineWidth',2);
axis tight
ylim([.1 1.7])
yticks([0.1,0.2:.2:1.7])
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('\textbf{Relative Error}','interpreter','latex')
legend('location','best','interpreter','latex');
% title('low freq')
set(gca,'FontSize',32)
set(gca,'FontName','Times')

% export_fig('pwd_lowfreq.pdf', '-pdf','-transparent');



%% pwd error high samp
figure(2);clf;
plot(compressRatio_Grid(1:9),mean(pwderror_canonical(1:9,:),2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid(1:9),mean(pwderror_gaussian(1:9,:),2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pwderror_fjlt(1:9,:),2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pwderror_haar(1:9,:),2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pwderror_uni,2),'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pwderror_sub1,2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pwderror_sub2(1:9,:),2),'--','DisplayName','particle compression','LineWidth',2);
axis tight
ylim([10^(-2.5) 1])
% yticks(0:.25:2.5)
xlabel('\textbf{Compress Ratio}','interpreter','latex')
ylabel('\textbf{Relative Error}','interpreter','latex')
legend('location','best','interpreter','latex');
set(gca,'FontSize',32)
set(gca,'FontName','Times')
set(gca,'yscale','log')

%% Wasserstein Distance Plot
figure('units','normalized','outerposition',[0 0 1 1]); clf;

subplot(2,2,1);
plot(compressRatio_Grid,mean(w2_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(w2_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(w2_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(w2_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(w2_uni,2),'--','DisplayName','naive uniform','LineWidth',3);
plot(compressRatio_Grid,mean(w2_sub1,2),'--','DisplayName','time compression','LineWidth',3);
plot(compressRatio_Grid,mean(w2_sub2,2),'--','DisplayName','particle compression','LineWidth',3);
% axis tight
xlim([0.01 0.1])
% xlim([0.1 0.9])
% xticks(0.1:0.1:0.9)
yticks(10.^[1,2,3])
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('$W_2$','interpreter','latex')
% legend('location','best','interpreter','latex');
set(gca,'FontSize',25)
set(gca,'FontName','Times')
set(gca,'yscale','log')


subplot(2,2,2); 
plot(compressRatio_Grid,mean(w1_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(w1_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(w1_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(w1_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(w1_uni,2),'--','DisplayName','naive uniform','LineWidth',3);
plot(compressRatio_Grid,mean(w1_sub1,2),'--','DisplayName','time compression','LineWidth',3);
plot(compressRatio_Grid,mean(w1_sub2,2),'--','DisplayName','particle compression','LineWidth',3);
% axis tight
xlim([0.01 0.1])
% xlim([0.1 0.9])
% xticks(0.1:0.1:0.9)
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('$W_1$','interpreter','latex')
% legend('location','best','interpreter','latex');
set(gca,'FontSize',25)
set(gca,'FontName','Times')
set(gca,'yscale','log')


subplot(2,2,3);
plot(compressRatio_Grid,mean(c2_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(c2_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(c2_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(c2_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(c2_uni,2),'--','DisplayName','naive uniform','LineWidth',3);
plot(compressRatio_Grid,mean(c2_sub1,2),'--','DisplayName','time compression','LineWidth',3);
plot(compressRatio_Grid,mean(c2_sub2,2),'--','DisplayName','particle compression','LineWidth',3);
% axis tight
xlim([0.01 0.1])
xticks(10.^([-1,0,1,2]))
% xlim([0.1 0.9])
% xticks(0.1:0.1:0.9)
ylim([1e-1 1e1])
yticks(10.^[-1,0,1])
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('$\textbf{CDF-Dist}_2$','interpreter','latex')
% legend('location','best','interpreter','latex');
set(gca,'FontSize',25)
set(gca,'FontName','Times')
set(gca,'yscale','log')


subplot(2,2,4);
plot(compressRatio_Grid,mean(linf_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(linf_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(linf_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(linf_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(linf_uni,2),'--','DisplayName','naive uniform','LineWidth',3);
plot(compressRatio_Grid,mean(linf_sub1,2),'--','DisplayName','time compression','LineWidth',3);
plot(compressRatio_Grid,mean(linf_sub2,2),'--','DisplayName','particle compression','LineWidth',3);
% axis tight
xlim([0.01 0.1])
% xlim([0.1 0.9])
% xticks(0.1:0.1:0.9)
xlabel('\textbf{Compression Ratio}','interpreter','latex')
ylabel('$\ell_\infty$ \textbf{distance}','interpreter','latex')
legend('location','best','interpreter','latex','NumColumns',7);
set(gca,'FontSize',25)
set(gca,'FontName','Times')
set(gca,'yscale','log')



%% 
export_fig('four_metrics_low.pdf', '-pdf','-transparent');

%%
figure(3); clf;
plot(compressRatio_Grid,mean(pkCount_canonical,2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid,mean(pkCount_gaussian,2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid,mean(pkCount_fjlt,2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid,mean(pkCount_haar,2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid,mean(pkCount_uni,2),'--','DisplayName','naive uniform','LineWidth',5);
plot(compressRatio_Grid,mean(pkCount_sub1,2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid,mean(pkCount_sub2,2),'--','DisplayName','particle compression','LineWidth',2);
plot(compressRatio_Grid,length(pks_orig)*ones(length(compressRatio_Grid),1),':','DisplayName','ground truth','LineWidth',2);
axis tight
ylim([0 13.5])
yticks([0:2:12,13])
xlabel('\textbf{Compress Ratio}','interpreter','latex')
ylabel('\textbf{Peak Counts}','interpreter','latex')
legend('location','best','interpreter','latex');
set(gca,'FontSize',32)
set(gca,'FontName','Times')

%% pkcount high samp rate
figure(3); clf;
plot(compressRatio_Grid(1:9),mean(pkCount_canonical(1:9,:),2),'DisplayName','Uniform-HD','LineWidth',2);
hold all;
plot(compressRatio_Grid(1:9),mean(pkCount_gaussian(1:9,:),2),'DisplayName','Gaussian','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pkCount_fjlt(1:9,:),2),'DisplayName','FJLT','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pkCount_haar(1:9,:),2),'DisplayName','Haar','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pkCount_uni,2),'--','DisplayName','naive uniform','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pkCount_sub1,2),'--','DisplayName','time compression','LineWidth',2);
plot(compressRatio_Grid(1:9),mean(pkCount_sub2(1:9,:),2),'--','DisplayName','particle compression','LineWidth',2);
plot(compressRatio_Grid(1:9),13*ones(9,1),':','DisplayName','ground truth','LineWidth',2);
axis tight
ylim([0 13.5])
yticks([0:2:12,13])
xlabel('\textbf{Compress Ratio}','interpreter','latex')
ylabel('\textbf{Peak Counts}','interpreter','latex')
legend('location','best','interpreter','latex');
set(gca,'FontSize',32)
set(gca,'FontName','Times')


%% autocorr plot , low samp rate regime
figure('units','normalized','outerposition',[0 0 1 1]); clf;
subplot(2,1,1)
plot(cc_haar(1:1e2),'LineWidth',1.5,'DisplayName','Haar');
hold all;
plot(cc_gaussian(1:1e2),'LineWidth',1.5,'DisplayName','Gaussian');
plot(cc_canonical(1:1e2),'LineWidth',1.5,'DisplayName','Uniform-HD');
plot(cc_fjlt(1:1e2),'LineWidth',1.5,'DisplayName','FJLT');
plot(cc_sub2(1:1e2),'LineWidth',1.5,'DisplayName','Particle compression')
plot(cc_sub1(1:1e2),'LineWidth',1.5,'DisplayName','Time compression')
plot(cc_uni(1:1e2),'LineWidth',1.5,'DisplayName','Naive uniform')
plot(cc_orig(1:1e2),'--k','LineWidth',1.5,'DisplayName','Ground Truth');
xlabel('\textbf{Lags}','interpreter','latex') % need to adjust units
ylabel('\textbf{Autocorrelation } $R_\tau$','interpreter','latex')
% title('Autocorrelation --- compress ratio 10%')
xticks(0:10:1e2)
axis tight
legend('location','best','interpreter','latex','NumColumns',8);
set(gca,'FontSize',25)
set(gca,'FontName','Times')

subplot(2,1,2)
plot(abs(cc_haar(1:1e2)-cc_orig(1:1e2)),'LineWidth',1.5,'DisplayName','Haar');
hold all;
plot(abs(cc_gaussian(1:1e2)-cc_orig(1:1e2)),'LineWidth',1.5,'DisplayName','Gaussian');
plot(abs(cc_canonical(1:1e2)-cc_orig(1:1e2)),'LineWidth',1.5,'DisplayName','Uniform-HD');
plot(abs(cc_fjlt(1:1e2)-cc_orig(1:1e2)),'LineWidth',1.5,'DisplayName','FJLT');
plot(abs(cc_sub2(1:1e2)-cc_orig(1:1e2)),'LineWidth',1.5,'DisplayName','Particle compression')
plot(abs(cc_sub1(1:1e2)-cc_orig(1:1e2)),'LineWidth',1.5,'DisplayName','Time compression')
plot(abs(cc_uni(1:1e2)-cc_orig(1:1e2)),'LineWidth',1.5,'DisplayName','Naive uniform')
xlabel('\textbf{Lags}','interpreter','latex') % need to adjust units
ylabel('\textbf{Absolute Error of } $R_\tau$','interpreter','latex')
% title('Autocorrelation --- compress ratio 10%')
% xticks(0:10:1e2)
% axis tight
% legend('location','best','interpreter','latex');
set(gca,'FontSize',25)
set(gca,'FontName','Times')
% set(gca,'yscale','log')

% subplot(3,1,3)
% plot((cc_haar(1:1e2)-cc_orig(1:1e2))./cc_orig(1:1e2),'LineWidth',2,'DisplayName','Haar');
% hold all;
% plot((cc_gaussian(1:1e2)-cc_orig(1:1e2))./cc_orig(1:1e2),'LineWidth',2,'DisplayName','Gaussian');
% plot((cc_canonical(1:1e2)-cc_orig(1:1e2))./cc_orig(1:1e2),'LineWidth',2,'DisplayName','Uniform-HD');
% plot((cc_fjlt(1:1e2)-cc_orig(1:1e2))./cc_orig(1:1e2),'LineWidth',2,'DisplayName','FJLT');
% plot((cc_sub2(1:1e2)-cc_orig(1:1e2))./cc_orig(1:1e2),'LineWidth',2,'DisplayName','Particle dim compression')
% plot((cc_sub1(1:1e2)-cc_orig(1:1e2))./cc_orig(1:1e2),'LineWidth',2,'DisplayName','Time dim compression')
% plot((cc_uni(1:1e2)-cc_orig(1:1e2))./cc_orig(1:1e2),'LineWidth',2,'DisplayName','Naive uniform')
% xlabel('\textbf{Lags}','interpreter','latex') % need to adjust units
% ylabel('\textbf{Relative Evaluation Error of } $R_\tau$','interpreter','latex')
% % title('Autocorrelation --- compress ratio 10%')
% % xticks(0:10:1e2)
% % axis tight
% legend('location','best','interpreter','latex');
% set(gca,'FontSize',12)
% set(gca,'FontName','Times')
% % set(gca,'yscale','log')

% export_fig('autocorr_eg.pdf', '-pdf','-transparent');


%% autocorr plot --- whole regime
figure(3);clf;
plot(cc_haar,'LineWidth',2,'DisplayName','Haar');
hold all;
plot(cc_gaussian,'LineWidth',2,'DisplayName','Gaussian');
plot(cc_canonical,'LineWidth',2,'DisplayName','Uniform-HD');
plot(cc_fjlt,'LineWidth',2,'DisplayName','fjlt');
plot(cc_orig,'--','LineWidth',1,'DisplayName','Ground Truth');
plot(cc_sub2,'LineWidth',2,'DisplayName','Particle dim compression')
plot(cc_sub1,'LineWidth',2,'DisplayName','Time dim compression')
plot(cc_uni,'LineWidth',2,'DisplayName','Naive uniform')
xlabel('Lags') % need to adjust units
title('Autocorrelation --- compress ratio 10%')
axis tight
legend;
set(gca,'FontSize',28)



%% pwr spectrum plot
figure(4);clf;
f = (1:1:5002)/50;
plot(f,pwd_haar,'LineWidth',1,'DisplayName','Haar');
hold all;
% plot(f,pwd_gaussian,'LineWidth',1,'DisplayName','Gaussian');
% plot(f,pwd_canonical,'LineWidth',1,'DisplayName','Uniform-HD');
% plot(f,pwd_fjlt,'LineWidth',1,'DisplayName','fjlt');
plot(f,pwd_orig,'--','LineWidth',1,'DisplayName','Ground Truth');
plot(f,pwd_sub1,'DisplayName','Time Compression','LineWidth',1);
% plot(f,pwd_sub2,'DisplayName','Particle Compression','LineWidth',1);
plot(f,pwd_uni,'DisplayName','Naive Uniform','LineWidth',1);
axis tight
legend('location','best','interpreter','latex');
xlabel('\textbf{Frequency} (THz)','interpreter','latex')
ylabel('\textbf{Power Spectral Density} (a.u.)','interpreter','latex')
set(gca,'FontSize',28)
set(gca,'FontName','Times')
set(gca,'YTick',[])

%% Autocorrelation computation with sketching and its accuracy versus total time length
% average over multiple repeats on fixed percentage
[T,N] = size(Vx');
TMax   = T; % TMax = 20001
percentage_grid = 10.^(linspace(-4,-1,30)); % T = 20 to 2000
subsample_ratio = .1;
repeat = 10;
[error_canonical_R1F, error_gaussian_R1F, error_haar_R1F, error_fjlt_R1F, error_baseline_R1F...
    ] = deal(zeros(length(percentage_grid),repeat));
[error_canonical_R12, error_gaussian_R12, error_haar_R12, error_fjlt_R12, error_baseline_R12...
    ] = deal(zeros(length(percentage_grid),repeat));
[error_canonical_R2F, error_gaussian_R2F, error_haar_R2F, error_fjlt_R2F, error_baseline_R2F...
    ] = deal(zeros(length(percentage_grid),repeat));
[error_canonical_R22, error_gaussian_R22, error_haar_R22, error_fjlt_R22, error_baseline_R22...
    ] = deal(zeros(length(percentage_grid),repeat));
[error_canonical_RinfF, error_gaussian_RinfF, error_haar_RinfF, error_fjlt_RinfF, error_baseline_RinfF...
    ] = deal(zeros(length(percentage_grid),repeat));
[error_canonical_Rinf2, error_gaussian_Rinf2, error_haar_Rinf2, error_fjlt_Rinf2, error_baseline_Rinf2...
    ] = deal(zeros(length(percentage_grid),repeat));

PRECOND_CANONICAL   = true;
SUBSAMPLE           = true;
NORMALIZATION       = false;
SCALE_N             = true;

errFun_1F   = @(cc,cc_orig,Sigma,SigmaHat) norm(cc-cc_orig(1:length(cc)),1)/norm(Sigma-SigmaHat,'fro');
errFun_12   = @(cc,cc_orig,Sigma,SigmaHat) norm(cc-cc_orig(1:length(cc)),1)/norm(Sigma-SigmaHat,2);
errFun_2F   = @(cc,cc_orig,Sigma,SigmaHat) norm(cc-cc_orig(1:length(cc)),2)/norm(Sigma-SigmaHat,'fro');
errFun_22   = @(cc,cc_orig,Sigma,SigmaHat) norm(cc-cc_orig(1:length(cc)),2)/norm(Sigma-SigmaHat,2);
errFun_infF = @(cc,cc_orig,Sigma,SigmaHat) norm(cc-cc_orig(1:length(cc)),inf)/norm(Sigma-SigmaHat,'fro');
errFun_inf2 = @(cc,cc_orig,Sigma,SigmaHat) norm(cc-cc_orig(1:length(cc)),inf)/norm(Sigma-SigmaHat,2);

Told = 0;
SCALE_ERR = 1;
for ind = 1:length(percentage_grid)
    fprintf('    Current compress grid point %u out of %u\n',ind,length(percentage_grid));
    T = ceil(TMax*percentage_grid(ind));
    if Told == 0
        data = X( 1:T, : );
        Sigma = data*data';
    else
        dataNew = X( Told+1:T, : );
        AB = data*dataNew';
        Sigma = [Sigma, AB; AB', dataNew*dataNew'];
        data = [data;dataNew];        
%         if ind==3
%             % do a check
%             dataNew = X(1:T,:);
%             SigmaCheck = dataNew*dataNew';
%             fprintf(2,'Checking: %g\n', norm( SigmaCheck - Sigma, 'fro' )/norm(SigmaCheck,'fro') );
%         end
    end
    Told = T;
    
    if SCALE_N
        scale = 1/N;
    else
        scale = 1;
    end
    
    cc_orig = autocorr(Sigma,'SCALEOPT', SCALEOPT,...
        'Normalization',NORMALIZATION,'scale',scale);

    repeat_ind = 1;
    while repeat_ind <= repeat
        fprintf('Time dim size: %u \t Subsample ratio: %.3f \t Repeat: %d \n', ...
            T, subsample_ratio, repeat_ind);
        % compute unbiased estimator for covariance matrix;
        [C_canonical,~,~] = covariance_compute(data, 'subsample_ratio',subsample_ratio,...
    'PRECOND',PRECOND_CANONICAL, 'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','canonical');
        fprintf('  Method %s done \n','canonical');
        
        [C_gaussian,~,~]  = covariance_compute(data,'subsample_ratio',subsample_ratio,...
    'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','gaussian');     
        fprintf('  Method %s done \n','gaussian');
        
        [C_fjlt,~,~]      = covariance_compute(data,'subsample_ratio',subsample_ratio,...
    'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','fjlt');      
        fprintf('  Method %s done \n','fjlt');
        
        [C_haar,~,~]      = covariance_compute(data,'subsample_ratio',subsample_ratio,...
    'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','haar');      
        fprintf('  Method %s done \n','haar');
        
        [C_baseline,~]    = covariance_compute(data,'subsample_ratio',subsample_ratio,...
    'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','baseline');
        fprintf('  Method %s done \n','baseline');
        
        cc_canonical = autocorr(C_canonical,'SCALEOPT', SCALEOPT,'Normalization',NORMALIZATION,'scale',scale);
        cc_gaussian  = autocorr(C_gaussian, 'SCALEOPT', SCALEOPT,'Normalization',NORMALIZATION,'scale',scale);
        cc_fjlt      = autocorr(C_fjlt,     'SCALEOPT', SCALEOPT,'Normalization',NORMALIZATION,'scale',scale);
        cc_haar      = autocorr(C_haar,     'SCALEOPT', SCALEOPT,'Normalization',NORMALIZATION,'scale',scale);
        cc_baseline  = autocorr(C_baseline, 'SCALEOPT', SCALEOPT,'Normalization',NORMALIZATION,'scale',scale);       
        
        error_canonical_R1F(ind,repeat_ind)= errFun_1F(cc_canonical,cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_canonical);
        error_gaussian_R1F(ind,repeat_ind) = errFun_1F(cc_gaussian, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_gaussian);
        error_fjlt_R1F(ind,repeat_ind)     = errFun_1F(cc_fjlt,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_fjlt);
        error_haar_R1F(ind,repeat_ind)     = errFun_1F(cc_haar,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_haar);
        error_baseline_R1F(ind,repeat_ind) = errFun_1F(cc_baseline, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_baseline);
        
        error_canonical_R12(ind,repeat_ind)= errFun_12(cc_canonical,cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_canonical);
        error_gaussian_R12(ind,repeat_ind) = errFun_12(cc_gaussian, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_gaussian);
        error_fjlt_R12(ind,repeat_ind)     = errFun_12(cc_fjlt,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_fjlt);
        error_haar_R12(ind,repeat_ind)     = errFun_12(cc_haar,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_haar);
        error_baseline_R12(ind,repeat_ind) = errFun_12(cc_baseline, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_baseline);
        
        error_canonical_R2F(ind,repeat_ind)= errFun_2F(cc_canonical,cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_canonical);
        error_gaussian_R2F(ind,repeat_ind) = errFun_2F(cc_gaussian, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_gaussian);
        error_fjlt_R2F(ind,repeat_ind)     = errFun_2F(cc_fjlt,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_fjlt);
        error_haar_R2F(ind,repeat_ind)     = errFun_2F(cc_haar,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_haar);
        error_baseline_R2F(ind,repeat_ind) = errFun_2F(cc_baseline, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_baseline);
        
        error_canonical_R22(ind,repeat_ind)= errFun_22(cc_canonical,cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_canonical);
        error_gaussian_R22(ind,repeat_ind) = errFun_22(cc_gaussian, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_gaussian);
        error_fjlt_R22(ind,repeat_ind)     = errFun_22(cc_fjlt,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_fjlt);
        error_haar_R22(ind,repeat_ind)     = errFun_22(cc_haar,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_haar);
        error_baseline_R22(ind,repeat_ind) = errFun_22(cc_baseline, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_baseline);

        error_canonical_RinfF(ind,repeat_ind)= errFun_infF(cc_canonical,cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_canonical);
        error_gaussian_RinfF(ind,repeat_ind) = errFun_infF(cc_gaussian, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_gaussian);
        error_fjlt_RinfF(ind,repeat_ind)     = errFun_infF(cc_fjlt,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_fjlt);
        error_haar_RinfF(ind,repeat_ind)     = errFun_infF(cc_haar,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_haar);
        error_baseline_RinfF(ind,repeat_ind) = errFun_infF(cc_baseline, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_baseline);
        
        error_canonical_Rinf2(ind,repeat_ind)= errFun_inf2(cc_canonical,cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_canonical);
        error_gaussian_Rinf2(ind,repeat_ind) = errFun_inf2(cc_gaussian, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_gaussian);
        error_fjlt_Rinf2(ind,repeat_ind)     = errFun_inf2(cc_fjlt,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_fjlt);
        error_haar_Rinf2(ind,repeat_ind)     = errFun_inf2(cc_haar,     cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_haar);
        error_baseline_Rinf2(ind,repeat_ind) = errFun_inf2(cc_baseline, cc_orig,SCALE_ERR*Sigma,SCALE_ERR*C_baseline);

        repeat_ind = repeat_ind + 1;
    end
end

error_canonical_1F = mean( error_canonical_R1F, 2 );
error_gaussian_1F  = mean( error_gaussian_R1F, 2 );
error_fjlt_1F      = mean( error_fjlt_R1F, 2 );
error_haar_1F      = mean( error_haar_R1F, 2 );
error_baseline_1F    = mean( error_baseline_R1F, 2 );

error_canonical_2F = mean( error_canonical_R2F, 2 );
error_gaussian_2F  = mean( error_gaussian_R2F, 2 );
error_fjlt_2F      = mean( error_fjlt_R2F, 2 );
error_haar_2F      = mean( error_haar_R2F, 2 );
error_baseline_2F  = mean( error_baseline_R2F, 2 );

error_canonical_12 = mean( error_canonical_R12, 2 );
error_gaussian_12  = mean( error_gaussian_R12, 2 );
error_fjlt_12      = mean( error_fjlt_R12, 2 );
error_haar_12      = mean( error_haar_R12, 2 );
error_baseline_12  = mean( error_baseline_R12, 2 );

error_canonical_22 = mean( error_canonical_R22, 2 );
error_gaussian_22  = mean( error_gaussian_R22, 2 );
error_fjlt_22      = mean( error_fjlt_R22, 2 );
error_haar_22      = mean( error_haar_R22, 2 );
error_baseline_22  = mean( error_baseline_R22, 2 );

error_canonical_infF = mean( error_canonical_RinfF, 2 );
error_gaussian_infF  = mean( error_gaussian_RinfF, 2 );
error_fjlt_infF      = mean( error_fjlt_R22, 2 );
error_haar_infF      = mean( error_haar_RinfF, 2 );
error_baseline_infF  = mean( error_baseline_RinfF, 2 );

error_canonical_inf2 = mean( error_canonical_Rinf2, 2 );
error_gaussian_inf2  = mean( error_gaussian_Rinf2, 2 );
error_fjlt_inf2      = mean( error_fjlt_Rinf2, 2 );
error_haar_inf2      = mean( error_haar_Rinf2, 2 );
error_baseline_inf2  = mean( error_baseline_Rinf2, 2 );

%% error plot: check the tightness of the bound

time_grid = TMax * percentage_grid;
ubd1F     = sqrt(1+log(time_grid))/N;
ubdinfF   = 1/N * ones(length(time_grid),1);

figure(1);clf;
subplot(2,1,1); 
plot(time_grid,error_canonical_1F,'linewidth',2,'DisplayName','Uniform-HD');
hold all;
plot(time_grid,error_fjlt_1F,'linewidth',2,'DisplayName','FJLT');
plot(time_grid,error_gaussian_1F,'linewidth',2,'DisplayName','Gaussian');
plot(time_grid,error_haar_1F,'linewidth',2,'DisplayName','Haar');
plot(time_grid,ubd1F,'r:','linewidth',3,'DisplayName','Upper Bound')
legend('Location','best','interpreter','latex');
xlabel('\textbf{Lags}','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_1}{\|\Sigma-\hat{\Sigma}\|_F}$','interpreter','latex')
title('$\mathbf{1\to F}$','interpreter','latex')
set(gca,'xscale','log')
set(gca,'FontSize',26,'FontWeight','bold');
set(get(gca,'YLabel'),'Rotation',90)
xlim([0 2200])
% axis tight

subplot(2,1,2); 
plot(time_grid,error_canonical_infF,'linewidth',2,'DisplayName','Uniform-HD');
hold all;
plot(time_grid,error_fjlt_infF,'linewidth',2,'DisplayName','FJLT');
plot(time_grid,error_gaussian_infF,'linewidth',2,'DisplayName','Gaussian');
plot(time_grid,error_haar_infF,'linewidth',2,'DisplayName','Haar');
plot(time_grid,ubdinfF,'r:','linewidth',3,'DisplayName','Upper Bound')
legend('Location','best','interpreter','latex');
xlabel('\textbf{Lags}','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_\infty}{\|\Sigma-\hat{\Sigma}\|_F}$','interpreter','latex')
title('$\mathbf{\infty\to F}$','Interpreter','latex')
set(gca,'xscale','log')
set(gca,'FontSize',28,'FontWeight','bold');
set(get(gca,'YLabel'),'Rotation',90)
xlim([0 2200])
% axis tight


%% error plot: check the tightness of the bound

time_grid = TMax * percentage_grid;
ubd1F     = sqrt(1+log(time_grid))/N;
ubdinfF   = 1/N * ones(length(time_grid),1);

figure(1);clf;
subplot(2,3,1); 
plot(time_grid,error_canonical_1F,'k-','linewidth',2);
hold all;
plot(time_grid,error_fjlt_1F,'b--','linewidth',2);
plot(time_grid,error_gaussian_1F,'m:','linewidth',2);
plot(time_grid,error_haar_1F,'c-.','linewidth',2);
plot(time_grid,ubd1F,'r-','linewidth',3)
legend('canonical','fjlt','gaussian','haar','TUB')
xlabel('Lags','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_1}{\|\Sigma-\hat{\Sigma}\|_F}$','interpreter','latex')
title('1F')
set(gca,'xscale','log')
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Rotation',90)
axis tight

subplot(2,3,2); 
plot(time_grid,error_canonical_2F,'k-','linewidth',2);
hold all;
plot(time_grid,error_fjlt_2F,'b--','linewidth',2);
plot(time_grid,error_gaussian_2F,'m:','linewidth',2);
plot(time_grid,error_haar_2F,'c-.','linewidth',2);
legend('canonical','fjlt','gaussian','haar')
xlabel('Lags','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_2}{\|\Sigma-\hat{\Sigma}\|_F}$','interpreter','latex')
title('2F')
set(gca,'xscale','log')
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Rotation',90)
axis tight

subplot(2,3,3); 
plot(time_grid,error_canonical_infF,'k-','linewidth',2);
hold all;
plot(time_grid,error_fjlt_infF,'b--','linewidth',2);
plot(time_grid,error_gaussian_infF,'m:','linewidth',2);
plot(time_grid,error_haar_infF,'c-.','linewidth',2);
plot(time_grid,ubdinfF,'r-','linewidth',3)
legend('canonical','fjlt','gaussian','haar','TUB')
xlabel('Lags','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_\infty}{\|\Sigma-\hat{\Sigma}\|_F}$','interpreter','latex')
title('infF')
set(gca,'xscale','log')
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Rotation',90)
axis tight

subplot(2,3,4); 
plot(time_grid,error_canonical_12,'k-','linewidth',2);
hold all;
plot(time_grid,error_fjlt_12,'b--','linewidth',2);
plot(time_grid,error_gaussian_12,'m:','linewidth',2);
plot(time_grid,error_haar_12,'c-.','linewidth',2);
legend('canonical','fjlt','gaussian','haar')
xlabel('Lags','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_1}{\|\Sigma-\hat{\Sigma}\|_2}$','interpreter','latex')
title('12')
set(gca,'xscale','log')
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Rotation',90)
axis tight

subplot(2,3,5);
plot(time_grid,error_canonical_22,'k-','linewidth',2);
hold all;
plot(time_grid,error_fjlt_22,'b--','linewidth',2);
plot(time_grid,error_gaussian_22,'m:','linewidth',2);
plot(time_grid,error_haar_22,'c-.','linewidth',2);
legend('canonical','fjlt','gaussian','haar')
xlabel('lags','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_2}{\|\Sigma-\hat{\Sigma}\|_2}$','interpreter','latex')
title('22')
set(gca,'xscale','log')
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Rotation',90)
axis tight

subplot(2,3,6); 
plot(time_grid,error_canonical_inf2,'k-','linewidth',2);
hold all;
plot(time_grid,error_fjlt_inf2,'b--','linewidth',2);
plot(time_grid,error_gaussian_inf2,'m:','linewidth',2);
plot(time_grid,error_haar_inf2,'c-.','linewidth',2);
legend('canonical','fjlt','gaussian','haar')
xlabel('lags','interpreter','latex')
ylabel('$\frac{\|\mathbf{L}[\Sigma-\hat{\Sigma}]\|_\infty}{\|\Sigma-\hat{\Sigma}\|_2}$','interpreter','latex')
title('inf2')
set(gca,'xscale','log')
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Rotation',90)
axis tight

%% Plot autocorr: check if correct autocorrelation is computed
figure(7); clf;
plot( cc_orig, 'k--','linewidth',10, 'DisplayName', 'Ground Truth');
hold all
plot( cc_canonical, 'linewidth',2, 'DisplayName','Canonical (precond)');
plot( cc_gaussian, 'linewidth',2, 'DisplayName','Gaussian');
plot( cc_fjlt, 'linewidth',2, 'DisplayName','FJLT');
plot( cc_haar, 'linewidth',2, 'DisplayName','Haar');
legend;

%%
time_grid = TMax * percentage_grid;
ubd = sqrt(1+log(time_grid));

figure(1);clf;
plot(time_grid,error_canonical_1F/N,'k-','linewidth',2);
hold all;
plot(time_grid,error_fjlt_1F/N,'b--','linewidth',2);
plot(time_grid,error_gaussian_1F/N,'m:','linewidth',2);
plot(time_grid,error_haar_1F/N,'c-.','linewidth',2);
plot(time_grid,ubd,'r-','linewidth',3)
legend('canonical','fjlt','gaussian','haar','Baseline')
xlabel('Time Size')
ylabel('Relative Error Ratio')
title('1F')
set(gca,'FontSize',25);


%% uniform plot
figure(3);clf;
plot(percentage_grid,error_canonical_1F ,'-s','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[1,0,0]);
hold all;
plot(percentage_grid , error_gaussian_1F ,'-<','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[0.5,1,0.5]);
plot(percentage_grid , error_fjlt_1F ,'-d','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[0,1,0]);
plot(percentage_grid , error_haar_1F ,'-^','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[0,0,1]);
plot(percentage_grid , error_sparse_1F ,'-p','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[1,1,0]);
% plot(sparsity_soft,error_soft,'linewidth',1.5,'color','r');
% plot(sparsity_SCAD,error_SCAD,'linewidth',1.5,'color','g');
title('Error threshold V.S. Subsampling Ratio')
legend('canonical','gaussian','fjlt','haar','sparse')
ylabel('Error')
xlabel('Subsampling Ratio')
set(gca,'fontsize',36)

%% nonuniform plot
figure(3);clf;
plot(error_canonical_1F,percentage_grid ,'-s','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[1,0,0]);
hold all;
plot(error_FD , percentage_grid ,'-<','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[0.5,1,0.5]);
plot(error_AHK06 , subRatio_AHK06 ,'-d','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[0,1,0]);
plot(error_AKL13 , percentage_grid ,'-^','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[0,0,1]);
plot(error_DZ11,subRatio_DZ11 ,'-p','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[1,1,0]);
plot(error_NDT09,subRatio_NDT09 ,'-v','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[0,1,1]);
plot(error_AM01,subRatio_AM01 ,'->','linewidth',2,'MarkerSize',10,...
    'MarkerFaceColor',[1,1,1]);
title('Error threshold V.S. Subsampling Ratio')
legend('canonical','FD','AHK06','AKL13','DZ11','NDT09','AM01');
% legend('canonical','canonicalBias');
legend('DZ11 -> s=100')
xlabel('Error')
ylabel('Subsampling Ratio')
set(gca,'fontsize',36)



% Several subsampling methods can be used to sparsify data. We test w.r.t.
% different subsampling methods the accuracy of autocorrelation when 



%%
Told = 0;
ind = 3;
%     ind = ceil(length(percentage_grid)/2);
T = ceil(TMax*percentage_grid(ind));
if Told == 0
    data = X( 1:T, : );
    Sigma = data*data';
else
    dataNew = X( Told+1:T, : );
    AB = data*dataNew';
    Sigma = [Sigma, AB; AB', dataNew*dataNew'];
    data = [data;dataNew];
    %         if ind==3
    %             % do a check
    %             dataNew = X(1:T,:);
    %             SigmaCheck = dataNew*dataNew';
    %             fprintf(2,'Checking: %g\n', norm( SigmaCheck - Sigma, 'fro' )/norm(SigmaCheck,'fro') );
    %         end
end
Told = T;
cc_orig = autocorr(Sigma,'SCALEOPT', SCALEOPT,...
    'Normalization',NORMALIZATION,'scale',1/N);
%     cc_orig2 = autocorr(Sigma/N,'SCALEOPT', SCALEOPT,...
%         'Normalization',NORMALIZATION,'scale',1); % this is equivalent
%     fprintf('norm(cc-cc2) is %g\n', norm( cc_orig - cc_orig2) );

repeat_ind = 1;

m = ceil( subsample_ratio * N );

fprintf('Time dim size: %u \t Subsample ratio: %.3f \t Repeat: %d \n', ...
    T, subsample_ratio, repeat_ind);

[C_canonical,~,~] = covariance_compute(data, 'subsample_absolute',m,...
    'PRECOND',PRECOND_CANONICAL, 'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','canonical');
fprintf('  Method %s done \n','canonical');

[C_gaussian,~,~] = covariance_compute(data,'subsample_absolute',m,...
    'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','gaussian');
fprintf('  Method %s done \n','gaussian');

%         m = ceil( subsample_ratio * T );
%         scale = 1/m;
scale   = 1/N;

cc_canonical = autocorr(C_canonical,'SCALEOPT', SCALEOPT,'Normalization',NORMALIZATION,'scale',scale);
cc_gaussian  = autocorr(C_gaussian,'SCALEOPT', SCALEOPT,'Normalization',NORMALIZATION,'scale',scale);


figure(2); clf;
plot( cc_orig, '--','linewidth',2, 'DisplayName', 'Ground Truth');
hold all
% plot( cc_canonical,  'linewidth',2, 'DisplayName','Canonical (precond)');
%   something is wrong
plot( cc_gaussian, 'linewidth',2, 'DisplayName','Gaussian');
legend

%% Check Unbiasedness
nReps   = 1e3;
errFcn  = @(CC) norm( CC - Sigma, 'fro' )/norm(Sigma,'fro');
errHist = zeros(nReps,1);
CC      = zeros(size(Sigma));
for rep = 1:nReps
    
%     [C_gaussian,~,~] = covariance_compute(data,'subsample_ratio',subsample_ratio,...
%         'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','gaussian');
    [C_sparse,~] = covariance_compute(data,'subsample_ratio',subsample_ratio,...
    'SUBSAMPLE',SUBSAMPLE,'SUBSAMPLE_METHOD','sparse');
    CC  = CC + C_sparse;
    errHist(rep) = errFcn( CC/rep );
end
%%
figure(3); clf;
semilogx( errHist )
