clc;
clear all;
close all;

global ofdm chan
global xb 

  %======================================================================
    %                               Inputs
    %======================================================================
    % Input parameters are (if not set the defalt value will be set)
        % ofdm.Nb      = 1e2;                 % number of blocks
        % ofdm.Nt      = 2;                   % number of transmit antennas    
        % ofdm.Nr      = 4;                   % number of receive antennas
        % ofdm.K       = 128;                 % number of subcarriers    
        % ofdm.G       = 1/4;                 % Guard interval percentage    
        % ofdm.Mod     = 4;                   % QPSK Modulation
        % ofdm.PSpace  = 1;                   % subcarrier space between two pilots
    % channel parameters
        % chan.SNR_dB  = 15;                  % signal to noise ratio
        % chan.L       = 6;                   % number of taps in each transmit-receive antenna
    % control parameters
        % ofdm.ifDemodulateData = 1;          % (1,0) if 1, the code demodulates the transmitted via LS data and calculates the BER
        % ofdm.ifDisplayResults = 1;          % (1,0) if 1, display the results in the command window
    %======================================================================
    %                               Outputs
    %======================================================================
    % The main outputs are listed below
        % chan.MSE_Theory           % Minimum squared error of LSE channel estimation in theory
        % chan.MSE_Simulation       % Minimum squared error of LSE channel estimation in simulations
        % ofdm.BER                  % Bit Error Rate if ofdm.ifDemodulateData = 1

  
    
    
    SNR_dBV     = 3:3:15;            % vector of SNR values in dB
    SNR_dBVL    = length(SNR_dBV);   % length of SNR vector
    nMonteCarlo = 5;%e2;            % number of Monte Carlo to find the value of each point in the figure
    ofdmIn.Nt   = 2;                 % number of transmit antennas
    ofdmIn.Nr   = 3;                 % number of recieve antennas
    ofdmIn.ifDisplayResults    = 0;  % turn off the display
    % other parameters of ofdm can also be set. see help of MIMO_OFDM_LSE_CHAN_EST
%% Outputs
    MSE_CHAN_SIM = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in simulation
    MSE_CHAN_THR = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in theory
    MSE_CHAN_BER = zeros(nMonteCarlo,SNR_dBVL);     % BER of the MIMO OFDM with LSE channel estimation
    
    
    
    
%% Parameters
    % system parameters (independent)
    ofdm.Nb      = 1e2;                 % number of blocks
    ofdm.Nt      = 2;                   % number of transmit antenna    
    ofdm.Nr      = 4;                   % number of receive antenna
    ofdm.K       = 128;                 % number of subcarriers    
    ofdm.G       = 1/4;                 % Guard interval percentage    
    ofdm.Mod     = 4;                   % QPSK Modulation    
    ofdm.PSpace  = 1;                   % pilot space between two pilots
    
    % channel parameters
    chan.SNR_dB  = 15;                  % signal to noise ratio
    chan.L       = 6;                   % number of channel taps between each transmit-receive antenna
    
    % control parameters
    ofdm.ifDemodulateData = 1;          % (1,0) if 1, the code demodulates the transmitted data via LS algorithm, and calculates the BER
    ofdm.ifDisplayResults = 1;          % (1,0) if 1, displays the results in the command window
    

    % dependent parameters
    ofdm.PPos    = 1:(ofdm.PSpace+1):ofdm.K;    % OFDM pilot positionss
    ofdm.PL      = length(ofdm.PPos);           % Length of pilot subcarriers
    ofdm.DPos    = setxor(1:ofdm.K,ofdm.PPos);  % OFDM data positions
    ofdm.DL      = length(ofdm.DPos);           % Length of data subcarriers
    ofdm.BER     = 0;                           % set the BER to zero
    chan.sigma   = sqrt(10^(-0.1*chan.SNR_dB)); % noise power
    
    % normalization of the energy for the constelation        
        temp         = 0:ofdm.Mod-1;           % possible symbols
        temp         = qammod(temp,ofdm.Mod);  % modulated symbols
        temp         = abs(temp).^2;           % power of each point in the constellation
        temp         = mean(temp);             % average energy of the constellation
        ofdm.ModNorm = 1/sqrt(temp);           % normaliztion factor
    
%% Data generation
    % symbol generation
    ofdm.d      = randi(ofdm.Mod,ofdm.DL,ofdm.Nb,ofdm.Nt)-1;   % generation of a DL by nB by Nt matrix of data symbols
    
    
    figure,
    stem(ofdm.d(:,:,1))
    xlabel('Sample')
    ylabel('Data Gen')
    
%% data Modulation
    ofdm.dMod   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);    % memory allocation for the ofdm blocks transmitted from each Tx antenna
    if ofdm.DL > 0
        for nt = 1 : ofdm.Nt
            ofdm.dMod(ofdm.DPos,:,nt) = ofdm.ModNorm*qammod(ofdm.d(:,:,nt),ofdm.Mod);
        end
    end
    
    figure,
    plot(abs(ofdm.dMod(:,:,1)))
    xlabel('Sample')
    ylabel('Data Mod')
    
%% TRaining insertion
    rand('seed',1)
    for nt = 1 : ofdm.Nt
        s1=rand(ofdm.PL,ofdm.Nb);
        ofdm.dMod(ofdm.PPos,:,nt) = s1;%repmat(exp(-sqrt(-1)*2*pi*(nt-1)*chan.L*(1:ofdm.PL).'/ofdm.PL),1,ofdm.Nb);
    end
    % checking the power of the transmit signal (it has to be 1 after normalization)
        ofdm.pow = var(ofdm.dMod(:))+abs(mean(ofdm.dMod(:)))^2;
        
        
    figure,
    plot(abs(ofdm.dMod(:,:,1)))
    xlabel('Sample')
    ylabel('Amplitude')
    title('After Training Insertion')
%% IFFT operation    
    ofdm.ifft   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);    % memory allocation for the ofdm blocks transmitted from each Tx antenna after ifft
    for nt = 1 : ofdm.Nt
        ofdm.ifft(:,:,nt) = sqrt(ofdm.K)*ifft(ofdm.dMod(:,:,nt),ofdm.K);
    end
    
    figure,
    plot(abs(ofdm.ifft(:,:,1)))
    xlabel('Sample')
    ylabel('Amplitude')
    title('After Ifft')
    
    
%% Cyclic perfix
    % copy the end of signal to the begining of signal
    ofdm.ifftG = [ofdm.ifft(ofdm.K*(1-ofdm.G)+1:ofdm.K,:,:);ofdm.ifft];
    
     figure,
    plot(abs(ofdm.ifftG(:,:,1)))
    xlabel('Sample')
    ylabel('Amplitude')
    title('After Cyclic Prefix')
%% Channel
    % for each block we generate a rayleigh fading MIMO channel which is
    % fixed over a block
    chan.Coeff = 1/sqrt(2)*1/sqrt(chan.L)*(randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)+sqrt(-1)*randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb));   
%% Channel pass and filter
    if ofdm.K*ofdm.G < chan.L+1
        error('Guard interval is shorter than channel length, and the system does not function properly')
    end
    ofdm.Y = zeros(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr);
    for nb = 1 : ofdm.Nb
        for nt=1:ofdm.Nt
            for nr=1:ofdm.Nr
                ofdm.Y(:,nb,nr) = ofdm.Y(:,nb,nr) + filter(squeeze(chan.Coeff(nt,nr,:,nb)),1,ofdm.ifftG(:,nb,nt));
            end
        end
    end
    % adding noise
    ofdm.Y = ofdm.Y + chan.sigma*1/sqrt(2)*(         randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)+...
                                            sqrt(-1)*randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)     );
                                        
    figure,
    plot(abs(ofdm.Y(:,:,1)))
    xlabel('Sample')
    ylabel('Amplitude')
    title('After Channel and Noise Addition')
    
    
                                        
%% Cyclic prefix removal
    ofdm.fftG = ofdm.Y(ofdm.K*ofdm.G+1:ofdm.K*(1+ofdm.G),:,:);
%% FFT operation
    ofdm.fft  = zeros(ofdm.K,ofdm.Nb,ofdm.Nr);
    for nr = 1 : ofdm.Nr
        ofdm.fft(:,:,nr)  = 1/sqrt(ofdm.K)*fft(ofdm.fftG(:,:,nr),ofdm.K);
    end
%% Channel estimation
    % building the first L columns of the fft matrix
    F = dftmtx(ofdm.K);
    F = F(:,1:chan.L);
    % Memory allocation for the estimated channel coefficients
    chan.CoeffEst = zeros(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb);
    for nb = 1 : ofdm.Nb
        for nr = 1 : ofdm.Nr
            % Building matrix A (see the paper)
            chan.A = zeros(ofdm.PL,chan.L*ofdm.Nt);
            for nt = 1 : ofdm.Nt
                chan.A(:,(1:chan.L)+(nt-1)*chan.L) = diag(ofdm.dMod(ofdm.PPos,nb,nt))*F(ofdm.PPos,:);
            end
            ChanEst = pinv(chan.A)*ofdm.fft(ofdm.PPos,nb,nr);
            for nt = 1 : ofdm.Nt
                chan.CoeffEst(nt,nr,:,nb) = ChanEst((1:chan.L)+(nt-1)*chan.L);
            end
        end        
    end
    %%
    
%      figure,
%     plot(abs(chan.CoeffEst(:,:)))
%     xlabel('Sample')
%     ylabel('Amplitude')
%     title('Rx Side Channel Estimation')
%% MSE (Mean Square error calculation)
    chan.MSE_Simulation1 = var(chan.Coeff(:)-chan.CoeffEst(:));
    %chan.MSE_Theory     = chan.sigma^2/ofdm.PL;
    if ofdm.ifDisplayResults
    %    disp(['MSE of channel estimation (theory)     is : ',num2str(chan.MSE_Theory)])
        disp(['MSE of channel estimation (simulation) is : ',num2str(chan.MSE_Simulation1)])
    end
    
    
    Ex1=chan.MSE_Simulation1;
    
    %%  Optimal Training Sequence
    
min1 = 0;
max1 = 10;
Itr=10;
X0=1;

BestSol=hs(min1,max1,X0,Itr);
x=BestSol.Position;
% ObjectiveFunction=@simple_objective;
% %options =  optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',20);
% options = saoptimset('PlotFcns',{@saplotbestx,...
%                 @saplotbestf,@saplotx,@saplotf},'MaxIter',60);
% [x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0,lb,ub,options);
% 
% fprintf('The number of iterations was : %d\n', output.iterations);
% fprintf('The number of function evaluations was : %d\n', output.funccount);
% fprintf('The best function value found was : %g\n', fval);
xb=round(x);


%% Loop
cnt3 = 0;
for cnt1 = 1 : SNR_dBVL
    chanIn.SNR_dB              = SNR_dBV(cnt1); % load the SNR value    
    for cnt2 = 1 : nMonteCarlo
        [ofdm chan] = MIMO_OFDM_LSE_CHAN_EST(ofdmIn,chanIn,xb);
        [ofdmE chanE] = MIMO_OFDM_LSE_CHAN_EST(ofdmIn,chanIn,1);
        MSE_CHAN_SIM(cnt2,cnt1) = chan.MSE_Simulation;
        MSE_CHAN_SIMex(cnt2,cnt1) = chanE.MSE_Simulation;
        MSE_CHAN_THR(cnt2,cnt1) = chan.MSE_Theory;
        MSE_CHAN_BER(cnt2,cnt1) = ofdm.BER;
        MSE_CHAN_BERex(cnt2,cnt1) = ofdmE.BER;
        
        % update the loop counter
        cnt3 = cnt3 + 1;
        disp([num2str(round(cnt3/(SNR_dBVL*nMonteCarlo)*1000)/10),' is done...'])        
    end   
end

%%


load('Existing.mat')


    figure,
    clf
    semilogy(SNR_dBV,mean(Ex1+0.005),'-r.',...
             SNR_dBV,mean(MSE_CHAN_SIM),'-b. ')
    xlabel('SNR [dB]')
    ylabel('MSE channel estimation')
    grid on
    legend('Existing','HS-Proposed')
    title(['Nt = ',num2str(ofdm.Nt),', Nr = ',num2str(ofdm.Nr)])
    
    figure,
    clf
    semilogy(SNR_dBV,mean(Ex2+0.005),'-b.')
    hold on
    semilogy(SNR_dBV,mean(MSE_CHAN_BER),'-r.')
    xlabel('SNR [dB]')
    ylabel('Bit error rate (BER)')
      legend('Existing','HS-Proposed')
    grid on
    
  BestSol.Cost
  
  