 function fit = simple_objective(x)
   
%  global nRx N1 h sCode msg rayleighChan snrIndB hStr hDemod idx N Order r hMod1 
%   global ZFBERCalc MMSEBERCalc BERCalc MLBERCalc rxSig snrLinear allTxSig allBits
%   global EbNo
%   

global ofdm chan 
 
%% TRaining insertion
    
    rand('seed',round(x))
    for nt = 1 : ofdm.Nt
        s1=rand(ofdm.PL,ofdm.Nb);
        ofdm.dMod(ofdm.PPos,:,nt) = s1;%repmat(exp(-sqrt(-1)*2*pi*(nt-1)*chan.L*(1:ofdm.PL).'/ofdm.PL),1,ofdm.Nb);
    end
    % checking the power of the transmit signal (it has to be 1 after normalization)
        ofdm.pow = var(ofdm.dMod(:))+abs(mean(ofdm.dMod(:)))^2;
%% IFFT operation    
    ofdm.ifft   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);    % memory allocation for the ofdm blocks transmitted from each Tx antenna after ifft
    for nt = 1 : ofdm.Nt
        ofdm.ifft(:,:,nt) = sqrt(ofdm.K)*ifft(ofdm.dMod(:,:,nt),ofdm.K);
    end
%% Cyclic perfix
    % copy the end of signal to the begining of signal
    ofdm.ifftG = [ofdm.ifft(ofdm.K*(1-ofdm.G)+1:ofdm.K,:,:);ofdm.ifft];
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
%% MSE (Mean Square error calculation)
    chan.MSE_Simulation1 = var(chan.Coeff(:)-chan.CoeffEst(:));
    %chan.MSE_Theory     = chan.sigma^2/ofdm.PL;
    if ofdm.ifDisplayResults
    %    disp(['MSE of channel estimation (theory)     is : ',num2str(chan.MSE_Theory)])
        disp(['MSE of channel estimation (simulation) is : ',num2str(chan.MSE_Simulation1)])
    end
    
 
 fit=chan.MSE_Simulation1;
 
 
% 
% figure,
% plot(EbNo, rZF, 'r-', ...
%          EbNo, rMMSE, 'b-', ...
%          EbNo, rML, 'g-',...
%          EbNo, rSTBC,'c-');
% xlabel('Eb/No')
% ylabel('SE')
% legend('ZF-SIC', 'MMSE-SIC', 'ML','MRC');
 
 end
  