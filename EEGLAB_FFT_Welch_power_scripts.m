%%  Perform the Fourier transform using Welch's method on multiple data stored in ALLEEG structure (EEGLAB)

% If individual data have different sampling rates, they should be all
% re-calculated to the same level, or the value of sRate/pNTS variables
% should be changed

sRate = ALLEEG(1).srate; 
pNTS = ALLEEG(1).pnts;

% 75% window overlapping, can be changed
winNoverlap = (pNTS*75)/100;

% pNNFT defines the number of discrete Fourier transform (DFT) points in the PSD estimation.
% The default value is a multiple of 256, or the next higher power of 2
% greater than the number representing the length of the segment to be analyzed.
% For example, for an EEG with a sampling rate of 500 Hz, pNFFT = 512.
   
pNFFT = 512;

h = waitbar(0,'Calculation started...');
waitX = size(ALLEEG,2);

      for ii = 1:size(ALLEEG,2)
         clear spectra freqs 
            for kk = 1:size(ALLEEG(ii).data,1)
                [FFT(ii).data(:,kk), freqs] = pwelch(ALLEEG(ii).data(kk,:), pNTS, winNoverlap, pNFFT, sRate,'psd');
            end
      waitbar(ii/waitX,h,sprintf('%d/%d datasets analyzed',ii, waitX))
      end
      
      close(h);
      
clear sRate pNTS winNoverlap pNFFT waitX ii jj kk h

%% An example image showing what the FFT looks like # here for dataset no. 1

for ii = 1:size(FFT,2)
      FFT(ii).MeanE=mean(FFT(ii).data,2);   
end 

figure;
hold on;
plot(freqs(1:size(freqs,1)),10*log10(FFT(1).MeanE(1:size(freqs,1))),'LineWidth', 3, 'Color', [1 0 0]); 
xlim([1 40]) % limit osi X - 0-40 Hz
set(gca,'fontsize',14)
grid on, grid minor, box off
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]')
hold off;

%% Bandpower: the average power of a signal in a specific frequency range

% Loop through all data stored in FFT
for ii = 1:size(FFT,2)

    % Loop through all electrodes
    for kk = 1:size(FFT(ii).data,2)

        % Total power of the signal needed for relative power calculation
        % (here for 1-30 Hz)   
        FFT(ii).sumPower(:,kk) = bandpower(FFT(ii).data(:,kk),freqs,[1 30],'psd');

        % Absolute power calculation for delta [1-3 Hz], theta [4-7 Hz], alpha [8-12 Hz] and beta [13-30 Hz] bands
        FFT(ii).delta(:,kk) = bandpower(FFT(ii).data(:,kk),freqs,[1 3],'psd');
        FFT(ii).theta(:,kk) = bandpower(FFT(ii).data(:,kk),freqs,[4 7],'psd');
        FFT(ii).alpha(:,kk) = bandpower(FFT(ii).data(:,kk),freqs,[8 12],'psd');
        FFT(ii).beta(:,kk) = bandpower(FFT(ii).data(:,kk),freqs,[13 30],'psd');
        
        % Relative power calculation
        FFT(ii).rel_delta(:,kk) = FFT(ii).delta(:,kk) / FFT(ii).sumPower(:,kk);
        FFT(ii).rel_theta(:,kk) = FFT(ii).theta(:,kk) / FFT(ii).sumPower(:,kk);
        FFT(ii).rel_alpha(:,kk) = FFT(ii).alpha(:,kk) / FFT(ii).sumPower(:,kk);
        FFT(ii).rel_beta(:,kk) = FFT(ii).beta(:,kk) / FFT(ii).sumPower(:,kk);           
             
    end
    
    FFT(ii).MeanDelta=mean(FFT(ii).delta);   
    FFT(ii).MeanTheta=mean(FFT(ii).theta); 
    FFT(ii).MeanAlpha=mean(FFT(ii).alpha);     
    FFT(ii).MeanBeta=mean(FFT(ii).beta);

    FFT(ii).Mean_rel_Delta=mean(FFT(ii).rel_delta);   
    FFT(ii).Mean_rel_Theta=mean(FFT(ii).rel_theta); 
    FFT(ii).Mean_rel_Alpha=mean(FFT(ii).rel_alpha);     
    FFT(ii).Mean_rel_Beta=mean(FFT(ii).rel_beta); 
  
       
end

%% An example of calculation of a band ratio, here alpha/theta

for ii = 1:size(FFT,2)
     FFT(ii).ratio_delta_theta_abs = FFT(ii).alpha ./ FFT(ii).theta;
end
