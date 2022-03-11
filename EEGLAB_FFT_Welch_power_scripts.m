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
waitX = length(ALLEEG);

      for ii = 1:length(ALLEEG)
         clear spectra freqs 
         for jj = 1:size(ALLEEG(ii).data,3)
            for kk = 1:size(ALLEEG(ii).data,1)
                [FFTResults(ii).data(kk,:,jj), freqs] = pwelch(ALLEEG(ii).data(kk,:,jj), pNTS, winNoverlap, pNFFT, sRate,'psd');
            end
         end
      waitbar(ii/waitX,h,sprintf('%d/%d datasets analyzed',ii, waitX))
      end
      
      close(h);
      
clear sRate pNTS winNoverlap pNFFT waitX ii jj kk h
%% Additional calculations that create a simple table of results ready for statistical analysis

% The average of epochs for each dataset
for i = 1:length(FFTResults)
      FFTResults(i).meanEpochs=mean(FFTResults(i).data,3);   
end  

% The average of all electrodes
for i = 1:length(FFTResults)
      FFTResults(i).meanEl=mean(FFTResults(i).meanEpochs,1);   
end 