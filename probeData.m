function data= probeData(varargin)
%Function plots raw data information: time domain plot, a frequency domain
%plot and a histogram. 
%
%The function can be called in two ways:
%   probeData(settings)
% or
%   probeData(fileName, settings)
%
%   Inputs:
%       fileName        - name of the data file. File name is read from
%                       settings if parameter fileName is not provided.
%
%       settings        - receiver settings. Type of data file, sampling
%                       frequency and the default filename are specified
%                       here. 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

% CVS record:
% $Id: probeData.m,v 1.1.2.7 2006/08/22 13:46:00 dpl Exp $

    
%% Generate plot of raw data ==============================================

if (nargin == 1)
    settings = deal(varargin{1});

else
    error('Incorect number of arguments');
end

    samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
      
                        % Read data for acquisition. 11ms of signal are needed for the
        % fine   ´Ó²¶»ñÖÐ¶ÁÈ¡Êý¾Ý11ms
        % frequency estimation
        %´ÓÎÄ¼þÖÐ¶ÁÈ¡ÐÅºÅ data = fread(fid, 11*samplesPerCode, settings.dataType)';
        %Generate the carrier frequency ÓÃµ¥ÔØ²¨ÐÅºÅÈ¡´úÔ­ÎÄ¼þÐÅºÅ -----------
           
      %%Prepare for generating signals  
     
            time=settings.msToProcess;  %1000ms   
    carrFreq_rcv= 15.58e6; %15.580030e6;
            blksize_rcv=time*samplesPerCode;  %about 1000ms's signal,there are 1000 c/a periods 50 base periods
            time_rcv    = (0:blksize_rcv) ./ settings.samplingFreq;
           % caCodeForSamples=zeros(1,samplesPerCode);
            
       
            
            %% Get the argument to sin/cos functions phase=0 ÏàÎ»Îª0
            trigarg_rcv = ((carrFreq_rcv * 2.0 * pi) .* time_rcv);
            
            % Finally compute the signal to mix the collected data to bandband
          %  data_0 =zeros(1,blksize_rcv,'int8');
            Am=sqrt(1);
            
            data_0 =Am* sin(trigarg_rcv(1:blksize_rcv));
            
  %--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % C/A chip period in sec
    %--- Generate CA code for given PRN -----------------------------------
  % PRN=25;
 %   caCode = generateCAcode(PRN);
 
    %=== Digitizing =======================================================
    
    %--- Make index array to read C/A code values -------------------------
    % The length of the index array depends on the sampling frequency -
    % number of samples per millisecond (because one C/A code period is one
    % millisecond).
    %codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);
    
    %--- Correct the last index (due to number rounding issues) -----------
    %codeValueIndex(end) = 1023;
    
    %--- Make the digitized version of the C/A code -----------------------
    % The "upsampled" code is made by selecting values form the CA code
    % chip array (caCode) for the time instances of each sample.

    %caCodeForSamples(:) = caCode(codeValueIndex);

    %multiple cacode and the carrier to generate the input--------------
           % samplesPerChip0=round(settings.samplingFreq/settings.codeFreqBasis);
           % samplesPerChip1=samplesPerChip0-1;
           %data_1=zeros(1,blksize_rcv,'int8');
            
           % for i=1:time
            %    data_1(((i-1)*samplesPerCode+1):(i*samplesPerCode))=data_0(((i-1)*samplesPerCode+1):(i*samplesPerCode)).*caCodeForSamples;
                      
            %end
  

       % probe message bit signals------------------------------------------------
        %data=zeros(1,blksize_rcv,'int8');
       %%     make a 62000*time signals
          for j=1:(time/20)
                    if mod(j,2)
                   data(((j-1)*20*samplesPerCode+1):j*20*samplesPerCode)= data_0(((j-1)*20*samplesPerCode+1):j*20*samplesPerCode);
                    else
                    data(((j-1)*20*samplesPerCode+1):j*20*samplesPerCode)= -data_0(((j-1)*20*samplesPerCode+1):j*20*samplesPerCode);
                    end
          end
           
          %add Gaussian noise
          data=awgn(data,10^(-5),'measured');    %,0.1,10^(-15),'linear');
    %--- Initialization ---------------------------------------------------
    figure(100);
    clf(100);
    
    timeScale = 0 : 1/settings.samplingFreq : 5e-3;    
    
    %--- Time domain plot -------------------------------------------------
    subplot(2, 2, 1);
    plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
         data(1:round(samplesPerCode/50)));
     
    axis tight;
    grid on;
    title ('Time domain plot');
    xlabel('Time (ms)'); ylabel('Amplitude');
    
    %--- Frequency domain plot --------------------------------------------
    subplot(2,2,2);
    pwelch(data-mean(data), 16384, 1024, 2048, settings.samplingFreq/1e6)
    
    axis tight;
    grid on;
    title ('Frequency domain plot');
    xlabel('Frequency (MHz)'); ylabel('Magnitude');
    
    %--- Histogram --------------------------------------------------------
    subplot(2, 2, 3.5);
    hist(data, -128:128)
    
    dmax = max(abs(data)) + 1;
    axis tight;
    adata = axis;
    axis([-dmax dmax adata(3) adata(4)]);
    grid on;
    title ('Histogram'); 
    xlabel('Bin'); ylabel('Number in bin');

end % if (fid > 0)
