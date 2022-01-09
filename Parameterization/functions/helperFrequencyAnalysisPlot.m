function helperFrequencyAnalysisPlot(F,Ymag,Yangle,NFFT,ttlMag,ttlPhase)
% Plot helper function for the FrequencyAnalysisExample
% F: Frequency []
% Ymag: Magnitude of Y []
% Yangle: Phase of Y
% NFFT: Length of Y [-]

% Copyright 2012 The MathWorks, Inc.

figure
subplot(2,1,1)
ind = fix(NFFT/2) + 1; 
plot(F(1:ind), 20*log10(Ymag(1:ind)));

if nargin > 4 && ~isempty(ttlMag)
  tstr = {'Magnitude response of the audio signal',ttlMag};
else
  tstr = {'Magnitude response of the audio signal'};
end

title(tstr)
xlabel('Frequency in Hz')
ylabel('dB')
grid on;
axis tight 


subplot(2,1,2)
plot(F(1:ind),Yangle(1:ind));
if nargin > 5
  tstr = {'Phase response of the audio signal',ttlPhase};
else  
  tstr = {'Phase response of the audio signal'};
end

title(tstr)
xlabel('Frequency in Hz')
ylabel('radians')
grid on;
axis tight