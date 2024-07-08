function [PSD] = getsegPSD(filtsegment, trainidxs,sampFreq)
%GETSEGPSD Generate Two-sided PSD from segment

%Load training segment
Seglen = length(filtsegment)/sampFreq;
trainseg = filtsegment(trainidxs(1):trainidxs(2));

%Welch's estimate (One-sided)
[pxx,f]=pwelch(trainseg, tukeywin(0.5*sampFreq),[],[],sampFreq);

%Interpolate PSD
[interPSD, ~] = createPSD(sampFreq, Seglen, pxx', f');

%Make two-sided
PSD = interPSD./2;
end

