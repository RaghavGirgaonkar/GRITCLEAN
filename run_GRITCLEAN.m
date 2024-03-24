function [glitch_status,candidate_idxs,abstoadiffs, abssnrdiffs,snrs,noise_positive_candidates,noise_negative_candidates] = run_GRITCLEAN(posfile,negfile,thresholds)
%RUN_GRITCLEAN Load Files and Thresholds and run GRITCLEAN
%   %-----------DO NOT CHANGE BELOW------------------------

positive_candidates = {};
negative_candidates = {};

noise_positive_candidates = [];
noise_negative_candidates = [];

poscands = textread(posfile,'%s','delimiter', '\n');
negcands = textread(negfile,'%s','delimiter', '\n');

%Gather candidates above detection thresholds
detection_threshold = thresholds(1);
param_thresholds = thresholds(2:end);

candidate_idxs = [];
noise_idxs = [];
for i = 1:length(poscands)
    posvals = str2num(poscands{i});
    negvals = str2num(negcands{i});

    %Ignore all candidates where the signal lies inside the overlap
    %region
    chirplength = getchirplength(posvals(2),posvals(3));
    if posvals(end-1) + chirplength > 512 - 64
        continue;
    end
    
    if posvals(end) > detection_threshold
        candidate_idxs = [candidate_idxs,posvals(1)];
        positive_candidates = [positive_candidates; posvals(2:end)];
        negative_candidates = [negative_candidates; negvals(2:end)];
    else
        noise_idxs = [noise_idxs, posvals(1)];
        noise_positive_candidates = [noise_positive_candidates;posvals(2:end)];
        noise_negative_candidates = [noise_negative_candidates;negvals(2:end)];
        continue;
    end
end

%Apply GRITCLEAN
[glitch_status, abstoadiffs, abssnrdiffs,snrs] = GRITCLEAN(positive_candidates, negative_candidates, param_thresholds);

end

