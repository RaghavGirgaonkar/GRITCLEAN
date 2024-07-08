function []=rungwpsosearch(segnum, runtype)
    %Function to run PSO-based matched-filtering search on a given segment
    %INPUT: segnum: Segment number
    %       runtype: 'P' for a run in S_P 'N' for a run in S_N

    
    %Load Segment data
    segment = loadsegment(segnum);

    %Load PSD file and PSD for given segment
    S = load('GVSsegPSDtrainidxs.mat');
    trainidxs = S.trainidxs{segnum}; %Load training indices for segment
    
    %High-Pass Filter segment data above fmin Hz
    rolloff = 0.125;
    fmin = 30;
    sampFreq = 4096;

    %Window and Highpass filter
    segwin = segment.*tukeywin(length(segment),rolloff*sampFreq/length(segment))';
    seghpass = highpass(segwin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    filtsegment = seghpass;
    Seglen = length(segment)/sampFreq; %Get segment length

    %Estimate two-sided PSD from filtered segment
    PSD = getsegPSD(filtsegment, trainidxs,sampFreq);

    %Whiten and Normalize segment using PSD
    [whtndseg, TFtotal] = segdatacond(filtsegment, PSD, sampFreq);

    
    %Check run type and set appropriate json file
    if runtype == 'P'
        paramfile = 'allparamfiles.json';
    else
        paramfile = 'allparamfiles_negative.json';
    end

    %Run PSO-based matched-filtering search
    runpso(segment, whtndseg, TFtotal, paramfile,segnum,Seglen);

end