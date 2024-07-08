function [glitch,varargout]=GRITCLEAN(posparams,negparams, thresholds)
    %Function to run veto on a set of positive and negative parameters
    %INPUTS: posparams: MX4 candidate cell array with each cell having a 1X4 array
    %                   containing positive run parameter estimates
    %                   [tau0_pos, tau1.5_pos, TOA_pos, SNR_pos]
    %        negparams: MX4 candidate cell array with each cell having a 1X4 array
    %                   containing negative run parameter estimates
    %                   [tau0_neg, tau1.5_neg, TOA_neg, SNR_neg]
    %       thresholds: 1X3 array containing complex-mass ratio,
    %                   |Delta t_a| and |Delta \rho| thresholds for the S_N veto steps 

    %OUTPUT: glitch: MX2 array of 1s and 0s indicating whether candidates
    %                are glitches (1s) or not (0s) and what type of veto they failed
    %                Each row would be [a,b] where a = 0 or 1 and b = 1,2,3
    %                or 0 to denote the number of the veto step that caught the
    %                glitch, 0 means no step catches it and it is
    %                classifies as not a glitch
    %        absdeltatoa: Optional Output, returns a 1XM array of absolute
    %                     TOA differences
    %        absdeltasnr: Optional Output, returns a 1XM array of absolute
    %                     relative SNR differences 
    %               snrs: Optional Output, returns a 1XM array of SNRs of
    %                     all candidates
    %-----------DO NOT CHANGE BELOW------------------------

    exception = MException('myComponent:inputError', 'Lengths of positive and negative run candidates do not match');
    if length(posparams) ~= length(negparams)
        throw(exception);
    end
    
    %Load Thresholds
    zeta = thresholds(1);
    delta_ta_threshold = thresholds(2);
    delta_snr_threshold = thresholds(3);

    glitch = [];
    absdeltatoa = [];
    absdeltasnr = [];
    snrs = [];
    for i = 1:length(posparams)
        params_positive = posparams{i};
        params_negative = negparams{i};

        %Run S_P veto steps
        % 1.) Check Chirp length
        tau0_pos = params_positive(1);
        tau1p5_pos = params_positive(2);
        TOA_pos = params_positive(3);
        SNR_pos = params_positive(4);

        TOA_neg = params_negative(3);
        SNR_neg = params_negative(4);
        
        delta_ta = abs(TOA_pos - TOA_neg);
        delta_snr = abs((SNR_pos - SNR_neg)/SNR_pos);
        absdeltatoa = [absdeltatoa,delta_ta];
        absdeltasnr = [absdeltasnr,delta_snr];
        snrs = [snrs; [SNR_pos,SNR_neg]];
        chirplength = getchirplength(tau0_pos,tau1p5_pos);

        

        if chirplength < 0
            glitch = [glitch; [1,1]];
            continue;
        end

        % 2.) Check Complex Mass ratio
        [m1,m2] = getmassestimates(tau0_pos, tau1p5_pos);
        ratio1 = abs(imag(m1)/real(m1)); ratio2 = abs(imag(m2)/real(m2));
        if ratio1 >= zeta || ratio2 >= zeta
            glitch = [glitch; [1,2]];
            continue;
        end

        % 3.) Run S_N veto step
        if delta_ta < delta_ta_threshold || delta_snr < delta_snr_threshold
            glitch = [glitch; [1,3]];
            continue;
        end

        %If all veto steps are bypassed return [0,0]
        glitch = [glitch; [0,0]];
    end

    if nargout > 1
        varargout{1}=absdeltatoa;
        varargout{2}=absdeltasnr;
        varargout{3}=snrs;
    end

end