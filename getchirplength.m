function [chirplength] = getchirplength(tau0,tau1p5)
%GETCHIRPLENGTH Function to calculate chirp-length based on a given set of
%chirp times

    %% Constants
    c = 3*10^8;
    Msolar = 1.989*10^30;
    fmin = 30;
    G = 6.6743*10^-11;
    M = (5/(32*fmin))*(tau1p5/(pi*pi*tau0))*(c^3/G);

    u = 1/(16*fmin*fmin)*(5/(4*pi^4*abs(tau0)*abs(tau1p5)^2))^(1/3)*(c^3/G);
    n = u/M;
    %% Calculate Chirp Times
    tau1 = (5/(192*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1))*(1/n)*((743/336)+ (11*n/4));    
    tau2 = (5/(128*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1/3))*(1/n)*((3058673/1016064) + (5429*n/1008) + (617*n*n/144));

    chirplength = tau0 - tau1p5 + tau1 + tau2;
end

