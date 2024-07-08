function [m1,m2]=getmassestimates(tau0,tau1p5)
    %Function to calculate mass estimates (in Solar Masses) based on a
    %given set of chirp times

    
    %% Constants
%     c = 299792458;
    c = 3*10^8;
    Msolar = 1.989*10^30;
    G = 6.6743*10^-11;
    fmin = 30;
    
    %Calculate Masses
    est_M = (5/(32*fmin))*(tau1p5/(pi*pi*tau0))*(c^3/G);
    est_u = (1/(16*fmin*fmin))*(5/(4*pi^4*tau0*tau1p5^2))^(1/3)*(c^3/G);
    
    m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/(2*Msolar);
    m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/(2*Msolar);
end