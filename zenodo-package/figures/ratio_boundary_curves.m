function [ratio_fig] = ratio_boundary_curves(ratio, linecolor, linesize)
%Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
%Frequency min max values
fmin = 30;
fmax = 700;

% ratio_thresholds = [0.1,0.25,0.5,0.9,1];

%Tau0 values 

tau0s = linspace(0,90,10000);


%Get boundary curve

a = ((5*fmin/8)*(1 + ratio^2))^3;
b = (5*pi^2/4)*tau0s.^2;
tau1p5s = (b./a).^(1/5);

% for i = 1:length(ratio_thresholds)
%     ratio_threshold = ratio_thresholds(i);
%     a = (5*fmin*(1+ratio_threshold^2)/8)^3;
%     b = 5*pi^2*tau0s^2;
%     tau1p5s = (b/a)^(1/5);
%     boundary_curves = [boundary_curves; tau1p5s];
% end

%Plot all boundary curves
% sz = 5;
% c = 'black';
ratio_fig =plot(tau0s, tau1p5s,"Color",linecolor,'LineWidth',linesize,'HandleVisibility','off');
end