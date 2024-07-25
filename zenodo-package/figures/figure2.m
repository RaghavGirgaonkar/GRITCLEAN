%Get image matrix
S = load('SP_imagvals.mat');
Fmatrix = S.fitvals;
% Fmatrix = Fmatrix./(max(Fmatrix(:)));
% Fmatrix(Fmatrix > 0) = 1;

%Make Plot
sigcolor = 'blue';
psocolor = 'black';
x = linspace(0, 90, 10000);
y = linspace(0, 2, 10000);
% x = linspace(0.00001, 90, 10000);
% y = linspace(0.00001, 2, 10000);
% x = linspace(0, 3, 10000);
% y = linspace(0, 2, 10000);
% x = log10(linspace(0.00001, 90, 10000));
% y = log10(linspace(0.00001, 2, 10000));
% x = linspace(-2, 0, 10000);
% y = linspace(-2, 0, 10000);

figure; hold on;
% xlim([0,90]);
% ylim([0,2]);
map = [153 204 255;
           255 255 133]./255;
% map = [255 255 133;
%     153 204 255]./255;
% map = [255 255 133;
%     255 255 255;
%     153 204 255;
%     255 204 255;
%        204 204 255;
%        153 255 153;
%        178 255 102;
%        255 255 51;
%        255 128 0;
%        204 0 0]./255;

% map = [255 255 255
%     204 229 255;
%     153 204 255;
%     102 178 255;
%     51 153 255;
%     0 128 255;
%     0 102 204;
%     0 76 153;
%     153 153 255;
%     102 102 255;
%     ]./255;

ratio_thresholds = [0,0.25,0.5,0.9];
colormap(map);
imagesc(x,y, Fmatrix'); axis xy;
hold on;
boundary_plot;
hold on;
for i = 1:length(ratio_thresholds)
    if ratio_thresholds(i) == 0
        ratio_boundary_curves(ratio_thresholds(i),'red',3);
        hold on;
    else
        ratio_boundary_curves(ratio_thresholds(i),'black',1);
        hold on;
    end
end
hold on;
scatter(0.6,1.8,70,'o','MarkerEdgeColor','black','LineWidth',1.5);
% hold on;
% scatter(0.001,1.4,70,'o','MarkerEdgeColor','black','LineWidth',1.5);
% for i = 1:203
%     scatter(chirptimes(i,1),chirptimes(i,2),70,'o','MarkerEdgeColor','black','LineWidth',1.5);
% end
% scatter(tau0, tau1p5, snrs + 30, psocolor,'DisplayName','PSO Estimated Locations');
% hold on;
% scatter(injectedtau0, injectedtau1p5, 70,'x',sigcolor,'DisplayName','Injected Signal Locations');
% xlim([-2,0]);
% ylim([-2,0]);
% xlim([0,90]);
% % ylim([0,2]);
xlim([0.0001,90]);
ylim([0.0001,2]);
% xlim([0,3]);
% ylim([0,2]);
% xlim(log10([0.0001,90]));
% ylim(log10([0.0001,2]));
hold off;
xlabel('\tau_0');
ylabel('\tau_{1.5}');
ax = gca;
ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
% legend('FontSize',20);
