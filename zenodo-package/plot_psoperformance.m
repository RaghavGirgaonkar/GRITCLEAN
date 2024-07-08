realsigfile = 'Output/psoperf_injsigs_realsnrs.txt';
lowsigfile = 'Output/psoperf_injsigs_lowsnrs.txt';
highsigfile = 'Output/psoperf_injsigs_highsnrs.txt';
massgapsigfile = 'Output/psoperf_massgap_snrs.txt';
highmassrealsigfile = 'Output/psoperf_25to40Msun_realsnrs.txt';
highmasslowsigfile = 'Output/psoperf_25to40Msun_lowsnrs.txt';
highmasshighsigfile = 'Output/psoperf_25to40Msun_highsnrs.txt';

realcandidates = textread(realsigfile, '%s', 'delimiter', '\n');
lowcandidates = textread(lowsigfile, '%s', 'delimiter', '\n');
highcandidates = textread(highsigfile, '%s', 'delimiter', '\n');
massgapcandidates = textread(massgapsigfile, '%s', 'delimiter', '\n');
highmassrealcandidates = textread(highmassrealsigfile, '%s', 'delimiter', '\n');
highmasslowcandidates = textread(highmasslowsigfile, '%s', 'delimiter', '\n');
highmasshighcandidates = textread(highmasshighsigfile, '%s', 'delimiter', '\n');

realogfitvals = [];
realpsofitvals = [];
lowogfitvals = [];
lowpsofitvals = [];
highogfitvals = [];
highpsofitvals = [];
massgapogfitvals = [];
massgappsofitvals = [];
highmassrealogfitvals = [];
highmassrealpsofitvals = [];
highmasslowogfitvals = [];
highmasslowpsofitvals = [];
highmasshighogfitvals = [];
highmasshighpsofitvals = [];

for i = 1:length(realcandidates)
    realvals = str2num(realcandidates{i});
    realogfitvals = [realogfitvals,realvals(2)];
    realpsofitvals = [realpsofitvals,realvals(3)];
end

for i = 1:length(lowcandidates)
    lowvals = str2num(lowcandidates{i});
    lowogfitvals = [lowogfitvals,lowvals(2)];
    lowpsofitvals = [lowpsofitvals,lowvals(3)];
end

for i = 1:length(highcandidates)
    highvals = str2num(highcandidates{i});
    highogfitvals = [highogfitvals,highvals(2)];
    highpsofitvals = [highpsofitvals,highvals(3)];
end

for i = 1:length(massgapcandidates)
    massgapvals = str2num(massgapcandidates{i});
    massgapogfitvals = [massgapogfitvals,massgapvals(2)];
    massgappsofitvals = [massgappsofitvals,massgapvals(3)];
end

for i = 1:length(highmassrealcandidates)
    highmassrealvals = str2num(highmassrealcandidates{i});
    highmassrealogfitvals = [highmassrealogfitvals,highmassrealvals(2)];
    highmassrealpsofitvals = [highmassrealpsofitvals,highmassrealvals(3)];
end
for i = 1:length(highmasslowcandidates)
    highmasslowvals = str2num(highmasslowcandidates{i});
    highmasslowogfitvals = [highmasslowogfitvals,highmasslowvals(2)];
    highmasslowpsofitvals = [highmasslowpsofitvals,highmasslowvals(3)];
end
for i = 1:length(highmasshighcandidates)
    highmasshighvals = str2num(highmasshighcandidates{i});
    highmasshighogfitvals = [highmasshighogfitvals,highmasshighvals(2)];
    highmasshighpsofitvals = [highmasshighpsofitvals,highmasshighvals(3)];
end

% ogfitvals = [realogfitvals,lowogfitvals,highogfitvals,massgapogfitvals,highmassrealogfitvals,highmasslowogfitvals,highmasshighogfitvals];
% psofitvals = [realpsofitvals,lowpsofitvals,highpsofitvals,massgappsofitvals,highmassrealpsofitvals,highmasslowpsofitvals,highmasshighpsofitvals];
ogfitvals = [realogfitvals,lowogfitvals,highogfitvals];
psofitvals = [realpsofitvals,lowpsofitvals,highpsofitvals];

percentage_reldifference_fitvals = (psofitvals - ogfitvals)*100./ogfitvals;

x = linspace(min(ogfitvals),max(ogfitvals));
y = zeros(1,length(x));
sz = 100;
figure;
% % plot(x,y);
hold on;
% scatter(ogfitvals,psofitvals - ogfitvals,'red','filled','MarkerEdgeColor','black');
scatter(ogfitvals,percentage_reldifference_fitvals,sz,'red','filled','MarkerEdgeColor','black');
hold off;
xline(9,'LineWidth',2,'Color','blue');
xlabel('$\hat{\rho}_{\rm true}$','Interpreter','latex');
ylabel('$[\hat{\rho}_{\rm PSO} - \hat{\rho}_{\rm true}]/\hat{\rho}_{\rm true} (\%)$','Interpreter','latex');
ax = gca;
ax.XAxis.FontSize = 60; ax.YAxis.FontSize = 60;
axis tight;

