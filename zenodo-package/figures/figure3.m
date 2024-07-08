%Figure 3
phase = 0;
ta = 88.7;
datalen = 512;
snr = 15;
sampFreq = 4096;
masses = [-0.6,-1.5];
type = 3;
signal1 = gensignal(masses, phase, ta, datalen, snr, sampFreq, type);
masses = [-1.5,-0.6];
type = 3;
signal2 = gensignal(masses, phase, ta, datalen, snr, sampFreq, type);
%signal1 is signal vector without swapping the chirp times. signal2 is
%signal vector for swapping. Chirp times = [-0.6,-1.5], TOA = 88.7 and SNR
% = 15, both vectors are 512 sec long. 

timeVec = (0:512*4096-1)*(1/4096);
range1 = 88*4096:(89.5)*4096;
range2 = 88*4096:(91)*4096;
trange1 = timeVec(range1);
trange2 = timeVec(range2);
[s1,f1,t1] = spectrogram(signal1(range1), 256, [],[], 4096);
logS1 = log10(abs(s1));
[s2,f2,t2] = spectrogram(signal2(range2), 256, [],[], 4096);
logS2 = log10(abs(s2));

fmax = find(f1 < 1000,1, 'last');

figure;
% [ha, pos] = tight_subplot(2, 2, [0.16, 0.16]);
t = tiledlayout(2, 2, "TileSpacing", "compact");
% axes(ha(1)); 
nexttile;
plot(trange1,signal1(range1)./5);
ax = gca(); set(gca, 'XTickLabel', []); ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
ylabel('h(t)');

% axes(ha(2));
nexttile;
plot(trange2,signal2(range2)./5);
ax = gca(); ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);

% axes(ha(3));
nexttile;
imagesc(trange1,f1(1:fmax),logS1); axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
c = colorbar;
c.FontSize = 20;
ax = gca(); ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;

% axes(ha(4));
nexttile;
imagesc(trange2,f2(1:fmax),logS2); axis xy;  
xlabel('Time (s)');
c = colorbar;
c.FontSize = 20;
ax = gca(); ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
set(gca, 'YTickLabel', []);
% axis tight;


