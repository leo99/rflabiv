% Final
clear;
clc;
close all;

% Save figures
fig = struct( ...
	'save',		true, ...
	'wid',		30, ...
	'hei',		20, ...
	'col',		'rgb',...
	'res',		300, ...
	'lax',		true, ...
	'fmo',		'fixed', ...
	'fsi',		12, ...
	'lmo',		'scaled', ...
	'fmt',		'pdf' ...
	);
interpFunc = 'pchip'; 	% Interpolation function
Navg = 3;

% [f, S11, Z0] = readSparamFile('../meas/0521/s11.s1p');
% [f, S21, Z0] = readSparamFile('../meas/0521/s21.s1p');
% [f, S12, Z0] = readSparamFile('../meas/0521/s12.s1p');
% [f, S22, Z0] = readSparamFile('../meas/0521/s22.s1p');

[f, S, Z0] = readSparamFile('../meas/0522/sparam.s2p');
S11 = savg(S(:, 1), Navg);
S21 = savg(S(:, 2), Navg);
S12 = savg(S(:, 3), Navg);
S22 = savg(S(:, 4), Navg);

h = [0 0 0 0];
figure(1);
clf;

hold on;
%h(3) = plot(f / 1e9, 20*log10(abs(S12)), 'linewidth', 2, 'color', 'g');
h(1) = plot(f / 1e9, 20*log10(abs(S11)), 'linewidth', 2, 'color', 'r');
h(4) = plot(f / 1e9, 20*log10(abs(S22)), 'linewidth', 2, 'color', 'b');
h(2) = plot(f / 1e9, 20*log10(abs(S21)), 'linewidth', 2, 'color', 'k');
hold off;

le = {'$$|S_{11}|$$', '$$|S_{21}|$$', '$$|S_{12}|$$', '$$|S_{22}|$$'};
i = [1 2 4];

grid on;
xlabel('Frequency $$f$$ [GHz]', 'interpreter', 'latex');
ylabel('Magnitude $$|S_{ij}|$$ [dB]', 'interpreter', 'latex');
title('Measured $$S$$-parameters of the amplifier', 'interpreter', 'latex');
legend(h(i), le(i), ...
	'location', 'northeast', ...
	'interpreter', 'latex', ...
	'orientation', 'horizontal');
ylim([-50 20]);

if fig.save
	orient landscape;
	exportfig(gcf, ['../meas/0522/sparam.', fig.fmt], ...
		'width',		fig.wid, ...
		'height',		fig.hei, ...
		'color',		fig.col, ...
		'resolution',	fig.res, ...
		'LockAxes',		fig.lax, ...
		'FontMode',		fig.fmo, ...
		'FontSize',		fig.fsi, ...
		'LineMode',		fig.lmo, ...
		'Format',		fig.fmt	);
end

f = fopen('../meas/0522/iip3.txt', 'r');
data = textscan(f, '%f %f %f');
fclose(f);

Gdiv = -3.3;
Gatt = mean([-3-(-27), -20-(-43.9)]);

data{1} = data{1} + Gdiv;
data{2} = data{2} + Gatt;
data{3} = data{3} + Gatt;

%P1 = polyfit(data{1}, data{2}, 1); %
P1 = [1 18.5];
%P3 = polyfit(data{1}, data{3}, 1); %
P3 = [3 12];
x = [-25 10];

TOI = solve(sprintf('%.10g * x + %.10g == %.10g * x + %.10g', ...
	P1(1), P1(2), P3(1), P3(2)), 'x');

h = [0 0 0 0 0];
i = [1 2 5];
le = {'Fundamental', 'Third-order', 'Fundamental', 'Third-order', 'TOI'};

figure(2);
clf;
hold on;
h(1) = plot(data{1}, data{2}, 'ro', 'linewidth', 2);
h(2) = plot(data{1}, data{3}, 'bo', 'linewidth', 2);
h(3) = plot(x, P1(1) * x + P1(2), 'r-', 'linewidth', 2);
h(4) = plot(x, P3(1) * x + P3(2), 'b-', 'linewidth', 2);
h(5) = plot(TOI, P1(1) * TOI + P1(2), 'ko', 'linewidth', 2, 'markerfacecolor', 'k');
hold off;
grid on;

xlabel('Input Power $$P_\mathrm{in}$$ [dBm]', 'interpreter', 'latex');
ylabel('Output Power $$P_\mathrm{out}$$ [dBm]', 'interpreter', 'latex');
title('$$\mathit{IIP}_3$$ measurements', 'interpreter', 'latex');
legend(h(i), le(i), ...
	'location', 'southeast', ...
	'interpreter', 'latex', ...
	'orientation', 'vertical');

if fig.save
	orient landscape;
	exportfig(gcf, ['../meas/0522/toi.', fig.fmt], ...
		'width',		fig.wid, ...
		'height',		fig.hei, ...
		'color',		fig.col, ...
		'resolution',	fig.res, ...
		'LockAxes',		fig.lax, ...
		'FontMode',		fig.fmo, ...
		'FontSize',		fig.fsi, ...
		'LineMode',		fig.lmo, ...
		'Format',		fig.fmt	);
end

T0 = 290;
ENR = 24.8;
Tc = T0;
Th = 290 * log10(10^(ENR/10)+1);

Y11 = 10 ^ ((-80.5 - (-96)) / 10);
Y12 = 10 ^ ((-96.5 - (-96.8)) / 10);
Y13 = 10 ^ ((-90.3 - (-100.8)) / 10);

Y21 = 10 ^ ((-104.0 - (-109.6)) / 10);
Y22 = 10 ^ ((-82.5 - (-100.0)) / 10);
Y23 = 10 ^ ((-95.2 - (-108.1)) / 10);

Te11 = (Th + Y11*Tc) / (Y11-1);
Te12 = (Th + Y12*Tc) / (Y12-1);
Te13 = (Th + Y13*Tc) / (Y13-1);
Te21 = (Th + Y21*Tc) / (Y21-1);
Te22 = (Th + Y22*Tc) / (Y22-1);
Te23 = (Th + Y23*Tc) / (Y23-1);

F11 = 10*log10(1 + Te11 / T0) %     3.2242
F12 = 10*log10(1 + Te12 / T0) %    17.0482
F13 = 10*log10(1 + Te13 / T0) %     3.6934
F21 = 10*log10(1 + Te21 / T0) %     5.2158
F22 = 10*log10(1 + Te22 / T0) %     3.1451
F23 = 10*log10(1 + Te23 / T0) %     3.4009

icp = csvread('../meas/0522/icp.csv');
P = savg(icp(:, 1), Navg) / 1e3;
G = savg(icp(:, 2), Navg);

idx = round(mean([find(G >= mean(G(G > 14)) - 1, 1, 'last'), ...
	find(G <= mean(G(G > 14)) - 1, 1, 'first')]));

h = zeros(1, 2);
i = [1 2];
le = {'Gain', 'ICP'};

figure(3);
clf;
hold on;
h(1) = plot(P, G, 'r-', 'linewidth', 2);
h(2) = plot(P(idx), G(idx), 'ko', 'linewidth', 2, 'markerfacecolor', 'k');
hold off;
grid on;

xlabel('Input Power $$P_\mathrm{in}$$ [dBm]', 'interpreter', 'latex');
ylabel('Gain $$G$$ [dBm]', 'interpreter', 'latex');
title('$$\mathit{ICP}$$ measurements', 'interpreter', 'latex');
legend(h(i), le(i), ...
	'location', 'southwest', ...
	'interpreter', 'latex', ...
	'orientation', 'vertical');
xlim([-25 0]);

if fig.save
	orient landscape;
	exportfig(gcf, ['../meas/0522/icp.', fig.fmt], ...
		'width',		fig.wid, ...
		'height',		fig.hei, ...
		'color',		fig.col, ...
		'resolution',	fig.res, ...
		'LockAxes',		fig.lax, ...
		'FontMode',		fig.fmo, ...
		'FontSize',		fig.fsi, ...
		'LineMode',		fig.lmo, ...
		'Format',		fig.fmt	);
end