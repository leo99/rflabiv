% Pre-study
clear;
clc;
%close all;

% Save figures
figs = struct( ...
	'global',	true, ...
	'mu',		true, ...
	'SLstab',	true, ...
	'Sstab',	false, ...
	'Lstab',	false, ...
	'UFM',		true, ...
	'DeltaG',	true, ...
	'MSG',		false, ...
	'GFcirc',	true, ...
	'GPcirc',	false ...
	);

dParam = struct( ...
	'f',		2.5e9, ...
	'Vds',		2, ...
	'Ids',		30e-3, ...
	'Gmin',		13, ...
	'Fmax',		0.8, ...
	'RLmin',	15, ...
	'h',		0.787e-3, ...
	'er',		2.3, ...
	'stab',		true ...
	);

orig = struct( ...
	'S',		[], ...
	'ABCD',		[], ...
	'Delta',	[], ...
	'K',		[], ...
	'mu',		[], ...
	'CL',		[], ...
	'RL',		[], ...
	'CS',		[], ...
	'RS',		[] ...
	);

stab = struct( ...
	'S',		[], ...
	'ABCD',		[], ...
	'Delta',	[], ...
	'K',		[], ...
	'mu',		[], ...
	'CL',		[], ...
	'RL',		[], ...
	'CS',		[], ...
	'RS',		[] ...
	);

interpFunc = 'pchip'; 	% Interpolation function

% S: S11 S21 S12 S22
[f1, S, Z0] = readSparamFile('2V_30mA_S.s2p');
[f2, Fm, Ro, Rn] = readNoiseParamFile('2V_30mA_F.s2p');
orig.S = S;

% Transmission parameters
ABCD = repmat(1 ./ (2 * S(:, 2)), 1, 4) .* [ ...
	(1 + S(:, 1)) .* (1 - S(:, 4)) + S(:, 2) .* S(:, 3), ...
	Z0 * ( (1 + S(:, 1)) .* (1 + S(:, 4)) - S(:, 2) .* S(:, 3) ), ...
	1/Z0 * ( (1 - S(:, 1)) .* (1 - S(:, 4)) - S(:, 2) .* S(:, 3) ), ...
	(1 - S(:, 1)) .* (1 + S(:, 4)) + S(:, 2) .* S(:, 3) ...
	];
orig.ABCD = ABCD;

% Stability
Delta = S(:, 1) .* S(:, 4) - S(:, 2) .* S(:, 3);
K = (1 - abs(S(:, 1)).^2 - abs(S(:, 4)).^2 + abs(Delta).^2) ./ ( 2 * abs(S(:, 2) .* S(:, 3)) );
mu = (1 - abs(S(:, 1)).^2) ./ ( abs(S(:, 4) - Delta .* conj(S(:, 1))) + abs(S(:, 2) .* S(:, 3)) );
orig.Delta = Delta;
orig.K = K;
orig.mu = mu;

Rollet = ((abs(Delta) < 1) & (K > 1));
UST = (mu > 1);
Fstab = f1(Rollet == 1);
Sstr = {'Potentially unstable', 'Unconditionally stable'};

figure(1);
clf;
plot(f1 / 1e9, mu, 'color', 'b', 'linewidth', 2);
grid on;
ylabel('Stability $$\mu$$', 'interpreter', 'latex');
xlabel('Frequency $$f$$ [GHz]', 'interpreter', 'latex');
title(sprintf('At $$f = %.2f$$ GHz: $$\\mu = %.2f$$, $$K = %.2f$$ and $$|\\Delta| = %.2f$$ $$\\Rightarrow$$ %s', ...
	dParam.f/1e9 , mu(f1 == dParam.f), K(f1 == dParam.f), abs(Delta(f1 == dParam.f)), Sstr{1 + UST(f1 == dParam.f)}), 'interpreter', 'latex');

% Save figure
if figs.global && figs.mu
exportfig(gcf, 'mu.eps', ...
	'width',		12, ...
	'height',		6, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end

% Stability circles
CL = conj( S(:, 4) - Delta .* conj(S(:, 1)) ) ./ ( abs(S(:, 4)).^2 - abs(Delta).^2 );
RL = abs( S(:, 2) .* S(:, 3) ./ ( abs(S(:, 4)).^2 - abs(Delta).^2 ) );
CS = conj( S(:, 1) - Delta .* conj(S(:, 4)) ) ./ ( abs(S(:, 1)).^2 - abs(Delta).^2 );
RS = abs( S(:, 2) .* S(:, 3) ./ ( abs(S(:, 1)).^2 - abs(Delta).^2 ) );

orig.CL = CL; 
orig.RL = RL; 
orig.CS = CS; 
orig.RS = RS;

figure(2);
clf;
initSmithChart(gca);
hold on;
h = [0 0];

C = CL(f1 == dParam.f);
R = RL(f1 == dParam.f);
p = angle(C) + linspace(0, 2 * pi, 4001);
p = p(1 : end-1);
z = C + R * exp(1i * p);
z = z(abs(z) <= 1.02);
if ~isempty(z)
	h(1) = plot(gca, real(z), imag(z), ...
		'color',		'b', ...
		'linewidth',	2);
end

C = CS(f1 == dParam.f);
R = RS(f1 == dParam.f);
p = angle(C) + linspace(0, 2 * pi, 4001);
p = p(1 : end-1);
z = C + R * exp(1i * p);
z = z(abs(z) <= 1.02);
if ~isempty(z)
	h(2) = plot(gca, real(z), imag(z), ...
		'color',		'r', ...
		'linewidth',	2);
end
hold off;

text(0.35, 0.75, 'Stable', ...
	'HorizontalAlignment', 'center', ...
	'color',	'r', ...
	'backgroundColor',	'white', ...
	'interpreter', 'latex');
text(-0.15, 0.85, 'Unstable', ...
	'HorizontalAlignment', 'center', ...
	'color',	'r', ...
	'backgroundColor',	'white', ...
	'interpreter', 'latex');
text(-0.3, 0.2, 'Stable', ...
	'HorizontalAlignment', 'center', ...
	'color',	'b', ...
	'backgroundColor',	'white', ...
	'interpreter', 'latex');
text(-0.65, 0.4, 'Unstable', ...
	'HorizontalAlignment', 'center', ...
	'color',	'b', ...
	'backgroundColor',	'white', ...
	'interpreter', 'latex');

legend(gca, h, ...
	{sprintf('Output: $$\\rho_\\mathrm{L}(|\\rho_\\mathrm{in}| = 1$$) ($$|S_{22}| = %.2f$$)', abs(S(f1 == dParam.f, 4))), ...
	sprintf('Input: $$\\rho_\\mathrm{S}(|\\rho_\\mathrm{in}| = 1$$) ($$|S_{11}| = %.2f$$)', abs(S(f1 == dParam.f, 1)))}, ...
	'interpreter',	'latex', ...
	'location',		'south');

xlim([-1 1]);
ylim([-1 1]);

if figs.global && figs.SLstab
exportfig(gcf, 'SLstab.eps', ...
	'width',		12, ...
	'height',		12, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end


figure(3);
clf;
initSmithChart(gca);
hold on;

for i = 1 : length(f1)
	lineStyle = repmat('-', 1, 2 - mod(i, 2));
	C = CL(i);
	R = RL(i);
	p = angle(C) + linspace(0, 2 * pi, 4001);
	p = p(1 : end-1);
	z = C + R * exp(1i * p);
	z = z(abs(z) <= 1.02);
	if ~isempty(z)
	plot(gca, real(z), imag(z), ...
		'color',		'b', ...
		'lineStyle',	lineStyle, ...
		'linewidth',	2);
	end
end
hold off;
xlim([-1 1]);
ylim([-1 1]);
text(0, -0.2, 'Output stability circles: $$|S_{22}| < 1$$', ...
	'HorizontalAlignment', 'center', ...
	'backgroundColor',	'white', ...
	'interpreter', 'latex');

if figs.global && figs.Lstab
exportfig(gcf, 'Lstab.eps', ...
	'width',		12, ...
	'height',		12, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end

figure(4);
clf;
initSmithChart(gca);
hold on;

for i = 1 : length(f1)
	lineStyle = repmat('-', 1, 2 - mod(i, 2));
	C = CS(i);
	R = RS(i);
	p = angle(C) + linspace(0, 2 * pi, 4001);
	p = p(1 : end-1);
	z = C + R * exp(1i * p);
	z = z(abs(z) <= 1.02);
	plot(gca, real(z), imag(z), ...
		'color',		'r', ...
		'lineStyle',	lineStyle, ...
		'linewidth',	2);
end
hold off;
xlim([-1 1]);
ylim([-1 1]);
text(0, -0.2, 'Input stability circles: $$|S_{11}| < 1$$', ...
	'HorizontalAlignment', 'center', ...
	'backgroundColor',	'white', ...
	'interpreter', 'latex');

if figs.global && figs.Sstab
exportfig(gcf, 'Sstab.eps', ...
	'width',		18, ...
	'height',		18, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end

% Forced stabilization

% Stabilization resistors
% SerRS >= 15.3, ParRS <= 80.9, SerRL >= 37.8, ParRL <= 6.27
% To specifications: SerRL = [69.8, 109]
SerRS =	[	15.3,	0,		0,		0		];
ParRS = [	Inf,	80.9,	Inf,	Inf		];
SerRL = [	0,		0,		100,	0		];
ParRL = [	Inf,	Inf,	Inf,	6.27	];

for j = 3 %1 : 4
	ABCDs = [1,		SerRS(j);	1/ParRS(j),		1];
	ABCDl = [1,		SerRL(j);	1/ParRL(j),		1];

	ABCD2 = zeros(size(ABCD));
	%KK = [0 0; 0 0];

	for i = 1 : length(f1)
		KK = ABCDs * [ABCD(i, 1), ABCD(i, 2); ABCD(i, 3), ABCD(i, 4)] * ABCDl;
		ABCD2(i,:) = [KK(1, 1), KK(1, 2), KK(2, 1), KK(2, 2)];
	end

	% Revert back to S-param
	SS = repmat(1 ./ ( ABCD2(:, 1) + ABCD2(:, 2) / Z0 + ABCD2(:, 3) * Z0 + ABCD2(:, 4) ), 1, 4) .* [ ...
		ABCD2(:, 1) + ABCD2(:, 2) / Z0 - ABCD2(:, 3) * Z0 - ABCD2(:, 4), ...
		repmat(2, size(ABCD, 1), 1), ...
		2 * (ABCD2(:, 1) .* ABCD2(:, 4) - ABCD2(:, 2) .* ABCD2(:, 3)), ...
		-ABCD2(:, 1) + ABCD2(:, 2) / Z0 - ABCD2(:, 3) * Z0 + ABCD2(:, 4) ...
		];

	Delta2 = SS(:, 1) .* SS(:, 4) - SS(:, 2) .* SS(:, 3);
	K2 = (1 - abs(SS(:, 1)).^2 - abs(SS(:, 4)).^2 + abs(Delta2).^2) ./ ( 2 * abs(SS(:, 2) .* SS(:, 3)) );
	mu2 = (1 - abs(SS(:, 1)).^2) ./ ( abs(SS(:, 4) - Delta2 .* conj(SS(:, 1))) + abs(SS(:, 2) .* SS(:, 3)) );

	fprintf('j = %d: mu = %.6f, G0 = %.2f dB\n', j, mu2(f1 == dParam.f), 20*log10(SS(f1 == dParam.f, 2)));
end

% Use stabilized or original design
if dParam.stab
	S = SS;
	Delta = Delta2;
	K = K2;
	mu = mu2;
	CL = conj( S(:, 4) - Delta .* conj(S(:, 1)) ) ./ ( abs(S(:, 4)).^2 - abs(Delta).^2 );
	RL = abs( S(:, 2) .* S(:, 3) ./ ( abs(S(:, 4)).^2 - abs(Delta).^2 ) );
	CS = conj( S(:, 1) - Delta .* conj(S(:, 4)) ) ./ ( abs(S(:, 1)).^2 - abs(Delta).^2 );
	RS = abs( S(:, 2) .* S(:, 3) ./ ( abs(S(:, 1)).^2 - abs(Delta).^2 ) );
	
	stab.S = SS;
	stab.Delta = Delta2;
	stab.K = K2;
	stab.mu = mu2;
	stab.CL = CL;
	stab.RL = RL;
	stab.CS = CS;
	stab.RS = CS;
end

% Unilateral figure of merit
U = abs( S(:, 1) .* S(:, 2) .* S(:, 3) .* S(:, 4) ) ./ ...
	( (1 - abs(S(:, 1)).^2) .* (1 - abs(S(:, 4)).^2) );
DiffLow = 1 ./ (1 + U).^2;
DiffHigh = 1 ./ (1 - U).^2;

figure(5);
clf;
plot(f1 / 1e9, U, 	...
	'color',		'b', ...
	'linewidth',	2);
grid on;
ylabel('Unilateral Figure of Merit $$U$$', 'interpreter', 'latex');
xlabel('Frequency $$f$$ [GHz]', 'interpreter', 'latex');
title(sprintf('At $$f = %.2f$$ GHz: $$U = %.2f$$', ...
	dParam.f/1e9 , U(f1 == dParam.f)), 'interpreter', 'latex');

if figs.global && figs.UFM
exportfig(gcf, 'UFM.eps', ...
	'width',		12, ...
	'height',		6, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end


figure(6);
clf;
hold on;
plot(f1 / 1e9, 10*log10(DiffLow), 	...
	'color',		'r', ...
	'linewidth',	2);
plot(f1 / 1e9, 10*log10(DiffHigh), 	...
	'color',		'b', ...
	'linewidth',	2);
hold off;
grid on;
ylabel('Error caused by unilater assumption $$\Delta G$$ [dB]', 'interpreter', 'latex');
xlabel('Frequency $$f$$ [GHz]', 'interpreter', 'latex');
title(sprintf('At $$f = %.2f$$ GHz: $$%.2f < G_\\mathrm{T}/G_\\mathrm{TU} < %.2f$$ [dB]', ...
	dParam.f/1e9 , 10*log10(DiffLow(f1 == dParam.f)), 10*log10(DiffHigh(f1 == dParam.f))), 'interpreter', 'latex');
xlim([0 18]);
%set(gca, 'XTick', 0:2:18);


legend(gca, ...
	{'Lower limit', 'Higher limit'}, ...
	'interpreter',	'latex', ...
	'location',		'northeast');

if figs.global && figs.DeltaG
exportfig(gcf, 'DeltaG.eps', ...
	'width',		18, ...
	'height',		10, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end

% MSG
MSG = abs(S(:, 2) ./ S(:, 3));

figure(7);
clf;
plot(f1 / 1e9, 10*log10(MSG), ...
	'color',		'b', ...
	'linewidth',	2);
grid on;
ylabel('Maximum Stable Gain $$\mathit{MSG}$$ [dB]', 'interpreter', 'latex');
xlabel('Frequency $$f$$ [GHz]', 'interpreter', 'latex');
title(sprintf('At $$f = %.2f$$ GHz: $$\\mathit{MSG} = %.1f$$ dB', ...
	dParam.f/1e9, 10*log10(MSG(f1 == dParam.f))), 'interpreter', 'latex');

if figs.global && figs.MSG
exportfig(gcf, 'MSG.eps', ...
	'width',		12, ...
	'height',		6, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end

% Constant Available Gain & Noise Figure Circles
figure(8);
clf;
initSmithChart(gca);
hold on;
h = [0 0 0];

C = CS(f1 == dParam.f);
R = RS(f1 == dParam.f);
p = angle(C) + linspace(0, 2 * pi, 4001);
p = p(1 : end-1);
z = C + R * exp(1i * p);
z = z(abs(z) <= 1.02);
if ~isempty(z)
	h(1) = plot(gca, real(z), imag(z), ...
		'color',		'k', ...
		'linewidth',	2);
end

%vGA = [12 13 14 15 16 17 18];
vGA = [10 11 12 13 13.5];
vFi = [0.331 0.4 0.6 0.8 1];

vGA = 10 .^ (vGA / 10);
vFi = 10 .^ (vFi / 10);

C1 = S(:, 1) - Delta .* conj(S(:, 4));

hh = zeros(size(vGA));
for i = 1 : length(vGA)
	GA = vGA(i);
	lineStyle = repmat('-', 1, 2 - mod(i, 2));
	%ga = ( 1 - abs(Rs).^2 ) ./ ( 1 - S(:, 4).^2 + abs(Rs).^2 .* ( abs(S(:, 1).^2 - abs(Delta).^2) - 2 * real(Rs .* C1) );
	ga = GA ./ abs(S(:, 2)).^2;

	Ca = ( ga .* conj(C1) ) ./ ( 1 + ga .* ( abs(S(:, 1)).^2 - abs(Delta).^2 ) );
	Ra = sqrt( 1 - 2 * K .* abs( S(:, 2) .* S(:, 3) ) .* ga + abs( S(:, 2) .* S(:, 3) ).^2 .* ga.^2 ) ./ ...
		abs( 1 + ga .* (abs(S(:, 1)).^2 - abs(Delta).^2 ) );
	
	C = Ca(f1 == dParam.f);
	R = abs(Ra(f1 == dParam.f));
	p = angle(C) + linspace(0, 2 * pi, 4001);
	p = p(1 : end-1);
	z = C + R * exp(1i * p);
	z = z(abs(z) <= 1.02);
	if ~isempty(z)
		hh(i) = plot(gca, real(z), imag(z), ...
			'color',		'b', ...
			'lineStyle',	lineStyle, ...
			'linewidth',	2);
	end
end
h(2) = hh(1);

hh = zeros(size(vFi));
for i = 1 : length(vFi)
	Fi = vFi(i);
	lineStyle = repmat('-', 1, 2 - mod(i, 2));
	Ni = (Fi - Fm) ./ (4 * Rn) .* abs(1 + Ro).^2;
	CFi = Ro ./ (1 + Ni);
	RFi = sqrt(Ni.^2 + Ni .* (1 - abs(Ro).^2)) ./ (1 + Ni);
	
	C = CFi(f1 == dParam.f);
	R = abs(RFi(f1 == dParam.f));
	p = angle(C) + linspace(0, 2 * pi, 4001);
	p = p(1 : end-1);
	z = C + R * exp(1i * p);
	z = z(abs(z) <= 1.02);
	if ~isempty(z)
		hh(i) = plot(gca, real(z), imag(z), ...
			'color',		'r', ...
			'lineStyle',	lineStyle, ...
			'linewidth',	2);
	end
end
hold off;
h(3) = hh(1);
if mu(f1 == dParam.f) > 1
	h = h(2:3);
end;

Fstr = sprintf(', %.2f', 10*log10(vFi));
Gstr = sprintf(', %.1f', 10*log10(vGA));

legend(gca, h, ...
	{sprintf('$$G_\\mathrm{A} = %s$$ dB', Gstr(3:end)), ...
	sprintf('$$F = %s$$ dB', Fstr(3:end)), ...
	sprintf('Input stability circle: $$|S_{11}| = %.2f$$', abs(S(f1 == dParam.f, 1)))}, ...
	'interpreter',	'latex', ...
	'location',		'south');

if figs.global && figs.GFcirc
exportfig(gcf, 'GFcirc.eps', ...
	'width',		12, ...
	'height',		12, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end


% Constant Operating Gain
figure(9);
clf;
initSmithChart(gca);
hold on;
h = [0 0];

C = CL(f1 == dParam.f);
R = RL(f1 == dParam.f);
p = angle(C) + linspace(0, 2 * pi, 4001);
p = p(1 : end-1);
z = C + R * exp(1i * p);
z = z(abs(z) <= 1.02);
if ~isempty(z)
	h(1) = plot(gca, real(z), imag(z), ...
		'color',		'k', ...
		'linewidth',	2);
end

%vGP = [12 13 14 15 16 17 18];
vGP = [10 11 12 13 13.5];
vGP = 10 .^ (vGP / 10);

C2 = S(:, 4) - Delta .* conj(S(:, 1));

hh = zeros(size(vGP));
for i = 1 : length(vGP)
	GP = vGP(i);
	lineStyle = repmat('-', 1, 2 - mod(i, 2));
	%ga = ( 1 - abs(Rs).^2 ) ./ ( 1 - S(:, 4).^2 + abs(Rs).^2 .* ( abs(S(:, 1).^2 - abs(Delta).^2) - 2 * real(Rs .* C1) );
	gp = GP ./ abs(S(:, 2)).^2;

	Cp = ( gp .* conj(C2) ) ./ ( 1 + gp .* ( abs(S(:, 4)).^2 - abs(Delta).^2 ) );
	Rp = sqrt( 1 - 2 * K .* abs( S(:, 2) .* S(:, 3) ) .* gp + abs( S(:, 2) .* S(:, 3) ).^2 .* gp.^2 ) ./ ...
		abs( 1 + gp .* (abs(S(:, 4)).^2 - abs(Delta).^2 ) );
	
	C = Cp(f1 == dParam.f);
	R = Rp(f1 == dParam.f);
	p = angle(C) + linspace(0, 2 * pi, 4001);
	p = p(1 : end-1);
	z = C + R * exp(1i * p);
	z = z(abs(z) <= 1.02);
	if ~isempty(z)
		hh(i) = plot(gca, real(z), imag(z), ...
			'color',		'b', ...
			'lineStyle',	lineStyle, ...
			'linewidth',	2);
	end
end
h(2) = hh(1);
if mu(f1 == dParam.f) > 1
	h = h(2);
end;

Gstr = sprintf(', %.1f', 10*log10(vGP));

legend(gca, h, ...
	{sprintf('$$G_\\mathrm{P} = %s$$ dB', Gstr(3:end)), ...
	sprintf('Output stability circle: $$|S_{22}| = %.2f$$', abs(S(f1 == dParam.f, 4)))}, ...
	'interpreter',	'latex', ...
	'location',		'south');

if figs.global && figs.GPcirc
exportfig(gcf, 'GFcirc.eps', ...
	'width',		12, ...
	'height',		12, ...
	'color',		'rgb', ...
	'resolution',	300, ...
	'LockAxes',		1, ...
	'FontMode',		'fixed', ...
	'FontSize',		9, ...
	'LineMode',		'scaled'	);
end

% Constant VSWR circles -> Needs Rin and Rout
% vRL = [10 15 20];
% Rab = 10^(-vRL/20);

% CVi = conj(Rin) .* (1 - Rab.^2) ./ (1 - abs(Rab .* Rin));
% RVi = Rab .* (1 - abs(Rin).^2) ./ (1 - abs(Rab .* Rin));
% CVo = conj(Rout) .* (1 - Rab.^2) ./ (1 - abs(Rab .* Rout));
% RVo = Rab .* (1 - abs(Rout).^2) ./ (1 - abs(Rab .* Rout));

% Gains & Noise figure for chosen Rs and Rl

% Is it stable inside or outside the stability circle
sInOut = ((abs(S(:, 1)) < 1) & abs(CS) < RS) | ((abs(S(:, 1)) > 1) & abs(CS) > RS);
lInOut = ((abs(S(:, 4)) < 1) & abs(CL) < RL) | ((abs(S(:, 4)) > 1) & abs(CL) > RL);
margin = 0.05;
YesNo = {'Yes', 'No'};

Fi = 10^(0.737/10);
Ni = (Fi - Fm) ./ (4 * Rn) .* abs(1 + Ro).^2;
CFi = Ro ./ (1 + Ni);
RFi = sqrt(Ni.^2 + Ni .* (1 - abs(Ro).^2)) ./ (1 + Ni);
C = CFi(f1 == dParam.f);
R = abs(RFi(f1 == dParam.f));

% GA = 10^(14.7/10);
% ga = GA ./ abs(S(:, 2)).^2;
% Ca = ( ga .* conj(C1) ) ./ ( 1 + ga .* ( abs(S(:, 1)).^2 - abs(Delta).^2 ) );
% Ra = sqrt( 1 - 2 * K .* abs( S(:, 2) .* S(:, 3) ) .* ga + abs( S(:, 2) .* S(:, 3) ).^2 .* ga.^2 ) ./ ...
% 	abs( 1 + ga .* (abs(S(:, 1)).^2 - abs(Delta).^2 ) );
% C = Ca(f1 == dParam.f);
% R = abs(Ra(f1 == dParam.f));

p = angle(C) + linspace(0, 2 * pi, 3601);
p = p(1 : end-1);
z = C + R * exp(1i * p);
z = z(abs(z) < 1);
i = find(f1 == dParam.f);

for Rs = z
	%Rs = 0.55 + 0.62 * 1i; %-0.3 + 0.8 * 1i; % Ro(f1 == dParam.f);
	Rout = S(:, 4) + (S(:, 2) .* S(:, 3) .* Rs) ./ (1 - S(:, 1) .* Rs);
	Rl = conj(Rout((f1 == dParam.f)));
	Rin = S(:, 1) + (S(:, 2) .* S(:, 3) .* Rl) ./ (1 - S(:, 4) .* Rl);

	Sstable = sqrt((real(CS) - real(Rs)).^2 + (imag(CS) - imag(Rs)).^2);
	Lstable = sqrt((real(CL) - real(Rl)).^2 + (imag(CL) - imag(Rl)).^2);

	Sstable = sInOut & (Sstable < (RS - margin)) | ~sInOut & (Sstable > (RS + margin));
	Lstable = lInOut & (Lstable < (RL - margin)) | ~lInOut & (Lstable > (RL + margin));

	Ga = abs(S(:, 2)).^2 .* (1 - abs(Rs).^2) ./ ( abs(1 - S(:, 1) .* Rs).^2 .* (1 - abs(Rout).^2) );
	Gt = abs(S(:, 2)).^2 .* (1 - abs(Rs).^2) .* (1 - abs(Rl).^2) ./ ...
		( abs(1 - Rs .* Rin).^2 .* abs(1 - S(:, 4) .* Rl).^2 );
	Gp = abs(S(:, 2)).^2 .* (1 - abs(Rl).^2) ./ ...
		( (1 - abs(Rin).^2) .* abs(1 - S(:, 4) .* Rl).^2 );
	GadB = 10*log10(Ga);
	GtdB = 10*log10(Gt);
	GpdB = 10*log10(Gp);

	F = Fm + 4 * Rn .* abs(Rs - Ro).^2 ./ ( (1 - abs(Rs).^2) .* abs(1 + Ro).^2 );
	FdB = 10*log10(F);

	RLi = 1 ./ (1 - Gt ./ Gp);
	RLo = 1 ./ (1 - Gt ./ Ga);
	RLi(abs(RLi) > 1e12) = Inf;
	RLo(abs(RLo) > 1e12) = Inf;
	RLidB = 10 * log10(RLi);
	RLodB = 10 * log10(RLo);
	
	if (RLidB(f1 == dParam.f) > dParam.RLmin && ...
			RLodB(f1 == dParam.f) > dParam.RLmin && ...
			GtdB(f1 == dParam.f) > dParam.Gmin && ...
			FdB(f1 == dParam.f) < dParam.Fmax) && ...
			abs(Rs) < 1 && Sstable(f1 == dParam.f) && ...
			abs(Rl) < 1 && Lstable(f1 == dParam.f) && ...
			MSG(f1 == dParam.f) > Gt(f1 == dParam.f)
		
		fprintf('Rs   = %+.3f %+.3fj (%.3f /_ %+.2f °)\n', real(Rs), imag(Rs), abs(Rs), angle(Rs) / pi * 180);
		fprintf('Rin  = %+.3f %+.3fj (%.3f /_ %+.2f °)\n', real(Rin(f1 == dParam.f)), imag(Rin(f1 == dParam.f)), abs(Rin(f1 == dParam.f)), angle(Rin(f1 == dParam.f)) / pi * 180);
		fprintf('Rl   = %+.3f %+.3fj (%.3f /_ %+.2f °)\n', real(Rl), imag(Rl), abs(Rl), angle(Rl) / pi * 180);
		fprintf('Rout = %+.3f %+.3fj (%.3f /_ %+.2f °)\n', real(Rout(f1 == dParam.f)), imag(Rout(f1 == dParam.f)), abs(Rout(f1 == dParam.f)), angle(Rout(f1 == dParam.f)) / pi * 180);
		
		%fprintf('Rs practical:   %s\n', YesNo{2 - (abs(Rs) < 1 && Sstable(f1 == dParam.f))});
		%fprintf('Rl practical:   %s\n', YesNo{2 - (abs(Rl) < 1 && Lstable(f1 == dParam.f))});
		%fprintf('Gain practical: %s\n', YesNo{2 - (MSG(f1 == dParam.f) > Gt(f1 == dParam.f))});
		fprintf('RLi  = %s dB (> 15)\n', num2str(RLidB(f1 == dParam.f), '%.2f'));
		fprintf('RLo  = %s dB (> 15)\n', num2str(RLodB(f1 == dParam.f), '%.2f'));
		fprintf('Ga   = %s dB\n', num2str(GadB(f1 == dParam.f), '%.2f'));
		fprintf('Gt   = %s dB (> 13)\n', num2str(GtdB(f1 == dParam.f), '%.2f'));
		fprintf('Gp   = %s dB\n', num2str(GpdB(f1 == dParam.f), '%.2f'));
		fprintf('F    = %s dB (< 0.8)\n', num2str(FdB(f1 == dParam.f), '%.3f'));
		fprintf('\n');
	end
	
% 	if (RLidB(f1 == dParam.f) > dParam.RLmin && ...
% 			RLodB(f1 == dParam.f) > dParam.RLmin && ...
% 			GtdB(f1 == dParam.f) > dParam.Gmin && ...
% 			FdB(f1 == dParam.f) < dParam.Fmax)
% 		fprintf('Specifications met!\n\n');
% 	else
% 		fprintf('Specifications NOT met!\n\n');
% 	end
end


% Stable: RL > 15 dB -> F > 1.7 dB or F < 0.8 dB -> RL < 9 dB
% Unstable: RL > 15 dB not possible or F < 0.8 -> RL < 1.1 dB

% Rs  = -0.208 +0.516j (0.556 /_ 112.00 °)
% Rl  = 0.356 +0.263j (0.443 /_ 36.40 °)
% RLi = 15.01 dB (> 15)
% RLo = Inf dB (> 15)
% Ga  = 14.74 dB
% Gt  = 14.74 dB (> 13)
% Gp  = 14.88 dB
% F   = 0.799 dB (< 0.8)
% SerRL = 69.8

% Rs  = -0.187 +0.443j (0.480 /_ 112.86 °)
% Rl  = 0.487 +0.152j (0.510 /_ 17.37 °)
% RLi = 15.01 dB (> 15)
% RLo = Inf dB (> 15)
% Ga  = 13.02 dB
% Gt  = 13.02 dB (> 13)
% Gp  = 13.16 dB
% F   = 0.722 dB (< 0.8)
% SerRL = 109

% Rs  = -0.191 +0.469j (0.506 /_ 112.14 °)
% Rl  = +0.440 +0.185j (0.478 /_ 22.81 °)
% RLi = 15.00 dB (> 15)
% RLo = Inf dB (> 15)
% Ga  = 13.61 dB
% Gt  = 13.61 dB (> 13)
% Gp  = 13.75 dB
% F   = 0.740 dB (< 0.8)
% SerRL = 93.3

% Rs  = -0.191 +0.457j (0.495 /_ 112.72 °)
% Rl  = +0.468 +0.166j (0.497 /_ 19.52 °)
% RLi = 15.25 dB (> 15)
% RLo = Inf dB (> 15)
% Ga  = 13.26 dB
% Gt  = 13.26 dB (> 13)
% Gp  = 13.39 dB
% F   = 0.735 dB (< 0.8)
% Ser-RL = 102.5 ohm

% Rs   = -0.191 +0.461j (0.499 /_ 112.57 °)
% Rin  = -0.307 -0.503j (0.589 /_ -121.45 °)
% Rl   = +0.461 +0.171j (0.492 /_ 20.36 °)
% Rout = +0.461 -0.171j (0.492 /_ -20.36 °)
% RLi  = 15.22 dB (> 15)
% RLo  = Inf dB (> 15)
% Ga   = 13.35 dB
% Gt   = 13.35 dB (> 13)
% Gp   = 13.49 dB
% F    = 0.737 dB (< 0.8)
% Par-RL = 100 ohm