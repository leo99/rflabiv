function [f, Fo, Ro, Rn] = readNoiseParamFile(fname)
% !FREQ   Fopt    GAMMA OPT       RN/Zo
% !GHZ     dB     MAG     ANG      -

fid = fopen(fname);
% Read actual data & close file
data = textscan(fid, '%f %f %f %f %f', 'headerlines', 3);
fclose(fid);

f = data{1} * 1e9;
Fo = 10 .^ (data{2} / 10);
Ro = data{3} .* exp(1i .* data{4} / 180 * pi);
Rn = data{5};

end