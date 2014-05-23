function [ out ] = savg( in, N )

if size(in, 2) > size(in, 1)
	in = in.';
end

temp = [ones(N, 1) * in(1); in; ones(N, 1) * in(end)];
out = zeros(size(in));

for i = -N : 1 : N
	out = out + temp(N + i + 1 : size(in, 1) + N + i);
end

out = out/numel(-N:1:N);


end

