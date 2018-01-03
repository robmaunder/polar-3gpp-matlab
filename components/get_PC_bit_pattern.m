function PC_bit_pattern = get_PC_bit_pattern(info_bit_pattern, Q_N, n_PC, n_PC_wm)
% GET_PC_BIT_PATTERN Obtain the Parity Check (PC) bit pattern, 
% according to Section 5.3.1.2 of 3GPP TS 38.212
%   PC_bit_pattern = get_PC_bit_pattern(info_bit_pattern, Q_N, n_PC, n_PC_wm)
%   obtains the PC bit pattern. 
%
%   info_bit_pattern should be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in info_bit_pattern having the value true should be I, where 
%   I = A+P+n_PC. These elements having the value true identify the 
%   positions of the information, CRC and PC bits within the input to the 
%   polar encoder kernal.
%
%   Q_N should be a row vector comprising N number of unique integers in the 
%   range 1 to N. Each successive element of Q_N provides the index of the
%   next most reliable input to the polar encoder kernal, where the first
%   element of Q_N gives the index of the least reliable bit and the last
%   element gives the index of the most reliable bit.
%
%   n_PC should be an integer scalar. It specifies the number of PC bits to
%   use, where n_PC should be no greater than I.
%
%   n_PC_wm should be an integer scalar. It specifies the number of PC bits
%   that occupy some of the most reliable positions at the input to the
%   polar encoder kernal. The remaining n_PC-n_PC_wm PC bits occupy some of 
%   the least reliable positions at the input to the polar encoder kernal.
%   n_PC_wm should be no greater than n_PC.
%
%   PC_bit_pattern will be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in PC_bit_pattern having the value true will be n_PC. 
%   These elements having the value true identify the positions of the 
%   PC bits within the input to the polar encoder kernal.

%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.


I = sum(info_bit_pattern);
N = length(info_bit_pattern);

if log2(N) ~= round(log2(N))
    error('N should be a power of 2.');
end
if length(Q_N) ~= N
    error('info_bit_pattern and Q_N have different lengths.');
end
if n_PC > I
    error('n_PC should be no greater than I.');
end
if n_PC_wm > n_PC
    error('n_PC_wm should be no greater than n_PC.');
end

Q_I = 1:N;
Q_N_I = intersect(Q_N, Q_I(info_bit_pattern), 'stable');

G_N = get_G_N(N);
w_g = sum(G_N,2);

Q_tilde_N_I = Q_N_I(n_PC+1:end); % This is what it says in TS 38.212
%Q_tilde_N_I = Q_N_I(n_PC-n_PC_wm+1:end); % I think that this would be slightly more elegant

Q_tilde_N_I_flip = fliplr(Q_tilde_N_I);
[w_g_sorted, indices] = sort(w_g(Q_tilde_N_I_flip));

Q_N_PC = [Q_N_I(1:n_PC-n_PC_wm), Q_tilde_N_I_flip(indices(1:n_PC_wm))];

PC_bit_pattern = false(1,N);
PC_bit_pattern(Q_N_PC) = true;

end
 