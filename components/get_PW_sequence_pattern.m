function Q_N = get_PW_sequence_pattern(N)
% GET_PW_SEQUENCE_PATTERN Get the Partial Weight (PW) sequence.
%   Q_N = GET_PW_SEQUENCE_PATTERN(N) obtains the PW sequence, as described 
%   in Section 2.2 of R1-167209...
%   http://www.3gpp.org/ftp/TSG_RAN/WG1_RL1/TSGR1_86/Docs/R1-167209.zip
%
%   N should be an integer scalar, which should be a power of 2. It 
%   specifies the number of bits in the input and output of the polar 
%   encoder kernal. 
%
%   Q_N will be a row vector comprising N number of unique integers in the 
%   range 1 to N. Each successive element of Q_N provides the index of the
%   next most reliable input to the polar encoder kernal, where the first
%   element of Q_N gives the index of the least reliable bit and the last
%   element gives the index of the most reliable bit.
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

n = log2(N);
if round(n) ~=n
    error('N must be a power of 2');
end

W = zeros(1,N);
for i = 0:N-1
    B = fliplr(dec2bin(i,n)-'0');
    W(i+1) = sum(B.*2.^((0:n-1)/4)); % +1 because indexing starts at 1 in Matlab
end
    
[~,Q_N] = sort(W);

end