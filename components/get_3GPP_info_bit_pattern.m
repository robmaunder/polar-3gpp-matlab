function info_bit_pattern = get_3GPP_info_bit_pattern(I, Q_N, rate_matching_pattern, mode)
% GET_3GPP_INFO_BIT_PATTERN Obtain the 3GPP information bit pattern, 
% according to Section 5.3.1.2 of 3GPP TS 38.212
%   info_bit_pattern = GET_3GPP_INFO_BIT_PATTERN(I, Q_N, rate_matching_pattern, mode)
%   obtains the information bit pattern info_bit_pattern.
%
%   I should be an integer scalar. It specifies the number of bits in the 
%   information, CRC and PC bit sequence. It should be no greater than N or E.
%
%   Q_N should be a row vector comprising N number of unique integers in the 
%   range 1 to N. Each successive element of Q_N provides the index of the
%   next most reliable input to the polar encoder kernal, where the first
%   element of Q_N gives the index of the least reliable bit and the last
%   element gives the index of the most reliable bit.
%
%   rate_matching_pattern should be a row vector comprising E number of
%   integers, each having a value in the range 1 to N. Each integer
%   identifies which one of the N outputs from the polar encoder kernal
%   provides the corresponding bit in the encoded bit sequence e.
%
%   mode should have the value 'repetition', 'puncturing' or 'shortening'.
%   This specifies how the rate matching has been achieved. 'repetition'
%   indicates that some outputs of the polar encoder kernal are repeated in 
%   the encoded bit sequence two or more times. 'puncturing' and 
%   'shortening' indicate that some outputs of the polar encoder kernal 
%   have been excluded from the encoded bit sequence. In the case of
%   'puncturing' these excluded bits could have values of 0 or 1. In the
%   case of 'shortening' these excluded bits are guaranteed to have values
%   of 0.
%
%   info_bit_pattern will be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in info_bit_pattern having the value true will be I. These elements 
%   having the value true identify the positions of the information and 
%   CRC bits within the input to the polar encoder kernal. The
%   information bit arrangement can be achieved according to
%   u(info_bit_pattern) = a.
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

N = length(Q_N);
n = log2(N);
E = length(rate_matching_pattern);

if n ~= round(n)
    error('N should be a power of 2');
end
if I > N
    error('polar_3gpp_matlab:UnsupportedBlockLength','I should be no greater than N.');
end
if I > E
    error('polar_3gpp_matlab:UnsupportedBlockLength','I should be no greater than E.');
end
if max(rate_matching_pattern) > N
    error('rate_matching_pattern is not compatible with N');
end
if strcmp(mode,'repetition') 
    if E < N
        error('mode is not compatible with E');
    end
elseif strcmp(mode,'puncturing')
    if E >= N
        error('mode is not compatible with E');
    end
elseif strcmp(mode,'shortening')
     if E >= N
        error('mode is not compatible with E');
    end
else
    error('Unsupported mode');
end

%% This is how the rate matching is described in TS 38.212
% P = [0 1 2 4 3 5 6 7 8 16 9 17 10 18 11 19 12 20 13 21 14 22 15 23 24 25 26 28 27 29 30 31];
% 
% J = zeros(1,N);
% for n=0:N-1
%     i=floor(32*n/N);
%     J(n+1) = P(i+1)*(N/32)+mod(n,N/32);
% end
% 
% Q_Ftmp_N = [];
% if E < N
%     if I/E <= 7/16 % puncturing
%         for n=0:N-E-1
%             Q_Ftmp_N = [Q_Ftmp_N,J(n+1)];
%         end
%         if E >= 3*N/4
%             Q_Ftmp_N = [Q_Ftmp_N,0:ceil(3*N/4-E/2)-1];
%         else
%             Q_Ftmp_N = [Q_Ftmp_N,0:ceil(9*N/16-E/4)-1];
%         end
%     else % shortening
%         for n=E:N-1
%             Q_Ftmp_N = [Q_Ftmp_N,J(n+1)];
%         end
%     end
% end

%% This is an equivalent but more flexible version, which also works with any rate matching pattern
Q_Ftmp_N = setdiff(1:N,rate_matching_pattern)-1;
if strcmp(mode,'puncturing')
    if E >= 3*N/4
        Q_Ftmp_N = [Q_Ftmp_N,0:ceil(3*N/4-E/2)-1];
    else
        Q_Ftmp_N = [Q_Ftmp_N,0:ceil(9*N/16-E/4)-1];
    end
end

%% 
Q_Itmp_N = setdiff(Q_N-1,Q_Ftmp_N,'stable'); % -1 because TS 38.212 assumes that indices start at 0, not 1 like in Matlab

if length(Q_Itmp_N) < I
    error('polar_3gpp_matlab:UnsupportedBlockLength','Too many pre-frozen bits.');
end    

Q_I_N=Q_Itmp_N(end-I+1:end);
%Q_F_N= setdiff(Q_N-1,Q_I_N,'stable'); % -1 because TS 38.212 assumes that indices start at 0, not 1 like in Matlab

info_bit_pattern = false(1,N);
info_bit_pattern(Q_I_N+1) = true;