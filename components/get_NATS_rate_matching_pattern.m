function [rate_matching_pattern, mode] = get_NATS_rate_matching_pattern(N,E)
% GET_NATS_RATE_MATCHING_PATTERN Get the Natural Shortening (NATS) 
% sequence, as specified in R1-1704318...
% http://www.3gpp.org/ftp/TSG_RAN/WG1_RL1/TSGR1_88b/Docs/R1-1704318.zip
%   [rate_matching_pattern, mode] = get_NATS_rate_matching_pattern(N,E)
%   obtains the rate matching sequence.
%
%   N should be an integer scalar, which should be a power of 2. It 
%   specifies the number of bits in the input and output of the polar 
%   encoder kernal. 
%
%   E should be an integer scalar, which should no greater than N. It 
%   specifies the number of bits in the encoded bit sequence. 
%
%   rate_matching_pattern will be a row vector comprising E number of
%   integers, each having a value in the range 1 to N. Each integer
%   identifies which one of the N outputs from the polar encoder kernal
%   provides the corresponding bit in the encoded bit sequence. Rate
%   matching can be achieved according to e = d(rate_matching_pattern).
%
%   mode will have the value 'repetition', 'puncturing' or 'shortening'.
%   This specifies how the rate matching has been achieved. 'repetition'
%   indicates that some outputs of the polar encoder kernal are repeated in 
%   the encoded bit sequence two or more times. 'puncturing' and 
%   'shortening' indicate that some outputs of the polar encoder kernal 
%   have been excluded from the encoded bit sequence. In the case of
%   'puncturing' these excluded bits could have values of 0 or 1. In the
%   case of 'shortening' these excluded bits are guaranteed to have values
%   of 0.
%
%   See also GET_BIVP_RATE_MATCHING_PATTERN, GET_BIVS_RATE_MATCHING_PATTERN 
%   and GET_NATP_RATE_MATCHING_PATTERN
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

if round(log2(N))~=log2(N)
    error('N should be a power of 2');
end
if E>N
    error('polar_3gpp_matlab:UnsupportedBlockLength','E should be no greater than N');
end

rate_matching_pattern = 1:E;
mode = 'shortening';