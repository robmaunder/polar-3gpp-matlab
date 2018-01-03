function channel_interleaver_pattern = get_3GPP_channel_interleaver_pattern(E)
% GET_3GPP_CHANNEL_INTERLEAVER_PATTERN Obtain the 3GPP channel 
% interleaver pattern, according to Section 5.4.1.3 of 3GPP TS 38.212 
%   channel_interleaver_pattern = GET_3GPP_CHANNEL_INTERLEAVER_PATTERN(E) obtains the channel interleaver
%   pattern.
%
%   E should be an integer scalar. It specifies the number of bits in the 
%   encoded bit sequence. 
%
%   channel_interleaver_pattern will be an integer row vector, compring E 
%   unique elements in the range 1 to E. Each integer identifies which one 
%   of the E encoded bits provides the corresponding interleaved bit. 
%   Interleaving can be implemented according to 
%   f = e(channel_interleaver_pattern), while deinterleaving is implemented 
%   accoring to e_tilde(channel_interleaver_pattern) = f_tilde.
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

T = 0;
while T*(T+1)/2 < E
    T = T+1;
end

v = zeros(T,T);

k=0;
for i = 0:T-1
    for j=0:T-1-i
        if k<E
            v(i+1,j+1) = k+1;
        else
            v(i+1,j+1) = NaN;
        end
        k = k+1;
    end
end

channel_interleaver_pattern = zeros(1,E);
k=0;
for j = 0:T-1
    for i=0:T-1-j
        if ~isnan(v(i+1,j+1))
            channel_interleaver_pattern(k+1) = v(i+1,j+1);
            k=k+1;
        end
    end
end