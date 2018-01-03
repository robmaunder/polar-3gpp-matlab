function Pi = get_3GPP_crc_interleaver_pattern(K)  
% GET_3GPP_CRC_INTERLEAVER_PATTERN Obtain the 3GPP Cyclic Redundancy Check 
% (CRC) interleaver pattern, according to Section 5.3.1.1 of 3GPP TS 38.212 
%   Pi = GET_3GPP_CRC_INTERLEAVER_PATTERN(K) obtains the CRC interleaver
%   pattern Pi.
%
%   K should be an integer scalar. It specifies the number of bits in the 
%   information and CRC bit sequence. 
%
%   Pi will be an integer row vector, comprising K unique elements in the
%   range 1 to K. Each integer identifies which one of the K information or 
%   CRC bits provides the corresponding interleaved bit. Interleaving
%   can be implemented according to c = b(Pi), while
%   deinterleaving is implemented accoring to 
%   b_hat(Pi) = c_hat.
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

% Note that TS38.212 defines a 164-bit interleaver, which can be shortened
% to interleave fewer bits. During the standardisation process, a 224-bit
% interleaver was designed, as is implemented below. When shortened to 164
% or fewer bits, it gives identical interleaver patterns to the 164-bit
% interleaver that was finally agreed.
Pi_IL_max = [0 2 3 5 6 8 11 12 13 16 19 20 22 24 28 32 33 35 37 38 39 40 41 42 44 46 47 49 50 54 55 57 59 60 62 64 67 69 74 79 80 84 85 86 88 91 94 102 105 109 110 111 113 114 116 118 119 121 122 125 126 127 129 130 131 132 136 137 141 142 143 147 148 149 151 153 155 158 161 164 166 168 170 171 173 175 178 179 180 182 183 186 187 189 192 194 198 199 200 1 4 7 9 14 17 21 23 25 29 34 36 43 45 48 51 56 58 61 63 65 68 70 75 81 87 89 92 95 103 106 112 115 117 120 123 128 133 138 144 150 152 154 156 159 162 165 167 169 172 174 176 181 184 188 190 193 195 201 10 15 18 26 30 52 66 71 76 82 90 93 96 104 107 124 134 139 145 157 160 163 177 185 191 196 202 27 31 53 72 77 83 97 108 135 140 146 197 203 73 78 98 204 99 205 100 206 101 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223];   

if K > length(Pi_IL_max)
    error('polar_3gpp_matlab:UnsupportedBlockLength','K should be no greater than 224.');
end

Pi = zeros(1,K);
k=0; 
for m=0:length(Pi_IL_max)-1
    if Pi_IL_max(m+1) >= length(Pi_IL_max)-K
        Pi(k+1) = Pi_IL_max(m+1)-(length(Pi_IL_max)-K)+1; % +1 because index start at 1 in Matlab
        k = k+1;
    end
end