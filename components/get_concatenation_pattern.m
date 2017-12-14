function concatenation_pattern = get_concatenation_pattern(G, C)
% GET_SEGMENTATION_PATTERN Get the segmentation pattern
%   segmentation_pattern = GET_SEGMENTATION_PATTERN(A, C) obtains the
%   segmentation pattern, as desribed in Section 5.2.1 of TS38.212
%   V1.1.1
%
%   A should be an integer scalar. It specifies the number of bits in the
%   information bit sequence.
% 
%   C should be an integer scalar. It specified the number of code block
%   segments to use.
%
%   segmentation_pattern will be an integer matrix, having C rows and
%   ceil(A/C) columns. The matrix will comprise A unique elements in the
%   range 1 to A, as well as additional zero-valued elements, if necessary
%   to fill the matrix. Each row in the matrix corresponds to a different
%   one of the codeblock segments. Each integer in each row identifies 
%   which one of the A information bits provides the corresponding bit for 
%   the codeblock segment. Zero-valued elements indicate that the
%   corresponding bit for the codeblock segment should be set to zero.
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

if mod(G,C) ~= 0
    error('polar_3gpp_matlab:UnsupportedBlockLength','G should be divisible by C.');
end

E = G/C;

g = 1:G;
concatenation_pattern = zeros(C,E);


k=0;
r=0;

while r<C
    j=0;
    while j<E
        concatenation_pattern(r+1,j+1) = g(k+1);
        k=k+1;
        j=j+1;
    end
    r=r+1;
end


