function segmentation_pattern = get_segmentation_pattern(A, C)
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

A_prime = ceil(A/C)*C;

a = 1:A;

a_prime = zeros(1,A_prime);
for i=0:A_prime-A-1
    a_prime(i+1) = 0;
end
for i=A_prime-A:A_prime-1
    a_prime(i+1) = a(i-(A_prime-A)+1);
end

segmentation_pattern = zeros(C,A_prime/C);
s=0;
for r=0:C-1
    for k=0:A_prime/C-1
        segmentation_pattern(r+1,k+1) = a_prime(s+1);
        s = s+1;
    end
end

