function update_bit(row, col)
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

global bits
global bits_updated

offset = size(bits,1)/2^(size(bits,2)-col+1);

for l = 1:size(bits,3)
    if mod(row-1,2*offset) >= offset
        if ~bits_updated(row,col-1)
            update_bit(row,col-1);
        end
        bits(row,col,l) = bits(row,col-1,l);
    else
         if ~bits_updated(row,col-1)
            update_bit(row,col-1);
        end
        if ~bits_updated(row+offset,col-1)
            update_bit(row+offset,col-1);
        end
        bits(row,col,l) = mod(bits(row,col-1,l)+bits(row+offset,col-1,l),2);
    end
end
bits_updated(row,col) = true;

   

