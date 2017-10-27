function update_llr(row, col)
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.


global llrs
global llrs_updated
global bits
global bits_updated

offset = size(bits,1)/2^(size(bits,2)-col);

for l = 1:size(bits,3)
    if mod(row-1,2*offset) >= offset
        if ~bits_updated(row-offset,col)
            update_bit(row-offset,col);
        end
        if ~llrs_updated(row-offset,col+1)
            update_llr(row-offset,col+1);
        end
        if ~llrs_updated(row,col+1)
            update_llr(row,col+1);
        end
        llrs(row,col,l) = (-1)^bits(row-offset,col,l)*llrs(row-offset,col+1,l)+llrs(row,col+1,l);
    else
        if ~llrs_updated(row,col+1)
            update_llr(row,col+1);
        end
        if ~llrs_updated(row+offset,col+1)
            update_llr(row+offset,col+1);
        end
        llrs(row,col,l) = minstar(llrs(row,col+1,l),llrs(row+offset,col+1,l));   
    end
end
llrs_updated(row,col) = true;
