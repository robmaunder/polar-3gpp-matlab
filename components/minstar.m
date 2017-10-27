function c = minstar(a,b)
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.


global approx_minstar
if approx_minstar || abs(a) == inf || abs(b) == inf
    c = sign(a)*sign(b)*min(abs(a),abs(b));
else
    c = 2*atanh(tanh(a/2)*tanh(b/2));
end

