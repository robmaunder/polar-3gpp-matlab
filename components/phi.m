function PM_i = phi(PM_iminus1,L_i,u_i)
% Implements equations (11b) and (12) from Balatsoukas-Stimming...
% http://dx.doi.org/10.1109/TSP.2015.2439211
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

global approx_minstar
if approx_minstar
    PM_i = PM_iminus1;
    flags = 0.5*(1-sign(L_i)) ~= u_i;
    PM_i(flags) = PM_i(flags) + abs(L_i(flags));   
else
    PM_i = PM_iminus1 + log(1+exp(-(1-2*u_i)*L_i));
end