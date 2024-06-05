%function [wx,wy] = besselgauss(th,ph)
%Copyright (C) 2018 Peter Munro
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU Lesser Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU Lesser Public License for more details.
%    
%You should have received a copy of the GNU Lesser Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [wx,wy] = besselgauss(th,ph)
    th_max = asin(1.4/1.518);
    beta = 1;
    
    L = exp(-beta^2*sin(th).^2/sin(th_max)^2).*besselj(1,2*beta*sin(th)/sin(th_max));
    
    wx = cos(ph).*L;
    wy = sin(ph).*L;
    
    
