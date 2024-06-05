%Calculate the field distribution that results when focussing a
%radially polarized Bessel-Gauss beam from a material with
%refractive index 1.518 into air. Follows:
%D. P. Biss and T. G. Brown, �Cylindrical vector beam focusing
%through a dielectric interface,� Optics Express 9, 490-7 (2001).
%
%and
%
%A. van de Nes, L. Billy, S. Pereira, and J. Braat, �Calculation of
%the vectorial field distribution in a stratified focal region of a
%high numerical aperture imaging system,� Optics Express 12,
%1281-1293 (2004).
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

% clear variable
clear all;
close all;
% add path to utility functions
addpath("../utils/")

% Parameters

% whether to recompute the results or not
% 0 means data is loaded from existing file,
% 1 means it is recomputed
recompute = 1;

% wavelength
lambda_wavelength = 400e-9;

% set up grid upon which to evaluate the field
num_grid_points = 101;
length_grid = 8 * lambda_wavelength;
x_grid = linspace(-length_grid / 2, length_grid / 2, num_grid_points);
y_grid = 0;
z_grid = linspace(-length_grid / 2, length_grid / 2, num_grid_points);

[x_ndgrid, y_ndgrid, z_ndgrid] = ndgrid(x_grid, y_grid, z_grid);
grid_vertices = [x_ndgrid( :), y_ndgrid( :), z_ndgrid( :)];

% refractive index and height vectors for the multilayer
% a single interface
n_vector = [1.518, 1.0];
height_vector = 0;

% numerical aperture
numerical_aperture = 1.4;

% number of integration points along theta and phi
num_points_theta = 500;
num_points_phi = 500;

% output file
output_directory = "artefacts_output_data/";

if ~exist(output_directory)
    mkdir(output_directory);
end;

output_filename = output_directory + "efields_lambda-" + lambda_wavelength + ...
                  "-num_grid_points-" + num_grid_points;
disp(output_filename);

if recompute
  [efield_forward, efield_backward] = focstratfield_general_pol(grid_vertices, n_vector, height_vector, numerical_aperture, lambda_wavelength, num_points_theta, num_points_phi, @besselgauss);
    save(output_filename, "efield_forward", "efield_backward")
else
    load output_filename;
end

% total field is sum of forward and backward fields
efield_total = efield_forward + efield_backward;

% plot results
figure(1); clf;
imagesc(x_grid/lambda_wavelength, z_grid/lambda_wavelength,...
    log10(abs(squeeze(reshape(efield_total(:,3),size(x_ndgrid))).').^2));
ca=caxis;
caxis([ca(2)-4 ca(2)]);
caz=caxis;
axis equal;
colormap hot;

figure(2);clf;
imagesc(x_grid/lambda_wavelength, z_grid/lambda_wavelength,...
    log10(abs(squeeze(reshape(efield_total(:,1),size(x_ndgrid))).').^2));
ca=caxis;
caxis([ca(2)-4 ca(2)]);
caxis(caz);
axis equal;
colormap hot;
