% In this tutorial we will construct the Zener Polaron model for
% Pr(Ca0.1Sr0.9)2Mn2O7 as described in:
% Ground State in a Half-Doped Manganite Distinguished by Neutron Spectroscopy
% G.E. Johnstone, T.G. Perring, O. Sikora, D. Prabhakaran, and A.T. Boothroyd
% Phys. Rev. Lett. 109 237202 (2012), https://arxiv.org/abs/1210.710
%
% The Zener polaron is a strongly coupled pair of Mn atoms which share an
% electron via the Zener double exchange mechanism. They act like a single
% S=7/2 spin, and will be modelled in SpinW as such with the position of
% the Mn_2^7+ polaron being the center of the pair of Mn-Mn atoms.

% Lattice parameters for Pr(Ca0.1Sr0.9)2Mn2O7
lat = [5.408 5.4599 19.266];
alf = [90 90 90];

% Like for the CE model in prcasrmn2o7.m we are forced to use a doubled
% "structural" unit cell to define the model in SpinW (which means that
% hk0 values should be halved to compare between SpinW and the paper).
zp = spinw;
zp.genlattice('lat_const', lat.*[2 2 1], 'angled', alf, 'spgr', 'x,y,-z');
% Note that you will get several warnings that SpinW does not understand
% what a Zener Polaron is and cannot assign it a form factor or scattering
% length - you can ignore this as we will just calculate the dispersion
% with this model.
zp.addatom('label', 'ZP-up', 'r', [0.375 0.375 0.1], 'S', 7/2, 'color', 'gold');
zp.addatom('label', 'ZP-up', 'r', [0.875 0.625 0.1], 'S', 7/2, 'color', 'gold');
zp.addatom('label', 'ZP-dn', 'r', [0.375 0.875 0.1], 'S', 7/2, 'color', 'gold');
zp.addatom('label', 'ZP-dn', 'r', [0.875 0.175 0.1], 'S', 7/2, 'color', 'gold');

% Set the magnetic structure
spin_up = cellfun(@isempty, strfind(zp.table('matom').matom, 'up'));
spin_dn = cellfun(@isempty, strfind(zp.table('matom').matom, 'dn'));
S0 = [1; 0; 0];
spin_up = find(~cellfun(@isempty, strfind(zp.table('matom').matom, 'up')));
spin_dn = find(~cellfun(@isempty, strfind(zp.table('matom').matom, 'dn')));
SS = zeros(3, numel(spin_up)+numel(spin_dn));
SS(:, spin_up) = repmat(S0, 1, numel(spin_up));
SS(:, spin_dn) = repmat(-S0, 1, numel(spin_dn));
zp.genmagstr('mode', 'direct', 'S', SS)

% Now assign the exchange interactions
% Parameters are scaled up from Ewings et al., Phys. Rev. B 94 014405 (2016)
% https://arxiv.org/abs/1501.01148
JFU = -1.6 * 1.5;
JFD = -0.11 * 1.5;
JA = 0.11 * 1.5;
Jperp = 1.05;

zp.gencoupling('forceNoSym', true)
zp.addmatrix('label', 'JFU', 'value', JFU, 'color', 'green') % FM between non parallel polarons in +b
zp.addmatrix('label', 'JFD', 'value', JFD, 'color', 'white') % FM between non parallel polarons in -b
zp.addmatrix('label', 'JA', 'value', JA, 'color', 'yellow') % AFM between parallel polarons
zp.addmatrix('label', 'Jperp', 'value', Jperp, 'color', 'blue')
zp.addcoupling('mat', 'Jperp', 'bond', 1)
zp.addcoupling('mat', 'JA', 'bond', 2)
zp.addcoupling('mat', 'JA', 'bond', 3)
zp.addcoupling('mat', 'JFD', 'bond', 4)
zp.addcoupling('mat', 'JA', 'bond', 5)
zp.addcoupling('mat', 'JFU', 'bond', 6, 'atom', 'ZP-up')
zp.addcoupling('mat', 'JFD', 'bond', 6, 'atom', {'ZP-up', 'ZP-dn'})
zp.addcoupling('mat', 'JFU', 'bond', 8)

zp.table('bond', 1:10)
plot(zp, 'range', [0 2; 0 2; -0.2 0.2]);

%%
% Optimises the structure
res = zp.optmagsteep()
plot(zp, 'range', [0 2; 0 2; 0 0.2])

%%
% Calculates the dispersion

spec = zp.spinwave({[0 0 0] [2 0 0] [2 0 2] [0 0 2] [0 0 0] 500}, 'hermit', false);
figure; sw_plotspec(spec);
specg = sw_egrid(spec, 'Evect', linspace(0,100,2000), 'imagChk', false);
figure; sw_plotspec(specg,'mode','color','dE',0.5);