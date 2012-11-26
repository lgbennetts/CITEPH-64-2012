function Res_mat = Freq_response_gen_multi(GeomDisks, TankDim, Sample, ...
    Param, Mesh)

% Res_mat = Freq_response_gen(GeomDisks,TankDim,Sample,Param,Mesh)
%
% Solution to the 3D wave tank problem over the frequency range defined by
% the sampling of the input signal. The inverse Fourier transform of the
% quantities computed here will provide the transient response of the
% system to a wave maker forcing. The structure array Res_mat contains the
% quantities of interest to produce results.


%% Initialisation of frequency-dependent variables

%%% Potential amplitudes:
Res_mat.Amp_1P = cell(1,Sample.Nfilt);
Res_mat.Amp_1M = cell(1,Sample.Nfilt);
Res_mat.Amp_3P = cell(1,Sample.Nfilt);
Res_mat.Amp_3M = cell(1,Sample.Nfilt);

%%% Wavenumbers in free water
Res_mat.k0 = cell(1,Sample.Nfilt);

%%% Frequency parameter
Res_mat.kappa = zeros(1,Sample.Nfilt);

%%% Beach reflection coefficient
Rb = beach_TF(Param.Rb, 2*pi*Sample.f);


%% Generation of the frequency response
disp(['Fourier spectrum conmposed of ' int2str(Sample.Nfilt) ...
    'frequencies'])
h_waiting_bar = waitbar(0, ['Fourier spectrum solution of 3D wave ' ...
    'tank problem. Please wait...']);
for i = 2:Sample.Nfilt
    waitbar(i / Sample.Nfilt,h_waiting_bar)
    
    %%% Frequencies
    Res_mat.kappa(i) = (2*pi*Sample.f(i))^2/Param.g;
    
    %%% Solution - disks in a channel
    [Rm,Tm,Rp,Tp,v_vec,u_vec,k0,weight_0,x_lim] = ...
    fn_MultiMode_MultiFloe([Param.rho_0 Param.rho Param.nu Param.E ...
    Param.rho_0], Param.N, Param.Mev, TankDim, Res_mat.kappa(i), ...
    Param.beta.', GeomDisks(:,4).', Param.draft.', GeomDisks(:,3).', ...
    [GeomDisks(:,1).';GeomDisks(:,2).'], Mesh.r_vec, Mesh.th_vec, ...
    Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, Param.res_green, ...
    Param.extra_pts);
end
close(h_waiting_bar)