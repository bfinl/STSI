%% Initialize and load data file
clear; close all; clc
subjID = 'Example_Patient1';
startup_patient;

cd(GridLoc_Folder)
Elec_loc_file = dir('*.xyz');
Elec_loc = readlocs(Elec_loc_file.name);

cd(Sz_Folder)
load([subjID '_Seizure_1.mat']);


Sig_st = 400;
Sig_end = 1600;
Noise_st = 1;
Noise_end = 200;

SigTimeStart = 8000;
Signal_Select = SzData(:,SigTimeStart:SigTimeStart+2000);



cd(Code_Folder)


sRate = 500;
channel = 17;
sig = Signal_Select(channel,:);

[wt,F] = cwt(sig,sRate,'FrequencyLimits',[0 50],'VoicesPerOctave',10);

S = abs(wt);

S_complex = wt;


Tensor = zeros(size(Signal_Select,1),size(Signal_Select,2),40);
sRate = 500;
inputcase = 2; % 1 is the amplitude , 2 is the real part, 3 is the imaginary part, 4 is complex tensor
for channel = 1:size(Signal_Select,1)
     sig = Signal_Select(channel,:);
    
                [wt,F] = cwt(sig,sRate,FrequencyLimits=[0 50]);
                S = abs(wt);
                S_complex  = wt;
             
                switch inputcase
                    case 1
                        S = sqrt(squeeze(S));
                    case 2 
                        S = squeeze(real(S_complex));
                    case 3
                        S = squeeze(imag(S_complex));
                    case 4
                        S = S_complex;
                end
               
               
    for freq = 1:length(F)
    
        for time = 1:size(Signal_Select,2)
            
            Tensor(channel,time,freq) = S(freq,time);
        end
    end
end


switch inputcase
    case 1
        R  = 10;
    case 2 
        R  = 10;
    case 3
        R  = 10;
    case 4
        R  = 10;
end



init = @cpd_gevd;
rng('default');
rng(0)

Uhat = cpd(Tensor,R,'Initialization', init);

%=====%
Topo_ICA = Uhat{1};
Act = Uhat{2};
Spectrum = Uhat{3};

if(inputcase ==  4)
    Topo_ICA = real(Topo_ICA);
    Act = real(Act);
    Spectrum = real(Spectrum);
end

maxfreqs = zeros(R,1);
for r = 1:R
    [~,maxfreqindex] = max(abs(Spectrum(:,r)));
    maxfreqs(r) = F(maxfreqindex);
end
[B,I_sorted] = sort(maxfreqs);
Topo_ICA = Topo_ICA(:,I_sorted);
Act = Act(:,I_sorted);
Spectrum = Spectrum(:,I_sorted);
% save all the components into figures
cd(Fig_Folder)

for r = 1:R
    fig1 = figure('visible','off');

    subplot (1,3,1); 
    topoplot(Topo_ICA(:,r),Elec_loc,'plotrad',0.55,'electrodes','on'); 
    
    subplot (1,3,2);
    plot(Act(:,r)); 
    set(gca,'FontSize',10)
    
    subplot (1,3,3);
    plot(F,Spectrum(:,r)); 
    set(gca,'FontSize',10)
    
    [~,maxfreqindex] = max(abs(Spectrum(:,r)));

    title(['main freq:' num2str(F(maxfreqindex))]);
    fig1.Position = [100 100 1000 400];
    
    saveas(fig1,[num2str(r) '.jpg']);

    % Use below for higher resolution
    %exportgraphics(fig1,[num2str(r) '.jpg'],'Resolution',300);
end

%%
%=====Reconstruct the spatialtemporal signal using certain components====%

component_selected = [3];
freq_select = [34];

Uhat_approx1 = Topo_ICA(:,component_selected); % spatial
Uhat_approx2 = Act(:,component_selected); % temporal
Uhat_approx3 = Spectrum(:,component_selected); % spectral

Uhat_approx2 = Uhat_approx2(:,:);
Uhat_approx3 = Uhat_approx3(freq_select,:);
Uhat_approx = {Uhat_approx1,Uhat_approx2,Uhat_approx3};

Tensor_Groundtruth = Tensor(:,:,freq_select);

Num_TBF = length(component_selected);

cd(GridLoc_Folder);

load('UnconstCentLFD.mat')


Norm_K = svds(K,1)^2;
load('Gradient.mat')
load('Laplacian.mat')
load('NewMesh.mat')
load('LFD_Fxd_Vertices.mat')
load('Edge.mat')
V  =  kron(V,eye(3));
load('Neighbor.mat')
Number_dipole = numel(K(1,:));
Location  = New_Mesh(1:3,:);
TRI = currytri.';
Vertice_Location = curryloc;
Number_edge = numel(V(:,1));

% Select an extended source and display it on the cortex
%  does not really have too much meaning
CortexViewInit;

figure
h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(0)); colorbar;
set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
hold on
hold off
light_position = [3 3 1];
light('Position',light_position);
light_position = [-3 -3 -1];
light('Position',light_position);
colorbar;
x_max = 1;
%x_max = 0.05;
caxis([-x_max x_max]);
view([-1,0.25,0.15])
grid off
axis off
%colorbar('off')
colormap(cmap)


position = Vertice_Location(:,100);
J = zeros(Number_dipole,1);
dist = norms(Location - position);
index_dist = find(dist<20);

J_abbas                 = reshape(repmat(J,1,Num_TBF), [3, Number_dipole/3, Num_TBF]);
J_abbas(:,index_dist,:) = 1;
J_abbas_2 = squeeze(reshape(J_abbas,[Number_dipole,1,Num_TBF]));
% % number of J is equal to number of triangles
J_abbas_1               = squeeze(norms(J_abbas,1));
% %================%
% J_abbas_1               = J_abbas_1.^2;%*var(TBF(:,:),[],2); % For this part can select a smaller time range
J_abbas_1 = sum(J_abbas_1.^2,2);
% %================%
J_init_plot                  = (J_abbas_1);

figure
h1 = trisurf_customized(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_init_plot)); colorbar;
set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
hold on
hold off
light_position = [3 3 1];
light('Position',light_position);
light_position = [-3 -3 -1];
light('Position',light_position);
colorbar;
x_max = max(sum(abs(J_init_plot),2));
if(size(J_init_plot,1) == 1)
    x_max = max(abs(J_init_plot));
end
%x_max = 0.05;
caxis([-x_max x_max]);
view(perspect)
grid off
axis off
%colorbar('off')
colormap(cmap)

%


Uhat_J_ini = J_abbas_2;
sum0 = sum(Uhat_J_ini,2);
index_choose_SourceSpatial = find(sum0~=0);

%Tensor;
N_SourceSpatial = Number_dipole; % Number of dipoles, if unconstrained, then is 3*number dipole
N_Temporal = size(Uhat_approx2,1); % Number of time points
N_Spectral = size(Uhat_approx3,1); % Number of frequencies included in the analysis
N_channel = size(Signal_Select,1);
SensorSpaceTensor = zeros(N_channel,N_Temporal,N_Spectral);

J1 = kron(eye(N_Spectral),J_abbas_2);
TBF = [];
for nn1 = 1:N_Spectral
    TBF = [TBF;Uhat_approx3(nn1,:)'.*Uhat_approx2'];
end
JA = J1*conj(TBF);

K1 = kron(eye(N_Spectral),K);
KJA = K1*JA;

KJA1 = KJA(1:N_channel,:);
KJA1_true= SensorSpaceTensor(:,:,1);

TBF_tensor = zeros(Num_TBF,N_Temporal,N_Spectral);
for nn1 = 1:N_Spectral
    TBF_tensor(:,:,nn1) = Uhat_approx3(nn1,:)'.*Uhat_approx2';
end

Phi_tensor = Tensor_Groundtruth;

num_it_new  = 0;

num_it_x    = 25;
num_it_y    = 25;
Number_iteration = 4;

parms = cell(N_Spectral,1);

for iter_TBF = 1:N_Spectral
    Phi = squeeze(Phi_tensor(:,:,iter_TBF));
    TBF = squeeze(TBF_tensor(:,:,iter_TBF));
    
    Phi_noisy = Phi(:,:);
    
    [Number_sensor,Number_Source] = size(Phi_noisy);
    
    % Estimating Noise and Initializing
    

    N_wht                  = Phi_noisy - Phi_noisy*TBF.'*pinv(TBF*TBF.')*TBF;
    Noise_only             = Phi_noisy(:,Noise_st:Noise_end);
    SNR_unused             = norm(Phi_noisy,'fro')^2/norm(Noise_only,'fro')^2; % a SNR that is not used
    
    % Noise Estimation try different options
    Sigma_inv_half = 1*diag(1./std(Noise_only(:,:),[],2));
    %==============%
    
    Sigma_inv_half = diag(repmat(mean(diag(Sigma_inv_half)),[Number_sensor 1]));
    
    %==============%
    Phi_noisy = Phi_noisy(:,Sig_st:Sig_end);
    TBF       = TBF (:,Sig_st:Sig_end);
    psi       = Phi_noisy*(TBF.')/(TBF*(TBF.'));
    
    % Initializing the Solution
    SNR_avg  = mean(std(Sigma_inv_half*Phi_noisy(:,1:end),[],2).^2)/mean(std(Sigma_inv_half*Noise_only,[],2).^2)
    D_W      = ((norms(K+10^-20).^-1));
    K_n      = (Sigma_inv_half*K).*repmat(D_W,[Number_sensor 1]);
    psi_n    = Sigma_inv_half*psi;
    [~, gam] = eig((Sigma_inv_half^1*K)*((Sigma_inv_half*K).'));
    %======%
    Smf      = 10^3; % smooth factor
    %======%
    J_ini    = (K_n.'*pinv(K_n*K_n.' + mean(diag(gam))*(Smf/SNR_avg)*eye(Number_sensor))*psi_n).*repmat(D_W',[1 Num_TBF]);
    %J_ini = randn(size(J_ini))*norm(J_ini,1);
   
    close all
    [Number_sensor,Number_Source] = size(Phi_noisy);
    epsilon = 0.1;   
    
    % Noise Power Parameter (Chi_2 Dist Theory)
    prob = 0.95;
    beta = icdf('chi2',prob,Number_sensor-1);  % Spike Avg
    
    Phi_norm = Sigma_inv_half*(Phi_noisy(:,:));
    K_norm   = Sigma_inv_half*(K);
    
    power = sum((Sigma_inv_half*Noise_only).^2,1);
 
    power_signal = sqrt(sum(Phi_noisy.^2,2))/size(Phi_noisy,2);
    power_noise = sqrt(sum(Noise_only.^2,2))/size(Noise_only,2);
    SNR_allchannel = power_signal./power_noise;
    SNR = mean(SNR_allchannel);
    disp(['SNR of the data is ' num2str(SNR)]);
    
    [X,B] = hist(power,100);
    % mean(power);
    % mean(power) + std(power);
    % figure; histogram(power,100)
    Sum_X = cumsum(X)/sum(X);
    % figure; plot(Sum_X)
    ind_90 = find(Sum_X > 0.90); B(ind_90(1,1));
    ind_95 = find(Sum_X > 0.95); B(ind_95(1,1));
    ind_50 = find(Sum_X > 0.50); B(ind_50(1,1));
    ind_30 = find(Sum_X > 0.30); B(ind_30(1,1));
    CF_min = beta/(norm(Phi_norm,'fro')^2/Number_Source); % norm(Phi_norm,'fro')^2/Number_Source = average power per time point 
    % for noise power as noise might not be whitened enough
    %==============%
    CF     = max(beta/B(ind_30(1,1)),CF_min) ;
    %==============%
    
    %=============%
    alpha_vec = 0.1;
    Num_alpha = numel(alpha_vec);
    %=============%
    % Variables to save results
    J_sol = zeros(Number_dipole,Num_TBF,Number_iteration,Num_alpha);
    Y_sol = zeros(Number_edge,Num_TBF,Number_iteration,Num_alpha);
    W_sol = zeros(Number_dipole,Num_TBF,Number_iteration,Num_alpha);
    W_d_sol = zeros(Number_edge,Num_TBF,Number_iteration,Num_alpha);
    J_sol_2 = J_sol;
    Y_sol_2 = Y_sol;
    
    % Weighting the depth - initializing weighting matrix
    M = size(K,2);
    N = size(V,1);
    W = ones(M,Num_TBF);
    W_d = ones(N,Num_TBF);
    
    % Initialize optimization parameters - set number of iterations, etc.
    Lambda_min = Norm_K*sqrt(Number_sensor)/max(norms(K.'*Phi_noisy));
    
    C_t         = (TBF*(TBF.'));
    TBF_norm    = TBF.'*(C_t\TBF);
    T_norm      = (TBF.')/C_t;
    alpha       = alpha_vec;
    lambda      = 6*max(1, Lambda_min);
    W_v         = (V.')*V;
    L_v         = 1.1*svds(W_v,1);
    x_0         = zeros(Number_dipole,Num_TBF);
    eps         = 10^-3;
    
    betta_tilda = (norm(Phi_norm,'fro')^2/SNR_avg)/CF;

    %betta_tilda = beta * Number_Source; %theoretical 
    
    Abbas       = K_norm*(K_norm.');
    [U, D, ~]   = svd(Abbas);
    K_U         = U.'*K_norm;
    
    parms{iter_TBF}.lambda = lambda;
    parms{iter_TBF}.betta_tilda = betta_tilda;
    parms{iter_TBF}.Phi_norm = Phi_norm;
    parms{iter_TBF}.TBF = TBF;
    parms{iter_TBF}.TBF_norm = TBF_norm;
    parms{iter_TBF}.T_norm = T_norm;
    parms{iter_TBF}.U = U;
    parms{iter_TBF}.D = D;
    parms{iter_TBF}.K_norm = K_norm;
    parms{iter_TBF}.K_U = K_U;

end


num_it_new  = 0;
stop_crt    = 10^-3;

% Initialize the optimization problem, Loop through till convergence,
% update the weights using the iterative re-weighting schema and solve
% agian until convergence or counter overflow.

stop_itr = 0;
weight_it = 0;
max_weight_itr_num  = Number_iteration;
t1 = tic;
res_obj_combined = cell(0);
for i_alpha = 1:Num_alpha
    alpha = alpha_vec(i_alpha);
    x_0 = J_ini;
    stop_itr = 0;
    weight_it = 0;
    
    W = ones(M,Num_TBF);
    W_d = ones(N,Num_TBF);
    
    while (~stop_itr) 
        weight_it = weight_it + 1;
        if weight_it > max_weight_itr_num
            stop_itr = 1;
        end
        
        [J,Y,res_obj] = FISTA_ADMM_IRES_Tensor (alpha, W_v, V, L_v,x_0, W, W_d, eps,num_it_x,num_it_y,num_it_new,parms);
        res_obj_combined{end+1} = res_obj;
        x_0 = J;
        J_sol(:,:,weight_it,i_alpha) = J;
        Y_sol(:,:,weight_it,i_alpha) = Y;
        
        J_n     = reshape(J, [3, Number_dipole/3, Num_TBF]);
        Ab_J    = squeeze(norms(J_n))+10^-20;
        if(size(Ab_J,1) == 1)
            Ab_J = Ab_J';
        end
        Y_n     = reshape(V*J, [3, Number_edge/3, Num_TBF]);
        Ab_Y    = squeeze(norms(Y_n))+10^-20;
        if(size(Ab_Y,1) == 1)
            Ab_Y = Ab_Y';
        end
        
        W_old   = W;
        W_tr    = 1./(0 + Ab_J./(repmat(max(Ab_J),[Number_dipole/3 1])) + (epsilon+10^-16) );
        W       = circshift(upsample(W_tr,3),[0 0]) + circshift(upsample(W_tr,3),[1 0]) + circshift(upsample(W_tr,3),[2 0]); 
        W_d_old = W_d;
        W_d_tr  = 1./(0 + Ab_Y./(repmat(max(Ab_Y),[Number_edge/3 1]))+(epsilon+10^-16));
        W_d     = circshift(upsample(W_d_tr,3),[0 0]) + circshift(upsample(W_d_tr,3),[1 0]) + circshift(upsample(W_d_tr,3),[2 0]); 

        W_sol(:,:,weight_it,i_alpha) = W;
        W_d_sol(:,:,weight_it,i_alpha) = W_d;
        if (norm(W-W_old)/norm(W_old) < stop_crt && norm(W_d-W_d_old)/norm(W_d_old) < stop_crt)
            stop_itr = 1;
        end
        clc
        figure
                h1 = trisurf_customized(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),sum((Ab_J),2)); colorbar
                set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
                 light_position = [3 3 1];
                 light('Position',light_position);
                 light_position = [-3 -3 -1];
                 light('Position',light_position);
                 colorbar;
                 x_max = max(abs(sum(J,2)));
                 if isnan(x_max)
                     x_max = 1;
                 end
                 caxis([-x_max x_max]);
                 view(232,-5)
                 
                 grid off
             colormap(cmap)
             cd(Fig_Folder)
         
                
                mkdir('STSI');
                cd('STSI')
        
             
         name_fig = ['TBF_1st_Iteration_',num2str(weight_it),'.fig'];
         saveas(gcf,name_fig)
         name_fig = ['TBF_1st_Iteration_',num2str(weight_it),'.jpeg'];
         saveas(gcf,name_fig)
         
    end
    
    close all; 
end
t_elaps = toc(t1);
disp(['total time elapsed in hours: ' num2str(t_elaps/3600)]);





    J                       = squeeze(J_sol(:,:,end,1));
    Num_TBF                 = size(J,2);

    
    J_abbas                 = reshape(J, [3, Number_dipole/3, Num_TBF]);
    % number of J is equal to number of triangles
    J_abbas_1               = squeeze(norms(J_abbas,1));
    %================%
  
    J_abbas_2               = J_abbas_1.^2*sum(squeeze(var(TBF_tensor(),[],2)),2);
    J_abbas_2 = J_abbas_2(:);
    %================%
  
    J_init                  = (J_abbas_2);
    
    figure
    h1 = trisurf_customized(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_init)); colorbar;
    set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
    hold on
    hold off
    light_position = [3 3 1];
    light('Position',light_position);
    light_position = [-3 -3 -1];
    light('Position',light_position);
    colorbar;
    x_max = max(sum(abs(J_init),2));
    if(size(J_init,1) == 1)
        x_max = max(abs(J_init));
    end
    %x_max = 0.05;
    caxis([-x_max x_max]);
   view(232,-5)
    grid off
    axis off
    %colorbar('off')
    colormap(cmap)
    
    % Now instead of plotting the whole J distribution, we find clusters (patches)
    % of J that has magnitude corssing a threshold
      
    J_col                   = J_abbas_2;
    %================%
    Thr                     = 1/64;
    %================%
    J_col(abs(J_col)<Thr*max(abs(J_col)))=0;
    IND                     = Find_Patch(Edge, 1, J_col);
    J_init                  = IND;

    J_abbas_2 = J_abbas_2.*IND;

  
    figure
    h1 = trisurf_customized(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_abbas_2));
    colorbar;
    set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
    light_position = [3 3 1];
    light('Position',light_position);
    light_position = [-3 -3 -1];
    light('Position',light_position);
    colorbar;
    x_max                   = max(sum(abs(J_abbas_2),2));
    caxis([-x_max x_max]);
   view(232,-5)
    grid off
    axis on
    colormap(cmap)

