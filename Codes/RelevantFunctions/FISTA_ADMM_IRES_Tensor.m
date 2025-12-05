% Code to solve the spatio=temporal IRES. Basically a block coordinate
% projected gradient descent algorithm where the edges and the solution
% itself will be updated sequentially.
% Xiyuan Jiang
% 11/17/2022

function [x, y,res_obj,J_intermediate,J_intermediate_count] = FISTA_ADMM_IRES_Tensor(alpha, W_v, V, L_v,x_0, W, W_d, eps,num_it_x,num_it_y,num_it_new,parms)

J_intermediate = {};
J_intermediate_count = [];
lambda = parms{1}.lambda;

Alpha = diag(sum(abs(x_0))/max(sum(abs(x_0))));
y_cond = 1;

y_it = 0;
x_tild = x_0;
x_old = x_0;
y = (V*x_tild);
u = zeros(size(y));

e_rel  = 10^-4;
e_abs  = 10^-8;
tau    = 2;
mu     = 10;
res_obj =[];
while (y_cond)
    
    y_it = y_it + 1;
    %disp(['FISTA solving at step ' num2str(y_it)]);
    
    if (y_it>num_it_y)
        y_cond = 0;
    end
    
    t_old = 1;
    x_it = 0;
    x_cond = 1;
    
    while (x_cond)
        
        x_it = x_it + 1;
        
        if (x_it>num_it_x)
            x_cond = 0;
        end

        %===start looping through all the J corresponding to each freq
        %Solve for the updated J
        x_new_tran = wthresh ((x_tild - (W_v*x_tild-V.'*(y+u))/L_v),'s',(W*Alpha)*(alpha/lambda/L_v));
        
        x_new_er_combined = cell(length(parms),1);

        for iter_TBF_inside = 1:length(parms)
            Phi =  parms{iter_TBF_inside}.Phi_norm;
            K = parms{iter_TBF_inside}.K_norm;
            TBF = parms{iter_TBF_inside}.TBF;
            TBF_norm = parms{iter_TBF_inside}.TBF_norm;
            T_norm = parms{iter_TBF_inside}.T_norm;
            U = parms{iter_TBF_inside}.U;
            betta_tilda = parms{iter_TBF_inside}.betta_tilda;
            D = parms{iter_TBF_inside}.D;
            K_U = parms{iter_TBF_inside}.K_U;

            % Prepare parameters for projection to hyper ellips
            phi_hat = Phi - K*x_new_tran*TBF;
            Phi_i = (norms ( ((U.')*phi_hat*TBF_norm).' ).').^2; 
            Tp_norm = phi_hat*T_norm;
            N_Phi =  phi_hat - phi_hat*TBF_norm;
            N_F   = norm(N_Phi,'fro')^2;
        
            if (norm(phi_hat,'fro')^2<=betta_tilda)
                x_new_er = zeros(size(x_new_tran));
            
            else
                x_new_er = Proj_Flat_Hyper_Ellips(betta_tilda-N_F, Phi_i, phi_hat, U, D, K_U, Tp_norm,num_it_new);
            end
            x_new_er_combined{iter_TBF_inside} = x_new_er;
       % figure;hist(x_new_er(:,1))
       % figure;hist(x_new_er(:,2))
            %===end looping through all the J corresponding to each freq
        end
        
        x_new = x_new_tran;
     
        % x_new  = x_new + x_new_er_combined{2};
        % for iter_TBF_inside = 1:length(parms)
        %       x_new_er_i = x_new_er_combined{iter_TBF_inside};
        %       x_new_er_j = zeros(size(x_new_tran));
        %       x_new_er_j(:,iter_TBF_inside) = x_new_er_i(:,iter_TBF_inside);
        % 
        %      x_new = x_new + x_new_er_j;
        % end
         for iter_TBF_inside = 1:length(parms)
            x_new = x_new + 1/length(parms)*x_new_er_combined{iter_TBF_inside};
        end

        t_new = 0.5*(1+sqrt(1+4*t_old^2));
        
        x_tild = x_new + ((t_old-1)/t_new)*(x_new-x_old);

        if norm(x_new-x_old,'fro')/norm(x_new,'fro') < eps
            x_cond =0;
            
        end
       
        x_old = x_new;
        t_old = t_new;

    end
    
    y_new = wthresh((V*x_new)-u,'s',W_d/lambda);
    
    
    if norm(y-y_new,'fro')/norm(y_new,'fro') < eps
            y_cond =0;
    end
    
    prim_res = norm(y_new - V*x_new,'fro');
    dual_res = norm(lambda*(V.')*(y - y_new),'fro');
    
    e_prim = sqrt(numel(y))*e_abs + ...
        e_rel*max(sum(norms(V*x_new)),sum(norms(y_new)));
    e_dual = sqrt(numel(x_new))*e_abs + ...
        e_rel*sum(norms(V.'*y_new));
    
    if (prim_res <= e_prim && dual_res <= e_dual)
        y_cond = 0;    
    end
        
    y = y_new;
    u = u + y - V*x_new;
    
%     if     (prim_res > mu*dual_res)
%         lambda = lambda*tau;
%         u      = u/tau;
%     elseif (dual_res > mu*prim_res)
%         lambda = lambda/tau;
%         u      = u*tau;
%     end
x11 = wthresh(x_new,'s',(W*Alpha)*(alpha/lambda/L_v));
y11 = wthresh((V*x_new)-u,'s',W_d/lambda);

res_obj = [res_obj, [alpha*sum(abs(W.*x11),'all')+sum(abs(W_d.*y11),'all');median(x11(x11>1e-10))]];

if(mod(y_it,20)==1)
    J_intermediate{end+1} = x11;
    J_intermediate_count = [J_intermediate_count,y_it];
    disp(['FISTA solving at step ' num2str(y_it)]);
end

end


%disp('FISTA solving this round finished');
x = wthresh(x_new,'s',(W*Alpha)*(alpha/lambda/L_v));
y = wthresh((V*x_new)-u,'s',W_d/lambda);
    
    