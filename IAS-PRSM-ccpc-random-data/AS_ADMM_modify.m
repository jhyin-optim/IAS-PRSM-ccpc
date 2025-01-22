function [k,x_ag,y_ag, tim, obj_err0, historz] = AS_ADMM_modify(A, b, paras)
%%%  %%%%%%%%%%%  Please cite: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%
%  Bai J., Hager W., Zhang H., Accelerated Stochastic ADMM for    
%  Separable Composite Convex Optimization, written in Agu 8,2019
%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%

%% Global parameters and defaults   
[n,m]  = size(A);  
%[~,m1] = size(b);
Ab     = A.*b ;
Abfu = -Ab ;
QUIET  = 0;
print_time = 0 ;
paras.barA =  A;
 
%%

const_batch_size = paras.const_batch_size ; 
s    = 1.618;     % stepsize of dual variable
c1   = 1;         %>0 constant in Collary 3.11
c2   = 1.e-2;     %>0                                            
v    = paras.v ;   % Lipschitz constant 
beta = 0.04;
mini_batch = paras.mini_batch ;

const_sbet = s*beta;
soft_para  = paras.lambda_1/beta;

M = 200 ;
temp_sigma = 2.e-5 ;

% initialization %
x  = zeros(m,1); y = zeros(n,1); lambda = zeros(n,1);   
x_ag =x; y_ag = y; 
f0 = Func_EvaluateAb(Ab, x, y, paras) ; 
gkeep_cal = 0 ;  erg_k = 0 ; time_spend = 0 ;
if ~QUIET 
    fprintf(' Start function value: %e\n', f0);
end
F_opt = paras.opt;             
obj_err0 = abs(F_opt-f0)/max(abs(F_opt),1);   

x_breve = x;
barAx = paras.barA*x ;
diff_xy = barAx - y;   %============================================
 
%% initial values to determine \C{M}_k
%rho = beta*max(diag(paras.barA'*paras.barA));      %---- >0 
rho = beta ;
eta = 1.1;      %---- >1 
rho_min= 1e-5;  %---- >0 

%t_start   = tic; 
for k = 1:paras.MAX_ITER
    tic;  
    x_old   = x;

    if const_batch_size == 1    
        c_t = c2*(1+k)^1.1 ;
        M_k = ceil(max(c_t,M)) ;      
        eta_k = min( 1./(M_k*(M_k+1)*v),1/(2*v));
        const_eta_k = 2/eta_k; 
    else
        c_t        = c2*(1+k)^1.1 ;
        mini_batch = min (max(paras.mini_batch, floor(c_t)), n) ; 
        M_k = M;
        eta_k       = min( 1./(M_k*(M_k+1)*v),1/(2*v));
        const_eta_k = 2/eta_k;         
    end
    
    rho_xold_h   = rho*x_old + paras.barA'*(lambda - beta* diff_xy) ;%==============
    for t = 1:M_k    
        beta_t   = 2/(t+1);  
        gamma_t  = const_eta_k/t;  
        tempt_bx = (1 - beta_t)*x;
        x_hat  = beta_t*x_breve + tempt_bx;  
        
       %% Generate integers drawn uniformly from 1:n  
        if (mini_batch == n)
            d_t = Grad_EvaluateAb(Ab, x_hat) ;
        else
            i  = randperm(n, mini_batch);
            Ai = Ab(i,:);
            
            if gkeep_cal == 1
                temp1   = 1./(1.+1./exp(-Ai*x_hat));
                temp2   = tkeep(i) ;
                if  mini_batch == 1
                    d_t = -Ai'*(temp1-temp2) + gkeep ; 
                else
                    d_t = -Ai'*(temp1-temp2)/mini_batch + gkeep ; 
                end
            else
                temp   = exp(-Ai*x_hat);
                if  mini_batch == 1                
                    d_t = -Ai'*((1./(1.+1./temp))) ; 
                else
                    d_t = -Ai'*((1./(1.+1./temp)))/mini_batch ;                    
                end                
            end
        end
        
       %% ------------ update x-subproblem -------------
        temp_gs = gamma_t*temp_sigma*x_breve + rho_xold_h - d_t; 
        x_breve = temp_gs/(gamma_t*temp_sigma + rho); 
        x       = beta_t*x_breve + tempt_bx;
    end
        
    %% y-update
    barAx_old = barAx ;
    barAx = paras.barA*x ;
    soft_y = barAx - lambda/beta;  %===========================  
    y      = wthresh(soft_y,'s', soft_para);
%    y      = W_thresh(soft_y,'s', soft_para);

    %% lambda(lagrange multiplier)-update  %          
    diff_xy = barAx - y;        %===============================      
    lambda  = lambda - const_sbet* diff_xy;
   
    %% update rho_k  %
    diff_x = x - x_old ;
    deta_1k  = norm( diff_x )^2;  
    deta_2k  =  norm(barAx-barAx_old)^2; %=====================
    if ( (rho*deta_1k - beta*deta_2k) < 0 )
        rho_min = eta * rho_min;
    end
    rho_BB= beta* deta_2k/deta_1k;
    rho   = max(rho_min, rho_BB);
    
    gkeep_cal = 0 ;     
    if  time_spend > (paras.Time_Budget/3) 
        erg_k = erg_k + 1 ;
        alpha_k = 1/erg_k;       
        x_ag = (1- alpha_k)*x_ag + alpha_k*x ;    
        y_ag = (1- alpha_k)*y_ag + alpha_k*y ;        
        if M_k/m > 1
%           gkeep= Grad_EvaluateAb(Ab, x_ag) ;   
            tkeep = 1./(1.+1./exp(Abfu*x_ag));
            gkeep = (1/n)*Abfu'*tkeep ;         
            gkeep_cal = 1 ;
       end       
    end  
       
   %% diagnostics, reporting, termination checks
    tt =  toc;
    if k==1
        historz.cpu(k) =  tt;
    else
        historz.cpu(k) =  historz.cpu(k-1)+ tt;
    end
    time_spend =  historz.cpu(k) ;  
    
% if erg_k == 1
%     fprintf(' start ergodic k: %3d\n', k) ;
%  end   
  
    if erg_k == 0
        x_ag = x ;
        y_ag = y ;        
    end
    barAx_ag = paras.barA*x_ag ;
    fk = Func_EvaluateAb(Ab, x_ag, barAx_ag, paras) ; 
    historz.equ(k) = norm(paras.barA*x_ag - y_ag);         
    historz.obj(k) = abs(F_opt-fk)/max(abs(F_opt),1);   
    
%     diff_time = historz.cpu(k) - print_time ;
%     if ~QUIET && diff_time >= 1.
%         ttt = Func_EvaluateAb(Ab,x,barAx,paras);
%         tt1 = norm(diff_xy); 
%         fprintf('%3d\t%e\t%e \t%e\t%e\t%e\n', k, historz.equ(k), historz.obj(k),ttt,tt1,historz.cpu(k));
%         print_time = historz.cpu(k) ;
%     end 
    
   % tim = toc(t_start) ;
    tim = historz.cpu(k) ;
    if ( tim >= paras.Time_Budget )
        break;
    end 

end
fprintf('Fk: %e  Euq_err: %e\n',fk,historz.equ(k));
fprintf('time: %4.4f\n',tim);
end


