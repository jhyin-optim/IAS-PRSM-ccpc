
%% Generate problem data--  
clc;
clear all;

% �趨������ӡ�
clear;
seed = 97006855;
ss = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(ss);

data = [1024 2048;2000 3000;2000 4000;3000 5000];

[row,~] = size(data);
  
%% ========= the following needs to be tuned and changed=========
paras.beta = 1e-3;  % penalty value
paras.opt = 2.179528e-02;
paras.gamma = 1e-5;
paras.theta_x = 0.1;
paras.theta_y = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state',0); 
paras.lambda_1   = 1e-5 ; % 1-norm's penalty 
paras.MAX_ITER = 1000000;  % outer max iter 
paras.mini_batch = 1;     % number of choosing samples 


% const_batch_size = 1 => mini_batch does not change, M_k may increase ;
% const_batch_size = 0 => mini_batch may increase, M_k does not change .
paras.const_batch_size = 1 ;

paras.batchsize  = 200; 

paras.Time_Budget = 120;

%% Solve problem--------------------------------------
 Tk= 10;%20;  
  
%% establish Latex table
titles = {'(m,n)=(1024,2048)','(m,n)=(2000,3000)','(m,n)=(2000,4000)','(m,n)=(3000,5000)'};

for index=1:row 
    m = data(index,1);
    n = data(index,2); 
    
for i=1:4

    p = 10/n;
    u = sprandn(n, 1, p);
    A = randn(m, n);
    b = A * u + 0.01*randn(m,1);

    %-------computing Lipschitz constant para.v----------
    Ab = A.*b ;
    AtA = full(Ab'*Ab); % Ab'*Ab ; 
    [~,bb]=eig(AtA);    bb=diag(bb); paras.v= sqrt(max(bb));
    
    %% start comparison %% 
  Averg_k3 =0;   Averg_Obj_val3 = 0;   Averg_Equ_err3 = 0; 
 for j= 1: Tk
     [k3, x3, y3, tim, f, history3] =   IAS_PRSM_ccpc(A, b, paras); %all samples 
      
     % compute the point in the average x-axis
     It_cpu3 = [0,history3.cpu];
     It_obj3 = [f,history3.obj];
     It_equ3 = [0,history3.equ];
     %----------------------------------------
     Nit_cpu = 0:1:paras.Time_Budget;  
     Nit_obj3 = interp1(It_cpu3, It_obj3 , Nit_cpu,'linear');%% linear,spline, bubic
     Nit_equ3 = interp1(It_cpu3, It_equ3 , Nit_cpu,'linear');
     
     Averg_k3   = Averg_k3 + k3;
     Averg_Obj_val3 =  Averg_Obj_val3 + Nit_obj3;
     Averg_Equ_err3 =  Averg_Equ_err3 + Nit_equ3;
 %end
 Averg_k3   = round(Averg_k3/Tk);
 Averg_Obj_val3 =  Averg_Obj_val3/Tk; 
 Averg_Equ_err3 =  Averg_Equ_err3/Tk;
 
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Averg_k2 =0;   Averg_Obj_val2 = 0;   Averg_Equ_err2 = 0; 
 %for j= 1: Tk
     [k2, x2, y2, tim, f, history2] =   SAS_ADMM(A, b, paras); %all samples 
      
     % compute the point in the average x-axis
     It_cpu2 = [0,history2.cpu];
     It_obj2 = [f,history2.obj];
     It_equ2 = [0,history2.equ];
     %----------------------------------------
     Nit_cpu = 0:1:paras.Time_Budget;  
     Nit_obj2 = interp1(It_cpu2, It_obj2 , Nit_cpu,'linear');%% linear,spline, bubic
     Nit_equ2 = interp1(It_cpu2, It_equ2 , Nit_cpu,'linear');
     
     Averg_k2   = Averg_k2 + k2;
     Averg_Obj_val2 =  Averg_Obj_val2 + Nit_obj2;
     Averg_Equ_err2 =  Averg_Equ_err2 + Nit_equ2;
% end
 Averg_k2   = round(Averg_k2/Tk);
 Averg_Obj_val2 =  Averg_Obj_val2/Tk; 
 Averg_Equ_err2 =  Averg_Equ_err2/Tk;
 
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Averg_k1 =0;   Averg_Obj_val1 = 0;   Averg_Equ_err1 = 0; 
 %for j= 1: Tk
     [k1, x1, y1, tim, f, history1] = AS_ADMM_modify(A, b, paras);% Algor 1:1~n-1 
  
     % compute the point in the average x-axis
     It_cpu1 = [0,history1.cpu];
     It_obj1 = [f,history1.obj];
     It_equ1 = [0,history1.equ];
     %--------------------------------------
     Nit_cpu = 0:1:paras.Time_Budget;  
     Nit_obj1 = interp1(It_cpu1, It_obj1 , Nit_cpu,'linear');%% linear,spline, bubic
     Nit_equ1 = interp1(It_cpu1, It_equ1 , Nit_cpu,'linear');
     
     Averg_k1   = Averg_k1 + k1;
     Averg_Obj_val1 =  Averg_Obj_val1 + Nit_obj1;
     Averg_Equ_err1 =  Averg_Equ_err1 + Nit_equ1;
 end
 Averg_k1   = round(Averg_k1/Tk);
 Averg_Obj_val1 =  Averg_Obj_val1/Tk; 
 Averg_Equ_err1 =  Averg_Equ_err1/Tk;
 
     %fprintf('Norm_IAS=%.2e, Itr_IAS=%d, Norm_SAS=%.2e, Itr_SAS=%d\n',out1.objec(end),NI1,out2.objec(end),NI2);
     %fprintf('Tcpu_ILM=%.3f, NCG_ILM=%d, Tcpu_OILM=%.3f, NCG_OILM=%d\n',T1,NinItr1,T2,NinItr2);
     %fprintf(fid_tex,'%s & %d & %d & %d/%d/%.3e & %.3f & %d/%d/%.3e & %.3f\\\\ \r\n',...
                %dataset{i},m,n,NI1,NinItr1,out1.objec(end),T1,NI2,NinItr2,out2.objec(end),T2);
     % plot
  h = figure(i);

   semilogy(Nit_cpu,  max(Averg_Obj_val1,Averg_Equ_err1), 'r-*',...
            Nit_cpu,  max(Averg_Obj_val2,Averg_Equ_err2), 'g-d',...
            Nit_cpu,  max(Averg_Obj_val3,Averg_Equ_err3), 'b-.>');
   xlabel('CPU time');ylabel('Opt\_err'); 
   legend('AS-ADMM','SAS-ADMM','IAS-PRSM-ccpc')
   title(titles{i})
end
end
 
 
