%% Generate problem data--  
clc;
clear all;

load('mnist_F.mat');    
load('mnist.mat');  
paras.batchsize  = 200;   

A = double(X); 
b = double(d);      % label         

%-------computing Lipschitz constant para.v----------
Ab = A.*b ;
AtA = full(Ab'*Ab); % Ab'*Ab ; 
[~,bb] = eig(AtA);    bb = diag(bb); paras.v = sqrt(max(bb));

%% ========= the following needs to be tuned and changed=========
paras.opt = 2.179528e-02; 
paras.beta = 1e-3;
paras.delta = 1e-5;
paras.theta_x = 0.1;
paras.theta_y = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state',0); 
paras.lambda_1   = 1e-5 ; % 1-norm's penalty 
paras.MAX_ITER = 1000000;  % outer max iter 
paras.mini_batch = 1;     % number of choosing samples 
paras.barA =  F;   %-------------------------------------------

% const_batch_size = 1 => mini_batch does not change, M_k may increase ;
% const_batch_size = 0 => mini_batch may increase, M_k does not change .
paras.const_batch_size = 1 ;

paras.Time_Budget = 200;

%% Solve problem--------------------------------------
 Tk = 10;%20;  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Averg_k3 = 0;   Averg_Obj_val3 = 0;   Averg_Equ_err3 = 0; 
 for j = 1: Tk
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
 end
 Averg_k3   = round(Averg_k3/Tk);
 Averg_Obj_val3 =  Averg_Obj_val3/Tk; 
 Averg_Equ_err3 =  Averg_Equ_err3/Tk;
 
   Averg_k2 = 0;   Averg_Obj_val2 = 0;   Averg_Equ_err2 = 0; 
 for j= 1: Tk
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
 end
 Averg_k2   = round(Averg_k2/Tk);
 Averg_Obj_val2 =  Averg_Obj_val2/Tk; 
 Averg_Equ_err2 =  Averg_Equ_err2/Tk;
 
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Averg_k1 =0;   Averg_Obj_val1 = 0;   Averg_Equ_err1 = 0; 
 for j= 1: Tk
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

 
 
