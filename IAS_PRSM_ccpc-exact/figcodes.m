 %% 
   h = figure;

 semilogy(Nit_cpu,  max(Averg_Obj_val1,Averg_Equ_err1), 'r-*',...
          Nit_cpu,  max(Averg_Obj_val2,Averg_Equ_err2), 'g-d',...
          Nit_cpu,  max(Averg_Obj_val3,Averg_Equ_err3), 'b-.>');
 xlabel('CPU time');ylabel('Opt\_err'); 
 legend('AS-ADMM','SAS-ADMM','IAS-PRSM-ccpc')