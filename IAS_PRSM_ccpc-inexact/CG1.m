function [y1,dn1] = CG1(step_par1,step_par2,rho,x,y,x_old,lambda_mid,baru,barv,u_ag,v_ag,vn1,method,paras)
% B = 2*mu*A'tA + (beta + eta1)*I
% b = beta*x + eta1*barv - lambda_mid;  % solve linear equations B*x - b=0
% the residual r(x) = B*x - b;

theta_x  = 0.9;
theta_y  = 0.9;
eta1 = paras.eta1;
beta = paras.beta;
lambda_1 = paras.lambda_1;

osos = (1 - step_par2)^2/(1 + step_par1);
constant_c0 = beta*(2 - step_par1 - step_par2 - osos);
constant_c2 = (1 - step_par1)/(1 + step_par1);
kappa = 1 - constant_c2*theta_y^2;

A = paras.barA;
AtA = A'*A;
r0 = 2*lambda_1*AtA*y + ( beta + eta1)*y - beta*x - eta1*barv + lambda_mid;
p0 = - r0;
r0tr0 = r0'*r0;
y1 = y;

switch method
    case 'IAS_PRSM_ccpc1'
     
    beta1 = paras.beta1;
    beta2 = paras.beta2;
    beta3 = paras.beta3;
    
    dn0 = r0;
    dn1 = r0;
    normdn1 = norm(dn1);
    normxbu_old = norm(x_old - u_ag);
    normxbu = norm(x - baru);
    Axu = A*(x_old - u_ag);
    normxbarA = norm(Axu);
    Axbu = A*(x - baru);
    normxbarAu = norm(Axbu);
    normxbv_old = norm(y1 - v_ag);
    normxbv = norm(y - barv);
   
    normMw = norm(x - y);
    wvtd = (y - vn1)'*dn1;
    yydd = (y - y1 )'*(dn1 - dn0);
    ybyy =  (y1 - v_ag) - (y - y1 );
    normybyy = norm(ybyy);
    awd = abs(wvtd);
    ayd = abs(yydd);
    asn = 2*awd + 2*constant_c2*ayd + beta*normdn1^2;
    btnn = beta1*(theta_x*rho*normxbu_old^2 - beta*theta_x*normxbarA^2 + rho*normxbu^2 - beta*normxbarAu^2 );
    btkkc = beta2*(theta_y*eta1*normxbv_old^2 + kappa*eta1*normxbv^2 + constant_c2*eta1*normybyy^2);
    bcn = beta3*constant_c0*normMw^2;
    bbb = btnn +  bcn + btkkc; 
    
    while asn > bbb  
            Ap0 = A*p0;
            Cp0 = 2*lambda_1*A'*Ap0 + ( beta + eta1)*p0;
            p0Cp0 = p0'*Cp0;
            alphak = r0tr0/p0Cp0;
            y1 = y1 + alphak*p0;
            r0 = r0 + alphak*Cp0;
            r0tr0_old = r0tr0;
            r0tr0 = r0'*r0;
            betak = r0tr0/r0tr0_old;
            p0 = -r0 + betak*p0;
            
    dn0 = r0;
    dn1 = r0;
    normdn1 = norm(dn1);
    normxbu_old = norm(x_old - u_ag);
    normxbu = norm(x - baru);
    Axu = A*(x_old - u_ag);
    normxbarA = norm(Axu);
    Axbu = A*(x - baru);
    normxbarAu = norm(Axbu);
    normxbv_old = norm(y1 - v_ag);
    normxbv = norm(y - barv);
   
    normMw = norm(x - y);
    wvtd = (y - vn1)'*dn1;
    yydd = (y - y1 )'*(dn1 - dn0);
    ybyy =  (y1 - v_ag) - (y - y1 );
    normybyy = norm(ybyy);
    awd = abs(wvtd);
    ayd = abs(yydd);
    asn = 2*awd + 2*constant_c2*ayd + beta*normdn1^2;
    btnn = beta1*(theta_x*rho*normxbu_old^2 - beta*theta_x*normxbarA^2 + rho*normxbu^2 - beta*normxbarAu^2 );
    btkkc = beta2*(theta_y*eta1*normxbv_old^2 + kappa*eta1*normxbv^2 + constant_c2*eta1*normybyy^2);
    bcn = beta3*constant_c0*normMw^2;
    bbb = btnn +  bcn + btkkc; 
     end
end
