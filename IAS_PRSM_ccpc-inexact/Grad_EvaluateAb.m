function g = Grad_EvaluateAb(Ab, x)
 [n,~]  = size(Ab); 
  
 temp = exp(-Ab*x);
 g = -1/n*Ab'*(1./(1.+1./temp)) ; 
 
end

