function palpha = operational_matrix( n , alpha ,bound, H)
%returns the operational matrix
% e.g.  here n is k in Chen's paper e.g. n=256,  
% bound is b in the paper A. Kilicman, Z.A.A. Al Zhour / Applied Mathematics and Computation 187 (2007) 250–265 


F=zeros(n,n);

for i=1:1:n
    for j=1:1:n
        if ((j-i)>0)
            F(i,j)=element_of_matrix_F(j-i,alpha);
            
        elseif(j==i)
            F(i,j)=1;
            
        else
            F(i,j)=0;
        end
        
        F(i,j)=F(i,j)/(gamma(alpha+2)*(n^alpha));
    end
end

palpha = (bound^alpha) * H * F/H ;
% using palpha=H*F_alpha*inverse(H)
% F_alpha defined on the interval [0,b] i.e. [0,bound] is equqtion (37) in
% the paper " A. Kilicman, Z.A.A. Al Zhour / Applied Mathematics and Computation 187 (2007) 250–265"

end

% A. Kilicman, Z.A.A. Al Zhour / Applied Mathematics and Computation 187 (2007) 250–265
% Y. Chen et al. / Journal of Computational Science 3 (2012) 367–373
