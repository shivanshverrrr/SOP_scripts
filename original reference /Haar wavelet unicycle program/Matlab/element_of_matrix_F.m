function value = element_of_matrix_F( k,alpha )
%%return elements of F which is used to calculate the operational matrix
% A. Kilicman, Z.A.A. Al Zhour / Applied Mathematics and Computation 187 (2007) 250–265
% Y. Chen et al. / Journal of Computational Science 3 (2012) 367–373

    value = ((k+1)^(alpha+1)) - 2*((k)^(alpha+1)) + ((k-1)^(alpha+1));
end

