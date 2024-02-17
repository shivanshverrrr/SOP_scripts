function y=haar_column_element(n,t)
 scaling_function=@(x) 0 + 1.*(x>=0&x<1) ;
 Mother_wavelet=@(x) 0 + 1.*(x>=0&x<0.5) + -1.*(x>=0.5&x<1);
 if(n==0)
y= scaling_function(t);
 elseif n>0
     [j,k] = haar_j_and_k_values_for_given_n(n);
 y=Mother_wavelet( ((2^(j))*t)  - k );
end