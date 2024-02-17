% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %  for a given value of n, the value of j and k can be determined with this program
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [y1,y2] = haar_j_and_k_values_for_given_n(n)
if n == 0
    j=0;
elseif n ~= 0   
    j = 0;
   while( (2^j) < ( n + 1))     
        p=2^j;
        j = j + 1;
   end
    j=log2(p);
    k = (n-p);
end
y1=j;
y2=k;
end
%  [jvalue,kvalue] = j_and_k_values_for_given_n(n)