function [y] = haar_column_vector(k, t)
A=zeros(k,1);
for i=1:k
A(i)=haar_column_element(i-1,t);
end
 y=A;
end

