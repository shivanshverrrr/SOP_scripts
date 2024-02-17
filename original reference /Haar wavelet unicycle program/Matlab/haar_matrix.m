function y = haar_matrix(k)
t=collocation(k);
B=zeros(k);
    for i=1:k
B(:,i)=haar_column_vector(k, t(i));
    end
y=B;
end
