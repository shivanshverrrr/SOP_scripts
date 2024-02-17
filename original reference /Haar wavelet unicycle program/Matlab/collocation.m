function f=collocation(k)
z=zeros(1,k)
for i=1:k
    z(i)=(i-(1/2))/k;
end
f=z;
end