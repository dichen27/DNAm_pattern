function sample_matrix=cd_sample_matrix(a,b,c)


sample_matrix=zeros(a,b,c);
j=1;
for fl=1:c
A=[];
for a=1:a
    A=[A;(j:j+(b-1))];
    j=j+b;  
end
sample_matrix(:,:,fl)=A;
end


end