function  d = cd_t2d(t,n1,n2,num_cov)  
% source
% https://www.polyu.edu.hk/mm/effectsizefaqs/effect_size_equations2.html

% t=3.6
% n1=423
% n2=72
% num_cov=22

% t=kan_tval_b_nor
% length(Z_3_pat_b),
% length(Z_3_nor),
% size(cov_b_nor,2)


    p=n1/(n1+n2);
    q=n2/(n1+n2);
    xi=sqrt(1/(p*q));
    df=(n1+n2)-2-num_cov;
    d=zeros(length(t),1);
for i=1:length(t)
    

 d(i)=xi*t(i)/sqrt(df);
end
    
   
% d=(t*(n1+n2))/sqrt((n1+n2-2)*(n1*n2));

end


