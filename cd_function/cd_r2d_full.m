function  d = cd_r2d_full(r,n1,n2)  
% source
% https://www.polyu.edu.hk/mm/effectsizefaqs/effect_size_equations2.html
    p=n1/(n1+n2);
    q=n2/(n1+n2);
    xi=sqrt(1/(p*q));

d=zeros(length(r),1);
for i=1:length(r)
    
d(i)=(xi*r(i))/sqrt(1-r(i)^2);
end

end
