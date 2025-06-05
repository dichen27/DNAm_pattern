function  d = cd_r2d(r)  
% source
% https://www.polyu.edu.hk/mm/effectsizefaqs/effect_size_equations2.html

d=zeros(length(r),1);
for i=1:length(r)
    
d(i)=(2*r(i))/sqrt(1-r(i)^2);
end

end

