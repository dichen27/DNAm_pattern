function p_3=cd_FDR(p_all)





p_1=reshape(p_all,size(p_all,1)*size(p_all,2),1);


p_2=mafdr(p_1,'BHFDR', true);

p_3=reshape(p_2,size(p_all,1),size(p_all,2));


end