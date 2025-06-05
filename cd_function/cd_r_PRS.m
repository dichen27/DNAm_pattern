
  function r_PRS=cd_r_PRS(PRS_ID_intersect,brain_intersect,cov_intersect)





r_PRS=[];
for i=1:size(brain_intersect,2)
    data=[brain_intersect(:,i),PRS_ID_intersect,cov_intersect];    
    data(any(isnan(data), 2),:) = [];% dele any row NaN

[r,p]=partialcorr(data(:,1),data(:,2),data(:,3:end));

r_PRS=[r_PRS;r];
end




end