function bootstat_adjust=cd_bootstat_adjust(bootstat,xi)


bootstat_cov=bootstat;
% xi=std_target/std(bootstat);

bootstat_adjust=[];
for i=1:length(bootstat_cov)
    
%     lin=(bootstat_cov(i)-mean(bootstat_cov))*xi;
    lin=(bootstat_cov(i)-mean(bootstat_cov))*xi+mean(bootstat_cov);
    
    bootstat_adjust=[bootstat_adjust;lin];
end


end