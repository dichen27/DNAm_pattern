
function [] = cd_circos(t_circos,p_circos,threshold)



net_data_aal2_brodmann_all = double(p_circos<threshold & p_circos>0);

net_data_aal2_brodmann_increase= t_circos>0 & net_data_aal2_brodmann_all>0;
net_data_aal2_brodmann_decrease = t_circos<0 & net_data_aal2_brodmann_all>0;
[sum(net_data_aal2_brodmann_decrease(:))/2 sum(net_data_aal2_brodmann_increase(:))/2];

edge_circos = -net_data_aal2_brodmann_decrease + net_data_aal2_brodmann_increase;

measure = [sum(abs(net_data_aal2_brodmann_increase),2), sum(abs(net_data_aal2_brodmann_decrease),2)];

linkfile = edge_circos;


parse_wei(measure, linkfile)


system('cd C:\softc\mat_tool\circos_wei')

system('perl C:\softc\mat_tool\circos_wei\parsemap -map matmap.txt -links links.txt -confdir C:\softc\mat_tool\circos_wei\etc -datadir C:\softc\mat_tool\circos_wei\data')
% system(['D:matlab_toolbox/circos-0.69-6/bin -conf C:/Users/jingnandu/Desktop/Drinking_and_Smoking/etc/circos.conf ' ...
%     '-outputdir C:/Users/jingnandu/Desktop/Drinking_and_Smoking/output'])

system('C:\softc\circos-0.69-6\bin\circos -conf C:\softc\mat_tool\circos_wei\etc\circos.conf -outputdir C:\Users\dell\Desktop')
end