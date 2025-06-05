function  accuracy = cd_accuracy( predict, true )  

if length(predict)~=length(true)
    error('cd27£ºnumber predict not equal to true !');   
end

num_rate=0;
for i=1:length(predict)
    if predict(i)==true(i)
num_rate=num_rate+1;
    end
end


accuracy=num_rate/length(predict);
end