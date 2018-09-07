data = [];
x=0;
for i = 1:51
    y = 1 - sqrt(x);
    data = [data;x y];
    x=x+0.02;
end
%  plot(data(:,1),data(:,2),'o')  
save('data','data');