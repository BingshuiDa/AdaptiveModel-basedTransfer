function funcval=LS(y,fhandles,L,U,no_of_objs)
    var = L + y.*(U-L);
    if no_of_objs == 2
        w1 = rand(1);
        w2 = (1-w1);
        f1=fhandles.fhandle1(var);
        f2=fhandles.fhandle2(var);
%         funcval = max([w1*f1,w2*f2]);
        funcval = w1*f1 + w2*f2;
    else
        w_1 = rand(1);
        w_2 = rand(1);
        w_3 = rand(1);
        w1 = w_1/(w_1 + w_2 + w_3);
        w2 = w_2/(w_1 + w_2 + w_3);
        w3 = w_3/(w_1 + w_2 + w_3);
        f1=fhandles.fhandle1(var);
        f2=fhandles.fhandle2(var);
        f3=fhandles.fhandle3(var);
%         funcval=max([w1*f1,w2*f2,w3*f3]);
        funcval = w1*f1 + w2*f2 + w3*f3;
    end
end