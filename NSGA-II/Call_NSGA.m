%% ZDT4-R
% clear all
% for rep=1:1      
%     f='ZDT4-R';
%     L=-5*ones(1,10);
%     U=5*ones(1,10);
%     L(1)=0;U(1)=1;
%     pop=100;
%     no_of_objs = 2;
%     store_ZDT4_R(rep).data=NSGA_II(f,no_of_objs,L,U,pop);
% end
% save('store_ZDT4_R','store_ZDT4_R')

%% ZDT4-RC
clear all
for rep=1:1
    f='ZDT4-RC';
    L=-5*ones(1,20);
    U=5*ones(1,20);
    L(1)=0;U(1)=1;
    pop=50;
    no_of_objs = 2;
    trial(rep).data=NSGA_II(f,no_of_objs,L,U,pop);
end
save('trial','trial')

%% ZDT4-A
% clear all
% for rep=1:1   
%     f='ZDT4-A';
%     L=-32*ones(1,10);
%     U=32*ones(1,10);
%     L(1)=0;U(1)=1;
%     pop=100;
%     no_of_objs = 2;
%     store_ZDT4_A(rep).data=NSGA_II(f,no_of_objs,L,U,pop);
% end
% save('store_ZDT4_A','store_ZDT4_A')
% 
%% ZDT4-G
% clear all
% for rep=1:31    
%     f='ZDT4-G';
%     L=-512*ones(1,10);
%     U=512*ones(1,10);
%     L(1)=0;U(1)=1;
%     pop=100;
%     no_of_objs = 2;
%     store_ZDT4_G(rep).data=NSGA_II(f,no_of_objs,L,U,pop);
% end
% save('store_ZDT4_G','store_ZDT4_G')

%% ZDT4-S
% clear all
% for rep=1:1    
%     f='ZDT4-S';
%     L=-500*ones(1,10);
%     U=500*ones(1,10);
%     L(1)=0;U(1)=1;
%     pop=100;
%     no_of_objs = 2;
%     store_ZDT4_S(rep).data=NSGA_II(f,no_of_objs,L,U,pop);
% end
% save('store_ZDT4_S','store_ZDT4_S')

%% ZDT1
% clear all
% for rep=1:11    
%     f='ZDT1';
%     L=0*ones(1,10);
%     U=1*ones(1,10);
%     L(1)=0;U(1)=1;
%     pop=100;
%     no_of_objs = 2;
%     store_ZDT1(rep).data=NSGA_II(f,no_of_objs,L,U,pop);
% end
% save('store_ZDT1','store_ZDT1')

%% ZDT3
% clear all
% for rep=1:1   
%     f='ZDT3';
%     L=-0*ones(1,30);
%     U=1*ones(1,30);
%     L(1)=0;U(1)=1;
%     pop=100;
%     no_of_objs = 2;
%     store_ZDT3(rep).data=NSGA_II(f,no_of_objs,L,U,pop);
% end