% 结构图
clc,clear
N=[0 0 0;1 -1 0;2 -1.5 0;3 -1 0;4 0 0;3 1 0;2 1.5 0;1 1 0]';N=N';
C_b_in = [1 2;2 3;3 4;4 5;2 8;3 7;4 6;1 8;8 7;7 6;6 5]; 
for i=1:11
    for j=1:2
        var1(:,j)=C_b_in(i,j);
    end
    var2=var1(1);var3=var1(2);
    var4=N(var2,:);var5=N(var3,:);
    var6=var4(1);var7=var4(2);var8=var5(1);var9=var5(2);
    x=[var6,var8];y=[var7,var9];
    plot(x,y,'k','LineWidth',2)
    hold on;
end
% 轴力图
1
