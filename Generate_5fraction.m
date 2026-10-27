clc;

clear;

interval=1;
n=100/interval;

a=0;
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                for m=1:n
                    if i+j+k+l+m==100/interval
                        %if ((i*interval>=52&&i*interval<=56)&&(j*interval>=14&&j*interval<=18))&&((k*interval>=12&&k*interval<=16)&&(l*interval>=8&&l*interval<=12))&&m*interval>=4&&m*interval<=8
                        a=a+1;
                        A(a)=i*interval;
                        B(a)=j*interval;
                        C(a)=k*interval;
                        D(a)=l*interval;
                        E(a)=m*interval;
                        %end
                    end
                end
            end
        end
    end
end

frac=[A;B;C;D;E];

N=length(frac);

% b=0;
% for i=1:N
%     if max(frac(:,i))==frac(1,i)
%         b=b+1;
%         fraction(:,b)=frac(:,i);
%     end
% end
id=frac(1,:)<30 | frac(1,:)>55 | frac(2,:)<10 | frac(2,:)>20|frac(3,:)<10|frac(3,:)>20|frac(4,:)<10|frac(4,:)>20|frac(5,:)<10|frac(5,:)>30;
frac(:,id)=[];
fraction=frac';