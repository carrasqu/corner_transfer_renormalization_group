function [ A,B,T,C ] = buildA( A,beta,T,C )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


%mmax=1; % kagome
mmax=2; % square ising model with 4 sites


for m=1:mmax
    
    for i1=1:2
        for i2=1:2
            for i3=1:2
                for i4=1:2
                    
                    A{m}(i1,i2,i3,i4)=exp(beta*( (2*i1-3)*(2*i2-3) + (2*i2-3)*(2*i3-3)+ (2*i3-3)*(2*i4-3)+(2*i4-3)*(2*i1-3)));
                   
                end
            end
        end
    end
    
end

z=[1.0 0; 0 -1.0];


%B{1}=scon({A{1},z,z},{[-1 -2 2 1],[-4 1],[-3 2]});

B{1}=scon({A{1},z},{[-1 -2 -3 1],[-4 1]});


xx=(max(A{1}(:)));

B{1}=B{1}/xx;

for m=1:mmax

    A{m}=A{m}/xx;
end



% corners


% for m=1:4
%     for i1=1:2
%         for i2=1:2
%             
%             C{m}(i1,i2)=exp(beta*((2*i1-3)*(2*i2-3)));
%             
%         end
%     end
%     
% end
% 
% % edges
% 
% for m=1:4
%     %for kk=1:3
%         
%         for i1=1:2
%             for i2=1:2
%                 for i3=1:2
%                     
%                    T{m}(i1,i2,i3)=exp(beta*((2*i1-3)*(2*i2-3) + (2*i1-3)*(2*i3-3)));
%                 end
%             end
%         end
%     %end
%     
% end
% 

% 
for i1=1:2
        for i2=1:2

aa(i1,i2)=exp(beta*( (2*i1-3)*(2*i2-3) ) );
        end 
end

[V,D]=eig(aa);

V1=V*sqrtm(D);
V2=sqrtm(D)*V';

D

    for i1=1:2
        for i2=1:2
            for i3=1:2
                for i4=1:2
                           A{1}(i1,i2,i3,i4)=V1(1,i1)*V1(1,i2)*V2(i3,1)*V2(i4,1)+V1(2,i1)*V1(2,i2)*V2(i3,2)*V2(i4,2);
                             B{1}(i1,i2,i3,i4)=V1(1,i1)*V1(1,i2)*V2(i3,1)*V2(i4,1)-V1(2,i1)*V1(2,i2)*V2(i3,2)*V2(i4,2);
                end 
            end 
        end
    end
    
A{2}=A{1};



end

