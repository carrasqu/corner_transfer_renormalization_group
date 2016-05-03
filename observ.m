
function [ mag,psi ] = observ( A,B,C,T,D,chi )

%observ computes the magnetization of the ising model in 2 d
%   Detailed explanation goes here

% 
% T{2,1}=T{1,1};
% T{4,1}=T{1,1};
% 
% T{2,2}=T{1,2};
% T{4,2}=T{1,2};
% 
% T{2,3}=T{1,3};
% T{4,3}=T{1,3};
% 
% C{2}=C{1};
% C{4}=C{1};

psi=scon({C{1},C{2},C{3},C{4},...
          T{1,1},T{1,2},...
          T{2,1},T{2,2},...
          T{3,1},T{3,2},...
          T{4,1},T{4,2},...
          A{2},A{1},A{2}},...
          {[4 1],[3 6  ],[17 20 ],[18 14],...
          [1 2 -1],[2 3 5],...
          [6 10 7],[10 17 13],...
          [20 19 16],[19 18 15],...
          [14 8 11],[8 4 -4],...
          [5 7 9 -2],[9 13 16 12],[-3 12 15 11]});


ma=scon({psi,B{1}},{[1 2 3 4],[ 1 2 3 4]});
zeta=scon({psi,A{1}},{[1 2 3 4],[1 2 3 4]}) ;
mag=ma/zeta;



end

