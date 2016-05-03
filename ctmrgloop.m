function [ C,T ] = ctmrgloop( C,T,A,D,chi,ii,each )
%cmtrgloop updates tensors C and T using the ctmrg procedure
%

 %disp('---------------------')
sh=cell(1,8);
% LEFT MOVE
% 
C1T=scon({C{1},T{1,1}},{[-1 1],[1 -3 -2]});
C1T=reshape(C1T,chi*D,chi);
C4T=scon({C{4},T{3,2}},{[1 -2],[-1 1 -3 ]});
C4T=reshape(C4T,chi,chi*D);

Q1T=scon( {C{1},T{1,1}, T{4,2},A{1} },{[4 1 ],[1 -4 2 ],[-1 4 3],[2 -3 -2 3]});
Q4T=scon( {C{4},T{4,1}, T{3,2},A{2}},{[4 1 ],[1 -3 2 ],[-1 4 3],[-4 -2 3 2]});
Q1T=reshape(Q1T,chi*D,D*chi);
Q4T=reshape(Q4T,chi*D,D*chi);

T42T=reshape(scon({T{4,2},A{1}},{[-1 -4 1],[-5 -3 -2 1]}),chi*D,D,chi*D);
T41T=reshape(scon({T{4,1},A{2}},{[-1 -4 1],[-5 -3 -2 1]}),chi*D,D,chi*D);
%T4T=reshape(scon({T{4},A{1}},{[-1 -4 1],[-5 -3 -2 1]}),chi*D,D,chi*D);

mat=Q1T*Q1T'+Q4T'*Q4T;
[W,DD] = eigs(mat,chi);
[DD, order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
W = W(:,order);


mat=C1T*C1T'+C4T'*C4T;
[V,DD] = eigs(mat,chi);
[DD, order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
V = V(:,order);


%T{4,2}=scon({T42T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{4,1}=scon({T41T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

%T{4,2}=scon({T42T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{4,1}=scon({T41T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

% T{4,1}=scon({T42T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
% T{4,2}=scon({T41T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});


T{4,1}=scon({T42T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
T{4,2}=scon({T41T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

%T{4}=scon({T4T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

C{1}=V'*C1T;
C{4}=C4T*V;

% % %Right MOVE

C3T=scon({C{3},T{3,1}},{[-1 1],[1 -3 -2]});
C3T=reshape(C3T,chi*D,chi);
C2T=scon({C{2},T{1,2}},{[1 -2],[-1 1 -3 ]});
C2T=reshape(C2T,chi,chi*D);

Q3T=scon( {C{3},T{3,1}, T{2,2},A{1}},{[4 1 ],[1 -4 2 ],[-1 4 3],[-2 3 2 -3]} );
Q2T=scon( {C{2},T{2,1}, T{1,2},A{2}},{[4 1 ],[1 -3 2 ],[-1 4 3],[3 2 -4 -2]});
Q3T=reshape(Q3T,chi*D,D*chi);
Q2T=reshape(Q2T,chi*D,D*chi);

T22T=reshape(scon({T{2,2},A{1}},{[-1 -4 1],[-2 1 -5 -3]}),chi*D,D,chi*D);
T21T=reshape(scon({T{2,1},A{2}},{[-1 -4 1],[-2 1 -5 -3]}),chi*D,D,chi*D);

mat=Q3T*Q3T'+Q2T'*Q2T;
[W,DD] = eigs(mat,chi);
[DD, order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
W = W(:,order);


mat=C3T*C3T'+C2T'*C2T;
[V,DD] = eigs(mat,chi);
[DD ,order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
V = V(:,order);

%T{2}=scon({T2T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});


%T{2,2}=scon({T22T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{2,1}=scon({T21T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

%T{2,2}=scon({T22T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{2,1}=scon({T21T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});


% T{2,1}=scon({T22T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
% T{2,2}=scon({T21T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});


T{2,1}=scon({T22T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
T{2,2}=scon({T21T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

C{3}=V'*C3T;

C{2}=C2T*V;
% 
% 
% % 
% % % 
% % %UP move
% % 
C2T=scon({C{2},T{2,1}},{[-1 1],[1 -3 -2]});
C2T=reshape(C2T,chi*D,chi);
C1T=scon({C{1},T{4,2}},{[1 -2],[-1 1 -3 ]});
C1T=reshape(C1T,chi,chi*D);

Q2T=scon( {C{2},T{2,1}, T{1,2},A{2} },{[4 1 ],[1 -4 2 ],[-1 4 3],[3 2 -3 -2]} );
Q1T=scon( {C{1},T{1,1}, T{4,2},A{1}}, {[4 1 ],[1 -3 2 ],[-1 4 3],[2 -4 -2 3]});
Q2T=reshape(Q2T,chi*D,D*chi);
Q1T=reshape(Q1T,chi*D,D*chi);


T12T=reshape(scon({T{1,2},A{2}},{[-1 -4 1],[1 -5 -3 -2 ]}),chi*D,D,chi*D);
T11T=reshape(scon({T{1,1},A{1}},{[-1 -4 1],[1 -5 -3 -2 ]}),chi*D,D,chi*D);

mat=Q2T*Q2T'+Q1T'*Q1T;
[W,DD] = eigs(mat,chi);
[DD, order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
W = W(:,order);

mat=C2T*C2T'+C1T'*C1T;
[V,DD] = eigs(mat,chi);
[DD, order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
V = V(:,order);

%T{1,2}=scon({T12T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{1,1}=scon({T11T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

%T{1,2}=scon({T12T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{1,1}=scon({T11T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
% 
% T{1,1}=scon({T12T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
% T{1,2}=scon({T11T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});


T{1,1}=scon({T12T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
T{1,2}=scon({T11T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

C{2}=V'*C2T;
C{1}=C1T*V;
% 
% % 
% % 
% % % 
% % % 
% % % 
% % % %down move
% 
C4T=scon({C{4},T{4,1}},{[-1 1],[1 -3 -2]});
C4T=reshape(C4T,chi*D,chi);
C3T=scon({C{3},T{2,2}},{[1 -2],[-1 1 -3 ]});
C3T=reshape(C3T,chi,chi*D);

Q4T=scon( {C{4},T{4,1}, T{3,2},A{2}},{[4 1 ],[1 -4 2 ],[-1 4 3],[-3 -2 3 2]});
Q3T=scon( {C{3},T{3,1}, T{2,2},A{1}},{[4 1 ],[1 -3 2 ],[-1 4 3],[-2 3 2 -4]});
Q4T=reshape(Q4T,chi*D,D*chi);
Q3T=reshape(Q3T,chi*D,D*chi);

T32T=reshape(scon({T{3,2},A{2}},{[-1 -4 1],[-3 -2 1 -5]}),chi*D,D,chi*D);
T31T=reshape(scon({T{3,1},A{1}},{[-1 -4 1],[-3 -2 1 -5]}),chi*D,D,chi*D);

mat=Q4T*Q4T'+Q3T'*Q3T;
[W,DD] = eigs(mat,chi);
[DD order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
W = W(:,order);


mat=C4T*C4T'+C3T'*C3T;
[V,DD] = eigs(mat,chi);
[DD order] = sort(diag(DD),'descend');  %# sort eigenvalues in descending order
V = V(:,order);


%T{3}=scon({T3T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{3,2}=scon({T32T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{3,1}=scon({T31T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

%T{3,2}=scon({T32T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
%T{3,1}=scon({T31T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
% 
% T{3,1}=scon({T32T,V,W'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
% T{3,2}=scon({T31T,W,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});


T{3,1}=scon({T32T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});
T{3,2}=scon({T31T,V,V'},{[1  -3 2 ],[ 2 -2],[-1 1 ]});

C{4}=V'*C4T;
C{3}=C3T*V;


% 
% % disp(sh);
% 
% %C{1}
% %T{3}


for j=1:4
   for i=1:2
   T{j,i}=T{j,i}/max(abs(T{j,i}(:)));
   end 
   C{j}=C{j}/max(abs(C{j}(:)));
   %C{j}={j}/max(abs(T{j}(:)));
  
end

 
 end
% 
