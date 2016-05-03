tic; 
%two d ising model with corner transfer matrix rg
clear;
%rng(2); % random number initialization

%boundary tensors
chi=20; %bond dimension of the boundary tensors
nco=4; % number of corners
nt=4; % number of edges
C=cell(1,nco); % corner tensors
nte=2;% number ot t tensor per edge
T=cell(nt,nte); % T tensors (edge tensors)

%bulk tensors
D=2; % bond dimension of the bulk tensors
na=2; %number of bulk tensor in the unit cell
A=cell(1,na);
itmax=1000; %maximum number of iterations in the CTMRG method
eps=0.00000000001; % criterion to stop iteration of the cmtrg
beta=log(1+sqrt(2))/2+0.1; % critical beta=log(1+sqrt(2))/2
beta=log(1+sqrt(2))/2-0.005;

format long
 % initializing tensors
for i=1:nco
  C{i}=1*ones(chi)+0.0000*rand(chi);
end
    
for i=1:nt
  for j=1:nte
     T{i,j}=1*ones(chi,chi,D)+0.0000000*rand(chi,chi,D);
   end
end


sh=cell(1,2);

Tr=2.0
each=100;

for Temp=Tr:Tr %2:0.01:2.4 %2.0:0.01:2.4  %2.26:0.001:2.28
    
    beta=1/Temp % inverse T
    [A,B,T,C]=buildA(A,beta,T,C); % building the bulk tensors and magnetization tensor
    
    %for i=1:2000
    delta=1;
    magold=10000;
    while(delta>eps)   
                      
        [C,T]=ctmrgloop(C,T,A,D,chi,i,each); % corner transfer main loop
        [mag,psi]=observ(A,B,C,T,D,chi); % magnetization measurement
         mag;
         delta=abs(mag-magold);
         magold=mag;
         
         M=(1-[(sinh(2.0*beta))^2 ]^-2 )^(1/8);
    end
    
                     
     [mag,psi]=observ(A,B,C,T,D,chi);
     [freee,e,beta_c]=ising2d_square_exact_results(beta);
     M=(1-[(sinh(2.0*beta))^2 ]^-2 )^(1/8);
     format long
     mag*2;
    
   %fprintf(fileID,'%16.12f %16.12f %16.12f \n',Temp,abs(mag),e);
   

    
end

%fclose(fileID);
mag % numerical result from CTMRG
M   % analytical  

TimeSpent = toc





