function [mse,mseconts]=vardecompnew(betaz,a0,nlags,errshk,steps)
%VARDECOMP calculates the how much the errshk error contributes
%to the 1 through steps period ahead mean squared forecast error
%betaz coefficients a0 variance covariance decomposition.
%output mse total mean squared error mseconsts contribution to mse by errshk
% 

fctz = inv(a0);
sigma = fctz*fctz';
nvars = length(sigma);
hascon = size(betaz,1)-nvars*nlags;

mseconts = zeros(nvars,nvars,steps);
mse = zeros(nvars,nvars,steps);
psiz = zeros(nvars,nvars,steps);

 %Fastest way? Maybe Not.
 %However this takes advantage of the relationship
 % between the impulse responses and the 
 % moving average representations

for zx=1:nvars;
   impz=mkimprep(betaz,eye(7),nlags,zx,steps);
	psiz(:,zx,:)=impz';   
end;

%Calculate Total Mean Squared Error
mse(:,:,1)=sigma;

for kr=2:steps;
   mse(:,:,kr)=mse(:,:,kr-1)+psiz(:,:,kr)*sigma*psiz(:,:,kr)';
end;

%Calculate the errshk contribution
zk = errshk;


aj=fctz(:,zk);
mseconts(:,:,1)=aj*aj';

for kr=2:steps;
   mseconts(:,:,kr)=mseconts(:,:,kr-1)+psiz(:,:,kr)*aj*aj'*psiz(:,:,kr)';   
end;

