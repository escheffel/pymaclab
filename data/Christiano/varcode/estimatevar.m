function [beta,sigma,uz,vcovbeta]=estimatevar(data,nlags,hasconstant);
% ESTIMATEVAR Estimates a Vector Autoregression. 
%	[Beta,Sigma,UZ]=ESTIMATEVAR(DATA,nlags,hasconstant);
%	Finds coefficients estimates BETA, variance covariance matrix SIGMA, and residuals UZ
%	that are the results of estimate by OLS a VAR on the matrix DATA.
% 	The matrix DATA  has each row being an observations
% from all the variables. NLAGS is the number of lags.  
% HASCONSTANT is equal to one if for a VAR with a constant and equal to zero otherwise. 


if nargin == 2;
   hasconstant=0; %the default value for hasconstant
end;



% A bunch of error checks
if and(nlags==0,hasconstant==0);
   error('WARNING: you are trying to estimate a VAR without a constant or lags');
 end;

if nlags>size(data,1);
   error('WARNING: you must have more observations than lags')   
end;


if size(data,1)<size(data,2);
   error('WARNING: You have more variables than observations or you have entered the data incorrectly.')
end;

nvars=size(data,2); % nvars equals number of columns

ibeg = nlags+1;  %first observation used as dependent variable
iend = size(data,1); %use all the data

% dependent variable
z=data(ibeg:iend,:);

% set up independent matrix 
% First see if we have any lags
if nlags > 0 
   xz = data(ibeg-1:iend-1,:);
   % adding on lags 
   if nlags>1;  
      for zlag=2:nlags;
	   	% if nlags equals one this step is skipped
   		xz=[xz data(ibeg-zlag:iend-zlag,:)];
		end;
	end;
else
   %  we just set xz to be an empty matrix.
   xz = [];
end;


if hasconstant 
   % we add in a constant term
   xz=[ones(iend-ibeg+1,1) xz];
end;



% I calculate the VAR using the matlab \ operator
beta=(xz'*xz)\(xz'*z);
uz = z-xz*beta;
sigma=uz'*uz/length(uz);
% for ease of use I return the transpose of beta.
beta = beta';


%variance covariance matrix for the parameters
%vcovbeta=kron(uz'*uz/(length(uz)-size(beta,2)),inv(xz'*xz));
vcovbeta=kron(uz'*uz/(length(uz)),inv(xz'*xz));


% this is equivalent to

% % now we declare the coefficient matrix
% betaz=zeros(nlags*nvars+hasconstant,nvars);
% % make the projection matrix
% ixzxz=inv(xz'*xz)*xz';
% for zv=1:nvars;
%   betaz(:,zv)=ixzxz*z(:,zv);
%   uz(:,zv)=[(z(:,zv)-xz*betaz(:,zv))];
% end;
% sigma=uz'*uz/length(uz);


% now to make the standard errors. 


