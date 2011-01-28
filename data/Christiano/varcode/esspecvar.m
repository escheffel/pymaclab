function hspec = esspecvar(ndata,nlags,lamj);
%Estimate the spectrum of a multivariate series using a Vector Autoregression

[beta,sig]=estimatevar(ndata,nlags);

AL = zeros(2,2,nlags+1);
nvars = size(ndata,2);

AL(:,:,1) = eye(nvars);

beg=0;

%In this loop we create the polynomial matrix A(L)such that
% A(L)y = u where y is the data and u is the error term.

for zv = 1:nlags;
   enz=beg+nvars;
   AL(:,:,zv+1)=-1*beta(:,beg+1:enz);
   beg=enz;
end;


%In this loop we create g(omg) = A(exp(-i*omg)). Note that we have
%to flip around the polynomial term so that it reads correctly.
% In other words we flip it so that the first term is the coefficient on L^nlags

for za=1:nvars;
   for zb=1:nvars;
      tmpz = squeeze(AL(za,zb,:));
      tmpz = tmpz(end:-1:1);
      g(za,zb,:)=polyval(tmpz,exp(-i*lamj));
   end;
end;

%for this loop we construct the spectrum as the product
% of inv(g)*sigma*inv(g)' 
% the inv(g) corresponds to the moving average representation of the data

for zm=1:length(lamj);
   ginv = inv(g(:,:,zm));
   hspec(:,:,zm)=(1/(2*pi))*ginv*sig*ginv';
end;
