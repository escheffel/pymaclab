function [azero]=mksimzha(x);
x=x';
azero=[x(1:7); 0  x(8:9) -x(8) 0 -x(8) 0; x(10:12)zeros(1,4); x(13) zeros(1,2) x(14:17);x(18) zeros(1,3) x(19:21) ;x(22) zeros(1,4) x(23:24);x(25) zeros(1,5) x(26)];
