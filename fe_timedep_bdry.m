function udot = fe_timedep_bdry(t, u)
% usage: udot = fe_timedep_bdry(t, u)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved

% extract problem data
global Pdata;
c = Pdata.c;
m = Pdata.m;
xspan = Pdata.xspan;
b1dot = Pdata.b1dot;
b2dot = Pdata.b2dot;

% initialize RHS terms
udot = zeros(m,1);

% interior index set
int = 2:m-1;
x = xspan(int);

% forcing term
udot(int) = c(1)*cos(x*1) + c(2)*cos(x*2) + c(3)*cos(x*3) + c(4)*cos(x*4) + c(5)*cos(x*5);

% left boundary
udot(1) = b1dot(t);

% right boundary
udot(m) = b2dot(t);

% end function
