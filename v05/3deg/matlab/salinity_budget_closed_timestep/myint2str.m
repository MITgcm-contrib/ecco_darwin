function s=myint2str(x,n);
%MYIN2STR(X,N)   convert integer to string with leading zeros
%                padds integer X with zeros to N locations

if nargin < 2, n=2;  end
if nargin < 1, help myint2str, return, end

%s=int2str(x);
%n=n-length(s);
%for i=1:n
%  s=['0' s];
%end

% replacement code contributed by M. Losch, November 28, 2008
 fmt = sprintf('%%0%iu',n);
 s=sprintf(fmt,x);
