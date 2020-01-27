function y=xpolate(x,m,n)

% xpolate(x,m,n)
%             extra/interpolate using tony's 2-D fft method
%             replaces nans in a 2-D array with sensible values
%
%       x     is the 2-D array to be filled
%       m,n   optional arguments to specify size of 2-D fft
%             if n is not specified, n=m;
%             if m and n are not specified, [m n]=size(x);
%
%             the algorithm truncates progressively fewer
%             high frequency components in order to achieve
%             smooth transitions to areas with no data

[M N]=size(x);

if     nargin==2, n=m;
elseif nargin==1, m=M; n=N;
end

y=nan*ones(m,n); y(1:M,1:N)=x;

[X,Y]=meshgrid(1:n,1:m); X=X-1; Y=Y-1; % distance from top left corner
ix=find(isnan(y));                   % location of nans
y(ix)=mmean(x)*ones(size(ix));       % replace nans by mean value
mc=fix(m/2)+1; nc=fix(n/2)+1;        % center value
w=zeros(m,n);
maxx=mmax(x); minx=mmin(x);

for i=1:(max(mc,nc)-1)
  
  % take 2-d fft of y
  
  f=fftshift(fft2(y));
  
  % truncate high frequency

  mx=max((mc-i),1):min((mc+i),m);
  nx=max((nc-i),1):min((nc+i),n);
  w(mx,nx)=ones(length(mx),length(nx));
  f=f.*w;
  
  % map back onto y
  
  f=real(ifft2(ifftshift(f)));       % revert to space domain
  y(ix)=f(ix);                       % replace missing values
  y(find(y>maxx))=maxx;
  y(find(y<minx))=minx;
  
end

y=y(1:M,1:N);
