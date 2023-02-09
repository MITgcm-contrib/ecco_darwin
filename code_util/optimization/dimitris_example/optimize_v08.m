% Read-in observations and model-equivalent structures.
clear all;
close all

cd('~/Google Drive/ECCO-Darwin/mat/')
load V4_model_obs

% Choose data constraints to include.

Type=fieldnames(observations);

Type=Type([1 3:8]);
n_rows=0;

for n=1:length(Type)
    eval(['n_rows=n_rows+length(observations.' Type{n} ');'])
end

% Cell array data_id identifies each column of matrix "data".

data_id{1}='observations';

tmp=fieldnames(model);

for n=1:length(tmp)
    data_id{n+1}=tmp{n};
end

n_columns=length(data_id);

% Construct matrix "data" of observations and sensitivity experiments.
% Structure data_info contains meta information about each row in "data".

data=zeros(n_rows,n_columns);
data_info.Type=zeros(n_rows,1);
column=1;
row=0;

for n=1:length(Type)
    eval(['tmp=observations.' Type{n} ';'])
    ir=1:length(tmp);
    data(row+ir,column)=tmp;
    data_info.Type(row+ir,column)=n;
    row=row+length(ir);
end

for column=2:n_columns
    row=0;
    for n=1:length(Type)
        eval(['tmp=model.' data_id{column} '.' Type{n} ';'])
        ir=1:length(tmp);
        data(row+ir,column)=tmp;
        row=row+length(ir);
    end
end

% Remove any rows that contain negative values.
ir=find(any(data'<=0));
data(ir,:)=[];
tmp=fieldnames(data_info);

for n=1:length(tmp)

    eval(['data_info.' tmp{n} '(ir)=[];'])

end

clear c* ir m* n* o* r* t*

%% 

% Define observations, baseline, and perturbation columns in "data".
d = 1;             % observations
b = 6;             % baseline
p = [3 2; 4 2; 5 2]; % perturbations pairs

% Build model-data difference vector, Eq. 6 in
% Menemenlis, Fukumori, and Lee, (2005), hereinafter MFT05.
y=data(:,d)-data(:,b);

% Build kernel matrix, Eq. 7 in MFT05.
G=[data(:,p(:,1))-data(:,p(:,2))];

% Create estimate for diagonal elements of error
% covariance matrix R in Eq. 5 of MFT05.
Robs=0*y;  % variance of observations
Rbase=0*y; % variance of baseline
Rdiff=0*y; % variance of observations-baseline difference

for n=1:length(Type)

    in=find(data_info.Type==n);
    Robs(in)=var(data(in,d));
    Rbase(in)=var(data(in,b));
    Rdiff(in)=var(data(in,b)-data(in,d));

end

R=max([Robs Rbase Rdiff]')';
clear n in Robs Rbase Rdiff

%% 

% Remove >3sigma outliers.

Mdata=0*data; % mean values for each data Type
Mdiff=0*data; % mean model-data difference

for m=1:length(data_id)
    
    for n=1:length(Type)
        
        in=find(data_info.Type==n);
        
        Mdata(in,m)=mean(data(in,m));
    
        Mdiff(in,m)=mean(data(in,m)-data(in,d));

    end
    
end

clf
subplot(311)
plot((data(:,d)-Mdata(:,d))./sqrt(R))
title('unbiased observations, normalized with sqrt(R)')

subplot(312)
plot((data(:,b)-Mdata(:,b))./sqrt(R))
title('unbiased baseline, normalized with sqrt(R)')

subplot(313)
plot((data(:,b)-data(:,d)-Mdiff(:,b))./sqrt(R))
title('unbiased observations-baseline difference, normalized with R')

%% 

ir=find(abs((data(:,d)-Mdata(:,d))./sqrt(R))>3 | ...
        abs((data(:,b)-Mdata(:,b))./sqrt(R))>3 | ...
        abs((data(:,b)-data(:,d)-Mdiff(:,b))./sqrt(R))>3);

data(ir,:)=[];

tmp=fieldnames(data_info);

for n=1:length(tmp)

    eval(['data_info.' tmp{n} '(ir)=[];'])

end

G(ir,:)=[];
R(ir)=[];
y(ir)=[];
clear M* i* m n t*

%% 

% Build uncertainty matrix, Eq. 9 in MFT05.
GtoR=G';

for i=1:length(p)

    GtoR(i,:)=GtoR(i,:)./R';

end

P = inv(GtoR*G);

%% 

% Estimate optimized parameters, Eq. 8 in MFT05.
e=P*GtoR*y;

%% 

% Display optimized linear comination.
exps=unique([b;p(:)]);
coef=0*exps;
coef(find(exps==b))=1;

for i=1:size(p,1)
    n=find(exps==p(i,1));
    coef(n)=coef(n)+e(i);
    n=find(exps==p(i,2));
    coef(n)=coef(n)-e(i);
end

disp('Optimized solution is:')
for n=1:length(exps)
    if n==1
        disp(['     ' data_id{exps(n)} ' * ' num2str(coef(n))])
    else
        disp(['   + ' data_id{exps(n)} ' * ' num2str(coef(n))])
    end
end

%% 

% Display cost (Eq. 5 in MFT05) for sensitivity
% experiments and optimized solution.
disp('Cost function:')

for n=1:length(exps)
    J(exps(n))=sum((data(:,exps(n))-data(:,d)).^2./R);
    disp(['   ' data_id{exps(n)} ' J = ' int2str(J(n))])
end

yopt=G*e; % optimized model-data difference
Jo=sum((data(:,b)+yopt-data(:,d)).^2./R);

disp(['   optimized J = ' int2str(Jo)])
disp(['Percent cost reduction: ' num2str(100*((J(b)-Jo)/Jo))])

%%
