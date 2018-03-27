function Model=coreVBCCGaussS(Xg,Xp, Xb, EZ,EVg,EVp, EVb, AZ0,AV0,AVW0, L0,M0,G0,S0,gamma1, gamma2, beta1, MaxS,M ,onefcluster)

%% Preprocessing
T = nanmax(nanmax(reshape(Xb, 1, [])), 1);% Max number of categories
EV = EVg; % changing name for Gaussian
X = Xg;



% Get size of data. If it is empty, it is replaced by N x 1 nan vector.
if ~isempty(X)
 [N,D]=size(X);
end
if ~isempty(Xp)
 [N,Dp]=size(Xp);
end
if ~isempty(Xb)
 [N, Db] = size(Xb);
end

if isempty(X)
  X = nan*ones(N, 1);
  D = 1;
end
if isempty(Xp)
  Xp = nan*ones(N, 1);
  Dp = 1;
end
if isempty(Xb)
  Xb = nan*ones(N, 1);
  Db = 1;
end


% Multiple
iniViewWeight = ones(1, M)/M;

% Append additional classes
% Subject
MaxG = cell(1, M);
for m=1:M
    EZ{m} = [EZ{m},zeros(size(EZ{m}))];
    MaxG{m}=size(EZ{m},2);
end

% Feature
% Gauss
MaxC = cell(1, M);
for m=1:M
    if onefcluster==0
        EV{m} = iniViewWeight(m)*[EV{m};zeros(size(EV{m}))];
    else
        EV{m} = iniViewWeight(m)*EV{m};
    end
    MaxC{m} =size(EV{m},1);
end

% Poiss
MaxCp = cell(1, M);
for m=1:M
    if onefcluster==0
        EVp{m} = iniViewWeight(m)*[EVp{m};zeros(size(EVp{m}))];
    else
        EVp{m} = iniViewWeight(m)*EVp{m};
    end
    MaxCp{m} =size(EVp{m},1);
end

% Bernouilli
MaxCb = cell(1, M);
for m=1:M
    if onefcluster==0
        EVb{m} = iniViewWeight(m)*[EVb{m};zeros(size(EVb{m}))];
    else
        EVb{m} = iniViewWeight(m)*EVb{m};
    end
    MaxCb{m} =size(EVb{m},1);
end


% Guassian
E1=ones(N,D);
E1(isnan(X))=0;
EX=X;
EX(isnan(X))=0;
EX2=EX.^2;

% Poisson
E1p = ones(N, Dp);
E1p(isnan(Xp))=0;
EXp = Xp;
EXp(isnan(Xp))=0;
EgammalnXp = gammaln(EXp+1);

% Bernouilli
EXb = cell(1, T);
for t=1:T
    EXb{t} = (Xb==t); % make one if the category is t.
    EXb{t}(isnan(Xb))=0;  % zeri fir Nan in Bernoulli
end

%% Main iteration of VB algorithm
sLE=zeros(1,MaxS);
oldLE=-inf;
cnv_f=0;

SZ = cell(1, M);
SV = cell(1, M);
SVp = cell(1, M);
SVb = cell(1, M);
S1Z = cell(1, M);
S1V = cell(1, M);
SZV1 = cell(1, M);
SZVX = cell(1, M);
SZVXX = cell(1, M);
S1Vp = cell(1, M);
SZV1p = cell(1, M);
SZVXp = cell(1, M);
SZVgammalnXp = cell(1, M);
S1Vb = cell(1, M);
SZVXb = cell(M, T);
RZP = cell(1, M);
AZP = cell(1, M);
ELRZ = cell(1, M);
EL1RZ = cell(1, M);
RVP = cell(1, M);
AVP = cell(1, M);
ELRV = cell(1, M);
EL1RV = cell(1, M);
RVPp = cell(1, M);
AVPp = cell(1, M);
ELRVp = cell(1, M);
EL1RVp = cell(1, M);
RVPb = cell(1, M);
AVPb = cell(1, M);
ELRVb = cell(1, M);
EL1RVb = cell(1, M);
LP = cell(1, M);
MP = cell(1, M);
GP = cell(1, M);
ELS = cell(1, M);
SP = cell(1, M);
gamma1p = cell(1, M);
gamma2p = cell(1, M);
beta1p = cell(M, T);

sumbeta1 = cell(1, M);
LEs = nan*ones(1, M);

LEV = cell(1, M);
LEVp = cell(1, M);
LEVb = cell(1, M);

EZ2{m} = cell(1, M);
EV2{m} = cell(1, M);
EV2p{m} = cell(1, M);
EV2b{m} = cell(1, M);
MaxG2{m} = cell(1, M);
MaxC2{m} = cell(1, M);
MaxC2p{m} = cell(1, M);
MaxC2b{m} = cell(1, M);


for s=1:MaxS
    % Starting iterations view by view
    for m=1:M
        SZ{m}=sum(EZ{m},1);
        [SZ{m},gidx]=sort(SZ{m},'descend');
        EZ{m}=EZ{m}(:,gidx);
        % Gaussian
        SV{m}=sum(EV{m},2);
        [SV{m},cidx]=sort(SV{m},'descend');
        EV{m}=EV{m}(cidx,:);
        % Poisson
        SVp{m}=sum(EVp{m},2);
        [SVp{m},cidxp]=sort(SVp{m},'descend');
        EVp{m}=EVp{m}(cidxp,:);
        % Bernouilli
        SVb{m}=sum(EVb{m},2);
        [SVb{m},cidxb]=sort(SVb{m},'descend');
        EVb{m}=EVb{m}(cidxb,:);
        
        
        % Compute the sufficient statistics
        SZ{m} = transpose(SZ{m});
        S1Z{m} = transpose(sum(1-EZ{m}*triu(ones(MaxG{m})),1));
        % Gaussian
        SV{m}=transpose(SV{m});
        S1V{m}=transpose(sum(1-tril(ones(MaxC{m}))*EV{m},2));
        SZV1{m}=EZ{m}'*E1*EV{m}';
        SZVX{m}=EZ{m}'*EX*EV{m}';
        SZVXX{m}=EZ{m}'*EX2*EV{m}';
        % Poisson
        SVp{m}=transpose(SVp{m});
        S1Vp{m}=transpose(sum(1-tril(ones(MaxCp{m}))*EVp{m},2));
        SZV1p{m}=EZ{m}'*E1p*EVp{m}';
        SZVXp{m}=EZ{m}'*EXp*EVp{m}';
        SZVgammalnXp{m}=EZ{m}'*EgammalnXp*EVp{m}';
        % Bernouilli
        SVb{m}=transpose(SVb{m});
        S1Vb{m}=transpose(sum(1-tril(ones(MaxCb{m}))*EVb{m},2));
        %SZV1b{m}=EZ{m}'*E1b*EVb{m}'; % not used
        
        for t=1:T
            SZVXb{m}{t}=EZ{m}'*EXb{t}*EVb{m}'; % total number of category t
        end
        
        % Update the posterior hyperparameters
        RZP{m}=1+SZ{m};
        AZP{m}=AZ0+S1Z{m};
        ELRZ{m}=psi(RZP{m})-psi(RZP{m}+AZP{m});
        EL1RZ{m}=psi(AZP{m})-psi(RZP{m}+AZP{m});
        % Gaussian
        RVP{m}=1+SV{m};
        AVP{m}=AV0+S1V{m};
        ELRV{m}=psi(RVP{m})-psi(RVP{m}+AVP{m});
        EL1RV{m}=psi(AVP{m})-psi(RVP{m}+AVP{m});
        % Poiss
        RVPp{m}=1+SVp{m};
        AVPp{m}=AV0+S1Vp{m}; % same concenraion AVO as Gaussian
        ELRVp{m}=psi(RVPp{m})-psi(RVPp{m}+AVPp{m});
        EL1RVp{m}=psi(AVPp{m})-psi(RVPp{m}+AVPp{m});
        % Bernouilli
        RVPb{m}=1+SVb{m};
        AVPb{m}=AV0+S1Vb{m}; % same concenraion AVO as Gaussian
        ELRVb{m}=psi(RVPb{m})-psi(RVPb{m}+AVPb{m});
        EL1RVb{m}=psi(AVPb{m})-psi(RVPb{m}+AVPb{m});
        
        % Gaussian
        LP{m}=L0+SZV1{m};
        MP{m}=(L0*M0+SZVX{m})./LP{m};
        GP{m}=G0+SZV1{m};
        SP{m}=(G0*S0+L0*(M0.^2)+SZVXX{m}-LP{m}.*(MP{m}.^2))./GP{m};
        ELS{m}=-log(SP{m})-log(GP{m}/2)+psi(GP{m}/2);
        
        % Poisson
        gamma1p{m} = gamma1 + SZVXp{m};
        gamma2p{m} = gamma2 + SZV1p{m};
        
        % Bernouilli
        for t=1:T
            beta1p{m}{t} = beta1 + SZVXb{m}{t};
            %beta2p{m} = beta2 + SZVXb2{m};
        end
    end
    
    % For views
    clear EVW  % SV for view
    clear EVWp;
    clear EVWb;
    for m=1:M
        EVW(m, :)= sum(EV{m}, 1);
        EVWp(m, :) = sum(EVp{m}, 1); % Poisson
        EVWb(m, :) = sum(EVb{m}, 1); % Bernouilli
    end
    
    
    EVW = [EVW EVWp EVWb]; % Gauss + Poisson + Bernouillil;which feature belongs to which view
    SVW = sum(EVW, 2);
    SVW = transpose(SVW);
    MaxCC = size(EVW, 1);
    S1VW=transpose(sum(1-tril(ones(MaxCC))*EVW,2));
    
    RVWP=1+SVW;
    AVWP=AVW0+S1VW;
    ELRVW=psi(RVWP)-psi(RVWP+AVWP);
    EL1RVW=psi(AVWP)-psi(RVWP+AVWP);
    
    HVW = 0; % This part is absorved in HV
    LE = sum((ELRVW+EL1RVW*triu(ones(MaxCC),+1))*EVW) + HVW;
    for m=1:M
        LL = 0;
        LL=LL+sum(EZ{m}*(ELRZ{m}+tril(ones(MaxG{m}),-1)*EL1RZ{m})) ...
            +sum((ELRV{m}+EL1RV{m}*triu(ones(MaxC{m}),+1))*EV{m})...
            +sum((ELRVp{m}+EL1RVp{m}*triu(ones(MaxCp{m}),+1))*EVp{m})... % Pois
            +sum((ELRVb{m}+EL1RVb{m}*triu(ones(MaxCb{m}),+1))*EVb{m}); %Ber
        
        %beta1p{m}
        %beta2p{m}
        
        LL=LL-(1/2)*( ...
            sum(sum(SZVXX{m}.*(1./SP{m}))) ...
            -2*sum(sum(SZVX{m}.*(MP{m}./SP{m}))) ...
            +sum(sum(SZV1{m}.*((MP{m}.^2)./SP{m}+(1./LP{m})-ELS{m}+log(2*pi)))))...
            +sum(sum(SZVXp{m}.*(psi(gamma1p{m})-log(gamma2p{m}))))... % Pois
            + sum(sum(SZV1p{m}.*(-gamma1p{m}./gamma2p{m})))... % Pois
            + sum(sum(-SZVgammalnXp{m}));% Pois
        sumbeta1{m} = 0;
        for t=1:T
            sumbeta1{m} =sumbeta1{m} + beta1p{m}{t};
        end
        for t=1:T
            LL = LL +  sum(sum(SZVXb{m}{t}.*(psi(beta1p{m}{t})-psi(sumbeta1{m}))));  % Bernoulli
        end
        
        HZ=-sum(sum(log(EZ{m}.^EZ{m}))); %%%%%%% % Objects
        HV=-sum(sum(log(EV{m}.^EV{m}))); %%%%%%% % Gaussian
        HP=-sum(sum(log(EVp{m}.^EVp{m}))); %%%%%%% % Pois
        HB=-sum(sum(log(EVb{m}.^EVb{m}))); %%%%%%% % Bernouilli
        
        DRZ=sum((RZP{m}-1).*ELRZ{m}+(AZP{m}-AZ0).*EL1RZ{m} ...
            +(gammaln(RZP{m}+AZP{m})-gammaln(RZP{m})-gammaln(AZP{m})) ...
            -(gammaln(1+AZ0)-gammaln(1)-gammaln(AZ0)));
        
        DRV=sum((RVP{m}-1).*ELRV{m}+(AVP{m}-AV0).*EL1RV{m} ...
            +(gammaln(RVP{m}+AVP{m})-gammaln(RVP{m})-gammaln(AVP{m})) ...
            -(gammaln(1+AV0)-gammaln(1)-gammaln(AV0)));
        
        DRVp=sum((RVPp{m}-1).*ELRVp{m}+(AVPp{m}-AV0).*EL1RVp{m} ...
            +(gammaln(RVPp{m}+AVPp{m})-gammaln(RVPp{m})-gammaln(AVPp{m})) ...
            -(gammaln(1+AV0)-gammaln(1)-gammaln(AV0))); % Pois
        
        DRVb=sum((RVPb{m}-1).*ELRVb{m}+(AVPb{m}-AV0).*EL1RVb{m} ...
            +(gammaln(RVPb{m}+AVPb{m})-gammaln(RVPb{m})-gammaln(AVPb{m})) ...
            -(gammaln(1+AV0)-gammaln(1)-gammaln(AV0))); % Ber
        
        
        DM=sum(sum((L0/2)*(((MP{m}-M0).^2)./SP{m}) ...
            +(1/2)*((L0./LP{m})-log(L0./LP{m})-1)));
        
        DS=sum(sum((G0/2).*((S0./SP{m})-log(S0./SP{m})-1)...
            +(G0/2).*(log(GP{m}./G0)+psi(G0/2)-psi(GP{m}/2)) ...
            +(gammaln(G0/2)+(G0/2)-(G0/2).*psi(G0/2)) ...
            -(gammaln(GP{m}/2)+(GP{m}/2)-(GP{m}/2).*psi(GP{m}/2))));
        % Poisson
        DBp = sum(sum(gamma1p{m}.*log(gamma2p{m})-gammaln(gamma1p{m})...
            + (gamma1p{m}-1).*(psi(gamma1p{m})-log(gamma2p{m}))-gamma1p{m}))...
            - (sum(sum((gamma1*log(gamma2)-gammaln(gamma1))...
            + (gamma1-1)*(psi(gamma1p{m})-log(gamma2p{m}))...
            -gamma2*gamma1p{m}./gamma2p{m})));
        % Bernouilli
        %         DBb=sum(sum((beta1p{m}-beta1).*(psi(beta1p{m})-psi(beta1p{m}+beta2p{m}))+(beta2p{m}-beta2).*(psi(beta2p{m})-psi(beta1p{m}+beta2p{m})) ...
        %             +(gammaln(beta1p{m}+beta2p{m})-gammaln(beta1p{m})-gammaln(beta2p{m})) ...
        %             -(gammaln(beta1+beta2)-gammaln(beta1)-gammaln(beta2))));
        %
        % Bernouilli
        DBb = 0;
        for t=1:T
            DBb = DBb + sum(sum((beta1p{m}{t}-beta1).*(psi(beta1p{m}{t})-psi(sumbeta1{m}))));
        end
        DBb = DBb + sum(sum(gammaln(sumbeta1{m}) - gammaln(T*beta1)));
        for t=1:T
            DBb = DBb-sum(sum(gammaln(beta1p{m}{t}) - gammaln(beta1)));
        end
        
        LEs(m)=(LL+HZ+HV+HP+HB)-(DRZ+DRV+DRVp+DRVb+DM+DS+DBp+DBb);
    end
    
    LE = LE + sum(LEs);
    % Draw learning curve
    sLE(s)=LE;
    %if(reports)
    %    figure(gcf);
    %    plot(1:s,sLE(1:s));
    %    drawnow;
    %end
    
    % Judge the convergence
    if((LE-oldLE)<(1e-10))
        cnv_f=cnv_f+1;
    else
        cnv_f=0;
    end
    
    oldLE=LE;
    
    % Update the posterior subject-class assignment
    for m=1:M
        LEZ=ones(N,1)*ELRZ{m}'+ones(N,1)*(tril(ones(MaxG{m}),-1)*EL1RZ{m})';
        LEZ=LEZ-(1/2)* ( ...
            EX2*EV{m}'*(1./SP{m})' ...
            -2*EX*EV{m}'*(MP{m}./SP{m})' ...
            +E1*EV{m}'*(((MP{m}.^2)./SP{m})+(1./LP{m})-ELS{m})')...
            +EXp*EVp{m}'*(psi(gamma1p{m})-log(gamma2p{m}))'... % Pois
            +E1p*EVp{m}'*(-gamma1p{m}./gamma2p{m})'... % Pois
            -EgammalnXp*EVp{m}'*ones(size(gamma1p{m}'));% Pois
        for t=1:T
            LEZ = LEZ + EXb{t}*EVb{m}'*(psi(beta1p{m}{t})-psi(sumbeta1{m}))'; % Ber
        end
        
        LEZ=LEZ-repmat(max(LEZ,[],2),[1 MaxG{m}]);
        EZ{m}=exp(LEZ);
        EZ{m}=EZ{m}./(sum(EZ{m},2)*ones(1,MaxG{m}));
    end
    
    % Update the posterior feature-class assignment
    clear LEV;
    clear LEVp;
    clear LEVb;
    for m=1:M
        % For DP of views: same values for all pairs (feature, f.cluster)
        if m==1
            LEV{m} =ELRVW(:, m)*ones(size(ELRV{m}))'*ones(1, D);
            % Poisson
            LEVp{m} =ELRVW(:, m)*ones(size(ELRVp{m}))'*ones(1, Dp);
            % Ber
            LEVb{m} =ELRVW(:, m)*ones(size(ELRVb{m}))'*ones(1, Db);
        else
            LEV{m} =ELRVW(:, m)*ones(size(ELRV{m}))'*ones(1, D) + ...
                + sum(EL1RVW(:, 1:(m-1)))*ones(size(ELRV{m}))'*ones(1,D);
            % Poisson
            LEVp{m} =ELRVW(:, m)*ones(size(ELRVp{m}))'*ones(1, Dp) + ...
                + sum(EL1RVW(:, 1:(m-1)))*ones(size(ELRVp{m}))'*ones(1,Dp);
            % Bernouilli
            LEVb{m} =ELRVW(:, m)*ones(size(ELRVb{m}))'*ones(1, Db) + ...
                + sum(EL1RVW(:, 1:(m-1)))*ones(size(ELRVb{m}))'*ones(1,Db);
        end
        
        LEV{m}=LEV{m} + ELRV{m}'*ones(1,D)+tril(ones(MaxC{m}),-1)*EL1RV{m}'*ones(1,D);
        LEV{m}=LEV{m}-(1/2)*( ...
            (1./SP{m})'*EZ{m}'*EX2 ...
            -2*(MP{m}./SP{m})'*EZ{m}'*EX ...
            +(((MP{m}.^2)./SP{m})+(1./LP{m})-ELS{m})'*EZ{m}'*E1);
        
        % Poisson
        LEVp{m}= LEVp{m} + ELRVp{m}'*ones(1,Dp)+tril(ones(MaxCp{m}),-1)*EL1RVp{m}'*ones(1,Dp);
        LEVp{m} = LEVp{m} + (psi(gamma1p{m})-log(gamma2p{m}))'*EZ{m}'*EXp...
            +(-gamma1p{m}./gamma2p{m})'*EZ{m}'*E1p...
            -ones(size(gamma1p{m}'))*EZ{m}'*EgammalnXp;
        % Bernouilli
        LEVb{m} = LEVb{m} + ELRVb{m}'*ones(1,Db)+tril(ones(MaxCb{m}),-1)*EL1RVb{m}'*ones(1,Db);
        for t=1:T
            LEVb{m} = LEVb{m} + (psi(beta1p{m}{t})-psi(sumbeta1{m}))'*EZ{m}'*EXb{t};
        end     
    end
    
    
    % Subtract the maximum for each feature using all m=1:M
    % Gaussian
    clear maxall;
    clear maxallp;
    clear maxallb;
    for m=1:M
        maxall(m, :) = max(LEV{m});
        maxallp(m, :) = max(LEVp{m}); % Pois
        maxallb(m, :) = max(LEVb{m}); % Ber
    end
    % size(maxall)
    if size(maxall, 1)>1
        maxall2 = max(maxall); % 1 x D matrix
    else
        maxall2 = maxall; % 1 x D matrix
    end
    
    if size(maxallp, 1)>1
        maxall2p = max(maxallp); % 1 x D matrix
    else
        maxall2p = maxallp; % 1 x D matrix
    end
    
    if size(maxallb, 1)>1
        maxall2b = max(maxallb); % 1 x D matrix
    else
        maxall2b = maxallb; % 1 x D matrix
    end
    
    for m=1:M
        EV{m}=exp(LEV{m}-(ones(MaxC{m},1)*maxall2));
        EVp{m}=exp(LEVp{m}-(ones(MaxCp{m},1)*maxall2p)); % Pois
        EVb{m}=exp(LEVb{m}-(ones(MaxCb{m},1)*maxall2b)); % Ber
    end
    
    % Normalize using all m
    sumall = zeros(size(sum(EV{1}, 1)));
    sumallp = zeros(size(sum(EVp{1}, 1))); %
    sumallb = zeros(size(sum(EVb{1}, 1))); % Pois
    for m=1:M
        % sum feature-wise sumall is 1 x  D vector
        sumall= sumall + sum(EV{m}, 1);
        sumallp= sumallp + sum(EVp{m}, 1);  % Pois
        sumallb= sumallb + sum(EVb{m}, 1);  % Ber
    end
    
    for m=1:M
        EV{m}=EV{m}./(ones(MaxC{m},1)*sumall);
        EVp{m}=EVp{m}./(ones(MaxCp{m},1)*sumallp); % Pois
        EVb{m}=EVb{m}./(ones(MaxCb{m},1)*sumallb); % Ber
    end
    
    % Changing the order of m: from large cluster size of view to small
    SVW=sum(EVW,2);
    [SVW,ccidx]=sort(SVW,'descend');
    EVW = EVW(ccidx,:);
    for m=1:M
        EZ2{m} = EZ{ccidx(m)};
        EV2{m} = EV{ccidx(m)};
        EV2p{m} = EVp{ccidx(m)};
        EV2b{m} = EVb{ccidx(m)};
        MaxG2{m} = MaxG{ccidx(m)};
        MaxC2{m} = MaxC{ccidx(m)};
        MaxC2p{m} = MaxCp{ccidx(m)};
        MaxC2b{m} = MaxCb{ccidx(m)};
    end
    
    for m=1:M
        EZ{m} = EZ2{m};
        EV{m} = EV2{m};
        EVp{m} = EV2p{m};
        EVb{m} = EV2b{m};
        MaxG{m} = MaxG2{m};
        MaxC{m} = MaxC2{m};
        MaxCp{m} = MaxC2p{m};
        MaxCb{m} = MaxC2b{m};
    end
    
    
    % Check
    add = 0;
    for m=1:M
        add = add + sum(EVb{m});
    end
     
    if(cnv_f==5)
        break;
    end
    
end

%% Post-processes
%Model.SVW = SVW;
%Model.AVWP = AVWP;

M2 = sum(SVW>=1e-5); % Only non-zero views

%Model.EVW = EVW(1:M2, :);
MapV = cell(1, M2);
MapVp = cell(1, M2);
MapVb = cell(1, M2);
for m=1:M2
    SZ{m}=sum(EZ{m},1);
    gidx=SZ{m}>=1e-10;
    %gidx=SZ{m}>=1e-10;
    Model.EZ{m}=EZ{m}(:,gidx);
    Model.SZ{m}=SZ{m}(:,gidx);
    [Model.MaxZ{m},Model.MapZ{m}]=max(EZ{m},[],2);
    
    SV{m}=sum(EV{m},2);
    cidx=SV{m}>=1e-10;
    Model.EV{m}=EV{m}(cidx,:);
    Model.SV{m}=SV{m}(cidx,:);
    [Model.MaxV{m}, MapV{m}]=max(EV{m},[],1);
    %Model.RZP{m}=RZP{m}(gidx,:);
    %Model.AZP{m}=AZP{m}(gidx,:);
    %Model.RVP{m}=RVP{m}(:,cidx);
    %Model.AVP{m}=AVP{m}(:,cidx);
    
    Model.LP{m}=LP{m}(gidx,cidx);
    Model.MP{m}=MP{m}(gidx,cidx);
    Model.GP{m}=GP{m}(gidx,cidx);
    Model.SP{m}=SP{m}(gidx,cidx);
    
    % Poisson
    SVp{m} = sum(EVp{m},2);
    cidx=SVp{m}>=1e-10;
    Model.EVp{m}=EVp{m}(cidx,:);
    %Model.SVp{m}=SVp{m}(cidx,:);
    [Model.MaxVp{m}, MapVp{m}]=max(EVp{m},[],1);
   
    Model.gamma1p{m}=gamma1p{m}(gidx,cidx);
    Model.gamma2p{m}=gamma2p{m}(gidx,cidx);
    
    % Berenouilli
    SVb{m} = sum(EVb{m},2);
    cidx=SVb{m}>=1e-10;
    Model.EVb{m}=EVb{m}(cidx,:);
    %Model.SVb{m}=SVb{m}(cidx,:);
    [Model.MaxVb{m}, MapVb{m}]=max(EVb{m},[],1);
    %Model.RVPb{m}=RVPb{m}(:,cidx);
    %Model.AVPb{m}=AVPb{m}(:,cidx);
    
    for t=1:T
        Model.beta1p{m}{t}=beta1p{m}{t}(gidx,cidx);
    end
end

% Gaussian
numf = length(Model.MaxV{1});
clear mapp;
for m=1:M2
    mapp(m, :) =  Model.MaxV{m};
end
for j=1:numf
    [~, maxnum] = max(mapp(:, j));
    for m=1:M2
        if m==maxnum
            mapp(m, j) = MapV{m}(j);
        else
            mapp(m, j) = 0;
        end
    end
end
for m=1:M2
    Model.MapV{m} = mapp(m, :);
end

% Poisson
numf = length(Model.MaxVp{1});
clear mapp;
for m=1:M2
    mapp(m, :) =  Model.MaxVp{m};
end
for j=1:numf
    [~, maxnum] = max(mapp(:, j));
    for m=1:M2
        if m==maxnum
            mapp(m, j) = MapVp{m}(j);
        else
            mapp(m, j) = 0;
        end
    end
end
for m=1:M2
    Model.MapVp{m} = mapp(m, :);
end

% Bernouilli
numf = length(Model.MaxVb{1});
clear mapp;
for m=1:M2
    mapp(m, :) =  Model.MaxVb{m};
end
for j=1:numf
    [~, maxnum] = max(mapp(:, j));
    for m=1:M2
        if m==maxnum
            mapp(m, j) = MapVb{m}(j);
        else
            mapp(m, j) = 0;
        end
    end
end
for m=1:M2
    Model.MapVb{m} = mapp(m, :);
end

% Remove unnecessary parts
for m=1:M2
    maxz = max(Model.MapZ{m});
    maxvb = max(Model.MapVb{m});
    for t=1:T
        Model.beta1p{m}{t}=Model.beta1p{m}{t}(1:maxz, 1:maxvb);
    end
    
    maxv = max(Model.MapV{m});
    Model.LP{m}= Model.LP{m}(1:maxz, 1:maxv);
    Model.MP{m}= Model.MP{m}(1:maxz, 1:maxv);
    Model.GP{m}= Model.GP{m}(1:maxz, 1:maxv);
    Model.SP{m}= Model.SP{m}(1:maxz, 1:maxv);
    
    maxvp = max(Model.MapVp{m});
    Model.gamma1p{m}= Model.gamma1p{m}(1:maxz, 1:maxvp);
    Model.gamma2p{m}= Model.gamma2p{m}(1:maxz, 1:maxvp);
end


Model.L0=L0;
Model.M0=M0;
Model.G0=G0;
Model.S0=S0;
Model.gamma1 = gamma1;
Model.gamma2 = gamma2;
Model.beta1 = beta1; % Symmetric Dirichlet

Model.LE=LE; % Maximum log-likeilhood


end

