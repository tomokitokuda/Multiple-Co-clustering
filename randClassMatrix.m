function Z=randClassMatrix(N)

% Assinge maximum number of clusters depending on sample size N
MaxG=ceil(min([log(N)+5,N]));

% Assign a random class label to each data.
Z=full(sparse(1:N,randi(MaxG,[1,N]),1,N,MaxG));

end

