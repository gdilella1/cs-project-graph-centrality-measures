function testAngle1(A, typeMeasure, m, i)

if (all(dfsearch(digraph(A), 1)) ~= 1) 
    error("Grafo non connesso");
end

%--- Inizializzazione ---
n=size(A,1);
theta=zeros(m,1)';
doubleCase="";

[eigenvectors, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);
[~, idx] = max(abs(eigenvalues));
dominantEigenvector = eigenvectors(:, idx);
dominantEigenvector = dominantEigenvector / norm(dominantEigenvector);
[eigenvectorsLeft, eigenvaluesLeft] = eig(A'); 
eigenvaluesLeft = diag(eigenvaluesLeft);
[~, idx] = max(abs(eigenvaluesLeft));
dominantEigenvectorLeft = eigenvectorsLeft(:, idx);
dominantEigenvectorLeft = dominantEigenvectorLeft / norm(dominantEigenvectorLeft);

switch true
    case i == 1
        Convergenza="$\beta\to 0^+$";
    case i == 2
        Convergenza="$\beta\to +\infty$";
    case i == 3
        Convergenza="$\beta\to \frac{1}{\lambda_1}$";
end

outdegree=zeros(1,n);
indegree=zeros(1,n);
for j=1:n
outdegree(j)= sum(A(j,:));
indegree(j)= sum(A(:,j));
end

%--- Scelta dei valori del parametro beta ---
betamax=30*log(10)/(norm(A));
for j=1:m
    if (i==1) beta(j)=1/j^(0.7); end
    if (i==2 || i==3) beta(j)=1-1/j; end
end
if (strcmp(typeMeasure, "RC") || strcmp(typeMeasure,"K")) scalar=1/max(abs(eigenvalues))-eps; end
if (strcmp(typeMeasure,"ESC") || strcmp(typeMeasure,"ETC")) scalar=betamax-eps; end
beta=scalar*beta;

%--- Grafico dell'andamento dell'angolo ---
if (i==1)
    ESC= @(x) diag(expm(x*A))-1;
    ETC= @(x) (expm(x.*A)-eye(n))*ones(n,1);
    RC= @(x) diag(inv(eye(n)-x*A))-1;
    K= @(x) (inv(eye(n)-x*A)-eye(n))*ones(n,1);
    TCr= @(x) (expm(x.*A')-eye(n))*ones(n,1);
    Kr= @(x) (inv(eye(n)-x*A')-eye(n))*ones(n,1);
end
if (i==2 || i==3)
    ESC= @(x) sqrt(diag(expm(x*A)));
    ETC= @(x) expm(x.*A)*ones(n,1); 
    RC= @(x) sqrt(diag(inv(eye(n)-x*A))); 
    K= @(x) (eye(n)-x*A)\ones(n,1);
    TCr= @(x) expm(x.*A')*ones(n,1); 
    Kr= @(x) (eye(n)-x*A')\ones(n,1);
end

%--- Scelta del vettore della centralità limite ---
if (isequal(A, A'))
    if (strcmpi(typeMeasure, "ESC")==1)
        Measure=ESC; 
        if (i==1) v=outdegree; end 
        if (i==2) v=dominantEigenvector; end
    end
    if (strcmpi(typeMeasure, "ETC")==1)
        Measure=ETC; 
        if (i==1) v=outdegree; end 
        if (i==2) v=dominantEigenvector; end
    end
    if (strcmpi(typeMeasure, "RC")==1)
        Measure=RC; 
        if (i==1) v=outdegree; end 
        if (i==3) v=dominantEigenvector; end
    end
    if (strcmpi(typeMeasure, "K")==1)
        Measure=K; 
        if (i==1) v=outdegree; end 
        if (i==3) v=dominantEigenvector; end
    end
else
    doubleCase=" (broadcast centrality)";
    if (strcmpi(typeMeasure, "ETC")==1)
        Measure=ETC; Measure2=TCr;
        if (i==1) v=outdegree; w=indegree; end 
        if (i==2) v=dominantEigenvector; w=dominantEigenvectorLeft; end
    end
    if (strcmpi(typeMeasure, "K")==1)
        Measure=K; Measure2=Kr;
        if (i==1) v=outdegree; w=indegree; end 
        if (i==3) v=dominantEigenvector; w=dominantEigenvectorLeft; end
    end
end

%--- Grafico dell'andamento dell'angolo ---
for k=1:m
    var=Measure(beta(k));
    var=var/norm(var);
    theta(k)=acos(dot(var,v)/(norm(var)*norm(v)));
end
figure(1);
plot(beta, theta, "bo-");
xlabel("Valori del parametro");
ylabel("Angolo $\theta$ (rad)", "Interpreter", "latex");
title("Convergenza della misura " + typeMeasure + " nel caso " + Convergenza + doubleCase, "Interpreter", "latex");

%--- Doppio grafico nel caso di grafi orientati ---
if strcmp(doubleCase," (broadcast centrality)") 
    for k=1:m
        var=Measure2(beta(k));
        var=var/norm(var);
        theta(k)=acos(dot(var,w)/(norm(var)*norm(w)));
    end
    figure(2);
    plot(beta, theta, "go-");
    xlabel("Valori del parametro");
    ylabel("Angolo $\theta$ (rad)", "Interpreter", "latex");
    title("Convergenza della misura " + typeMeasure + " nel caso " + Convergenza + " (authority centrality)", "Interpreter", "latex");
end
end
