function testColorful(A, typeMeasure, m, i)

if (all(dfsearch(digraph(A), 1)) ~= 1) 
    error("Grafo non connesso");
end
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

outdegree=zeros(1,n);
indegree=zeros(1,n);
for j=1:n
outdegree(j)= sum(A(j,:));
indegree(j)= sum(A(:,j));
end

switch true
    case i == 1
        Convergenza="$\beta\to 0^+$";
    case i == 2
        Convergenza="$\beta\to +\infty$";
    case i == 3
        Convergenza="$\beta\to \frac{1}{\lambda_1}$";
end

ESC= @(x) diag(expm(x*A));
ETC= @(x) expm(x.*A)*ones(n,1);
RC= @(x) diag(inv(eye(n)-x*A));
K= @(x) (eye(n)-x*A)\ones(n,1);
TCr= @(x) expm(x.*A')*ones(n,1); 
Kr= @(x) (eye(n)-x*A')\ones(n,1);

%--- Scelta dei valori del parametro beta ---
beta=zeros(1,m+1);
betamax=30*log(10)/(norm(A));
for j=1:m
    if (i==1) beta(j)=1/j^(1.7); end
    if (i==2 || i==3) beta(j)=1-1/j^(1.7); end
end
barrier(1)=0; barrier(2)=betamax; barrier(3)=1/max(abs(eigenvalues));

%--- Scelta del vettore della centralità limite ---
if (strcmp(typeMeasure, "RC") || strcmp(typeMeasure,"K")) scalar=1/max(abs(eigenvalues))-eps; end
if (strcmp(typeMeasure,"ESC") || strcmp(typeMeasure,"ETC")) scalar=betamax-eps; end
if (i==1) beta(m+1)=beta(m)+(beta(m)-beta(1))/5; end
if (i==2 || i==3) beta(m+1)=beta(m)+(beta(m)-beta(1))/5; end
beta=scalar*beta;

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

%--- Grafico dell'andamento della classifica ---
B=zeros(n,m+1);
for j=1:m
    u=Measure(beta(j))';
    [~, rankedIndices] = sort(u, 'descend');
    for k=1:n
        B(rankedIndices(k),j)=1-(k-1)/n;
    end
end
B(:,m+1)=B(:,m);

[~, rankedIndices] = sort(v, 'descend');
l=1; colore=zeros(1,n);
for k=1:(n-1)
    colore(rankedIndices(k))=1-(l-1)/n;
    if (v(rankedIndices(k))~=v(rankedIndices(k+1))) l=l+1; end
end
colore(rankedIndices(n))=1-(l-1)/n;
colore =(colore-min(colore))/(max(colore)-min(colore));

cmap = jet(256);
figure(1);
xlabel("Valori del parametro");
ylabel("Evoluzione dei rankings");
title("Convergenza della misura " + typeMeasure + " nel caso " + Convergenza + doubleCase, "Interpreter", "latex");
hold on
for j=1:n
    colorIndex = round(colore(j) * (size(cmap, 1) - 1)) + 1;
    color = cmap(colorIndex, :);
    plot(beta, B(j,:), 'Color', color, 'LineWidth', 1);
end
plot([barrier(i),barrier(i)], [-0.05,1.05], 'black', 'LineWidth', 1.5);
ylim([-0.05,1.05]);

if (doubleCase == " (broadcast centrality)")
    figure(2);
    hold on
    for j=1:n
        colorIndex = round(colore(j) * (size(cmap, 1) - 1)) + 1;
        color = cmap(colorIndex, :);
        plot(beta, B(j,:), 'Color', color, 'LineWidth', 1);
    end
    plot([barrier(i),barrier(i)], [-0.05,1.05], 'black', 'LineWidth', 1.5);
    ylim([-0.05,1.05]); 
    xlabel("Valori del parametro");
    ylabel("Angolo dei vettori di ranking");
    title("Convergenza della misura " + typeMeasure + " nel caso " + Convergenza + " (authority centrality)", "Interpreter", "latex");
end
end
