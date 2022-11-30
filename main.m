% Programme Analyse en Composantes Principales

prompt='Entrer la matrice de taille nxp correspondant au tableau de données '
X=input(prompt);

prompt='Entrer la liste contenant le nom de tous les individus '
Nom_individus=input(prompt); 

prompt='Entrer la liste contenant le nom de toutes les variables '
Nom_variables=input(prompt); 

ligne=size(X,1); %nombre de lignes de la matrice X
colonne=size(X,2); %nombre de colonnes de la matrice X

prompt='Entrer la matrice des poids (si les poids sont égaux on entrera : eye(ligne)/ligne '
Poids=input(prompt) 

% Tout d'abord on doit centrer et réduire la matrice X, on nottera Y cette matrice.

Barycentre=diag(Poids)'*X; % barycentre de la matrice X
Y=X-Barycentre; % Y est la matrice X centrée

V=Y'*Poids*Y; %matrice de variance-covariance de Y

% On réduit maintenant la matrice Y, on obtient alors la matrice centrée réduite de X  
Y=Y/sqrt(diag(diag(V))); 

M=eye(colonne); %matrice identité de taille nxn
[U,lambda]=eig(V*M); % U stocke les vecteurs propres du tableau,et lambda stocke les vecteurs 
%propres associées

bar(sort(lambda,'descend')); % trace les valeurs propres dans un graphique
title('Ebouli des valeurs propres') % On donne un nom au graphique tracé
print('Ebouli_des_valeurs_propres','-djpeg') % On enregistre le graphique 

%Cette représentation permet de déterminer la proportion d'information 
% contenue dans un plan. Il suffit d'additionner le pourcentage
% d'inertie des deux axes considérés

Contribution_axe_k=diag(lambda./(sum(diag(lambda)))); %vecteur des contributions
% de l'axe k à l'inertie  totale du nuage

pareto(Contribution_axe_k) %diagramme de Pareto de la contribution de chaque axe à l inertie du nuage
title('Diagramme de Pareto de la contribution relative des axes à l inertie du nuage')
print('Diagramme_de_Pareto_de_la_contribution_relative_des_axes_à_l_inertie_du_nuage','-djpeg')

C=Y*M*U; %matrice des composantes principales

% On va construire une matrice qui contient les indices ponctuels de
% qualité

Ind_ponct_qual=zeros(ligne,colonne); %on initialise le tableau de données de l'indice ponctuel de qualité

for i=1:ligne
 for j=1:colonne
 Ind_ponct_qual(i,j)=(C(i,j)/(norm(C(i,:))))^2;
 end
end


% On définit ici, la qualité de représentation
% des individus sur le plan factoriel

Qual_indiv=Ind_ponct_qual(:,1)+Ind_ponct_qual(:,2);

%matrice des coordonnées des variables Y dans le plan principal

S=U*sqrt(lambda);

Contribution_var=U.^2; %contribution relative des variables aux axes
Qual_var=S.^2; %qualité de représentation des variables

Contribution_indiv_axe_k=Poids*C.^2/lambda; %contribution relative d'un individu i à l'axe k
Inertie_nuage=sum(sum(Poids*Y.^2)); %calcul de l'inertie totale du nuage

Contribution_indiv_nuage=sum(Poids*Y.^2/Inertie_nuage,2); %contribution relative 
% des individus à l'intertie du nuage


%représentation du nuage des individus selon deux axes factoriels Ck 
% qui représentent le mieux qualitativement les données
scatter(C(:,1),C(:,2),'filled') 
text(C(:,1),C(:,2),Nom_individus)
title('Représentation des individus')
print('Représentation_des_individus','-djpeg')

%représentation des variables (vecteurs) selon deux composantes principales réduites
compass(S(:,1),S(:,2),'blue')
text(S(:,1),S(:,2),Nom_variables)
title('Représentation des variables')
print('Représentation_des_variables','-djpeg')

%représentation simultanée des individus et des variables
scatter(C(:,1),C(:,2),'filled')
compass(S(:,1),S(:,2),'blue')
text(C(:,1),C(:,2),Nom_individus)
text(S(:,1),S(:,2),Nom_variables)
title('Représentation simultanée des individus et des variables')
print('Représentation_simultanée_des_individus_et_des_variables','-djpeg')

%cercle des corrélations
compass(U(:,1),U(:,2),'blue')
text(U(:,1),U(:,2),Nom_variables)
title('Cercle Des Correlations')
print('Cercle_des_correlations','-djpeg')
