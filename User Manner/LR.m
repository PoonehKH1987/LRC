function [W]=LR(X)
y=-0.15.*X(:,1)-0.04.*X(:,2)-0.09*X(:,3)+0.002*X(:,4)-0.02*X(:,5)-0.105*X(:,6)+0.195.*X(:,7)-0.348;
%y=-0.15.*X(:,3)-0.04.*X(:,13)-0.09*X(:,16)+0.002*X(:,18)-0.02*X(:,19)-0.10
%5*X(:,20)+0.195.*X(:,21)-0.348; 
W=exp(y)./(1+exp(y));  