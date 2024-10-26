function [all_coeffs sig_resid]=sigmoid_doseResponse(dose,response)
%hill equation sigmoid
sigmoid=@(beta,x) beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));
%calculate some rough guesses for initial parameters
minResponse=min(response);
maxResponse=max(response);
midResponse=mean([minResponse maxResponse]);
minDose=min(dose);
maxDose=max(dose);

%fit the curve and compute the values
[coeffs,r,J]=nlinfit(dose,response,sigmoid,[minResponse maxResponse midResponse 1]);

all_coeffs=coeffs;
sig_resid=r;