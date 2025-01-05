%This Matlab function was developed to generate simulation results to:
%Unnikrishnan Kunnath Ganesan, Emil Bjrnson and Erik G. Larsson (2021), 
%[1] "Clustering Based Activity Detection Algorithms for Grant-Free Random Access in Cell-Free Massive MIMO", IEEE Transactions in Communications
function [c] = compute_poly_coefficients(am,bm)
    M = length(am) ; 
    c = zeros(1,2*M) ;
    for m=1:1:M
        temp = 1 ; 
        for k=1:1:M
            if k==m
                continue;
            end
            temp = conv(temp,[am(k)^2 2*am(k) 1]) ;
        end
        temp = conv(temp,[am(m)^2 am(m)-bm(m)]) ;
        c = c +  temp ; 
    end
   
end