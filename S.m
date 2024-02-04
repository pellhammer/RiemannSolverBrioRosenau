function Y = S(uv,a,mu)

nu = 2*((-uv(2)*mu.^2-uv(1)*mu*(a-1)+uv(2))./((a-2)*mu.^2+1));

Y = -2*uv(2)*mu + 2*U(uv,a,mu);


end
