function Y = U(uv,a,mu)
%V value of Hugoniot
nu = 2*((-uv(2)*mu.^2-uv(1)*mu*(a-1)+uv(2))./((a-2)*mu.^2+1));

Y = uv(1) + mu.*nu;

end
