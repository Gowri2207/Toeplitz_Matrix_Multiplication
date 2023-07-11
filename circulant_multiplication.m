function P = circulant_multiplication(C,G)
    n1=size(C,1);
    v=zeros(n1,1);
    n=size(G,2);
    c=C(:,1);
    C_hat=fft(c);
    P=zeros(n,n);
    for i=1:n
        v(1:n)=G(:,i);
        Fn_v=fft(v);
        Y=Fn_v .* C_hat;
        r = ifft(Y);
        P(:,i)=r(1:n);
    end
end

