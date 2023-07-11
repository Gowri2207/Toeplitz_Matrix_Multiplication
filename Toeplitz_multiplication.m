
E_n1=[];
E_n2=[];
E_n3=[];
size=[];
T_3=[];
T_2=[];
T_1=[];
T_4=[];
itr=[];
 %size of the matrix
for n=1:10:500
    E_11=[];
    E_12=[];
    E_13=[];
    for j=1:2*n %running the iteration for 100 times for same size matrix
        % Generating random toeplitz matrix
        c=randn(n,1);
        r=randn(1,n);
        c(1)=r(1);
        T=toeplitz(c,r);
        G=randn(n,n);
    
    
        % Algo-1 circular decomposition , Taking only first circulant component
        R=zeros(n,1);
        R(1)=c(1);
        for i=2:n
            R(i)=(1/n)*(((i-1)*(r(n-i+2)))+((n-i+1)*(c(i))));   %O(n)
        end
        % R is the first row of first circulant component of toeplitz matrix
        %Generating Circulant matrix from row vector
        C_11 = gen_circ(R); % O(n^2) (Actually there is no need to find complete matrix only the first column is sufficient)
        M_11 = circulant_multiplication(C_11,G); %O(n^2 log n)
        M = T*G;
        norm_M = norm(M, 'fro');
        e_11=norm(abs(M_11 - M), 'fro');
        e_11=e_11/norm_M;
        E_11=[E_11,e_11];
        %Error for first circulant decomposition algorithm

    
        %Algo-2 converting n-Toeplitz to 2n-circulant matrix
        v=zeros(2*n,1); %O(2n)
        v(1:n)=T(1,:);
        for i = n+2 : 2*n
            v(i,1) = T((2*(n+1)-i),1);  %O(n)
        end
        v(n+1,1)=randn;
        C_12 = gen_circ(v);  % circulant matrix of 2n is generated from toeplitz matrix %O(n^2)
        M_12 = circulant_multiplication(C_12,G); %O(n^2 log 2n)
    
        e_12=norm(abs(M_12 - M), 'fro');
        e_12=e_12/norm_M;
        E_12=[E_12,e_12];   
        %Error for second circulant decomposition algorithm
   
    
        
    
        %Algo-3 Random Sampling
        c=n/2;
        M_13 = Randomsampling(T, G, c);
        e_13=norm(abs(M_13 - M), 'fro');
        e_13=e_13/norm_M;
        E_13=[E_13,e_13];      
        
        itr=[itr,j];
    end
    t_1= n^2 * log(n);
    T_1=[T_1,t_1];
    t_2= 4*n^2 * log(2*n);
    T_2=[T_2,t_2]; 
    t_3= c*n^2;
    T_3=[T_3,t_3]; 
    e_n1=mean(E_11);
    E_n1=[E_n1,e_n1];
    e_n2=mean(E_12);
    E_n2=[E_n2,e_n2];
    e_n3=mean(E_13);
    E_n3=[E_n3,e_n3];  
    t_4=n^3;
    T_4=[T_4,t_4];
    size=[size,n];
    disp(n);
end
subplot(2,1,1)
plot(size,E_n1);
hold on;
plot(size, E_n2);
hold on;
plot(size, E_n3);
hold off;
xlabel('size of matrix')
ylabel('Error')
legend('Circulant comp','circulant 2n','random sampling(c=n/2)');
title('Error in Matrix multiplication');

subplot(2,1,2)
plot(size,T_1);
hold on;
plot(size, T_2);
hold on;
plot(size, T_3);
hold on;
plot(size, T_4);
hold off;
xlabel('Size of matrix')
ylabel('Time')
legend('Circulant comp','circulant 2n','random sampling(c=n/2)','Navie Multiplication');
title('Error in Matrix multiplication');
disp('Mean error: ');
disp(mean(E_n1));
disp(mean(E_n2));
disp(mean(E_n3));
