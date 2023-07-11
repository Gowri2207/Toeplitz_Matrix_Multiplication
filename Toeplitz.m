e_circ_comp=[];
e_cir_2n=[];
e_random_samp_100=[];
e_random_samp_10=[];
e_random_samp_50=[];

length=[];
TC1=[];
TC2=[];
TC3_100=[];
TC3_10=[];
TC3_50=[];
nor=[];
min_val=0;
max_val=9;
for n = 3:250
    %generating a general matrix
    A=randi([min_val,max_val],n);
    
    c = randi([0,9], 1,n); %generate a random vetor of length n with elements between 0 1nd 10 for col in toeplitz
    r = randi([0,9], 1,n); %generate a random vector of length n with elements between 0 and 10 for row in toeplitz
    c(1)=randi([0,9],1);
    r(1)=c(1); %fixing the first elements of two vectors
    T = toeplitz(c,r); % generate toeplitz matrix of size n using vectors c and r
    R_True=T*A; % Exact solution of T*A
    
    % General matrix x first circulant component of toeplitz matrix
    
    R_1=Algo1(T,A);
    Error_1=(R_1-R_True);
    E_1=norm(Error_1,'fro');
    tc1=(n^3)+((n^2)*(log(n)));
    e_circ_comp=[e_circ_comp,E_1];
    TC1=[TC1,tc1];
    %fprintf('Error from algo 1 %d\n',E_1);
    
    % General matrix x circulant matrix of size 2n generated from toeplitz matrix
    R_2=Algo_2n(T,A);
    Error_2=(R_2 - R_True);
    E_2=norm(Error_2, 'fro');
    tc2=(n^2)*log(2*n);
    %fprintf('Error from algo 2 %d',E_2);
    e_cir_2n=[e_cir_2n,E_2];
    TC2=[TC2,tc2];

    %From random sampling of the matrices
    
    R_3=Randomsampling(T,A,10);
    Error_3=(R_3-R_True);
    E_3=norm(Error_3,'fro');
    e_random_samp_10=[e_random_samp_10,E_3];
    tc3_10=10 * n^2;
    TC3_10=[TC3_10,tc3_10];

    R_4=Randomsampling(T,A,50);
    Error_3=(R_4-R_True);
    E_4=norm(Error_3,'fro');
    e_random_samp_50=[e_random_samp_50,E_4];
    tc3_50=50 * n^2;
    TC3_50=[TC3_50,tc3_50];

    R_5=Randomsampling(T,A,100);
    Error_5=(R_5-R_True);
    E_5=norm(Error_5,'fro');
    e_random_samp_100=[e_random_samp_100,E_5];
    tc3_100=100 * n^2;
    TC3_100=[TC3_100,tc3_100];
    length=[length , n];
    nor=[nor,n^3]

end
fprintf('Time complexity of Toeplitz x General when first circulant component of Toeplitz matrix is considered\n is O(n^3+n^2log n');

fprintf('\n');
fprintf('Time complexity in multiplication of Toeplitz x General by generating 2n circulant matrix from Toeplitz matrix is O(n^2 log2n)\n')
%fprintf('%.2f ',e_cir_2n);
fprintf('\n');
fprintf('Time complexity in multiplication random sampling of matrix is O(c * n^2)\n')
%fprintf('%.2f ',e_random_samp_10);
fprintf('\n');
subplot(2,1,1)
plot(length, e_circ_comp);
hold on;
plot(length, e_cir_2n);
hold on;
plot(length, e_random_samp_100);
hold off;
title('Length vs Error');
xlabel('Length of matrix');
ylabel('Error');
legend('Circulant comp','circulant 2n','random samplingc100');
subplot(2,1,2)
plot(length, TC1);
hold on;
plot(length, TC2);
hold on;
plot(length, TC3_100);
hold on;
plot(length, nor);
hold off;
title('Length vs Time Complexity');
xlabel('Length of matrix');
ylabel('Time Complexity');
legend('Circulant comp','circulant 2n','random samplingc100','Navie Mult');
