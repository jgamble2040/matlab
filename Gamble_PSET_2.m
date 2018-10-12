%% LDLT Decomposition
% Construct an algorithm that turns A into LDLT

A = randn(10);
[n,n] = size(A);
A=triPosDef(n);
L=eye(n);


    
%inner loop generates L_{1-(n-1)}
    for i=1:(n-1)
        L(2,1)=-A(2,1)/A(1,1);
        Aprime=L*A*L';
        L(1,:)=[];
        L(:,1)=[];
        Aprime(1,:)=[];
        Aprime(:,1)=[];
        L(2,1)=-Aprime(2,1)/Aprime(1,1);
        M=1;
        M(2:10,2:10)=L(i,i);
    end
    end

L_final=M*M_1*M_2*M_3*M_4*M_5*M_6*M_7*M_8;


%Construct matrix L out of smaller matrices.

L(2,1)=-A(2,1)/A(1,1);
M=L;

L(2,1)=-A(2,1)/A(1,1);
Aprime=L*A*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_1=1;
M_1(2:10,2:10)=L;


Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_2=1;
M_2(2:9,2:9)=L;
M_2(2:10,2:10)=M_2;

Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_3=1;
M_3(2:8,2:8)=L;
M_3(2:9,2:9)=M_3;
M_3(2:10,2:10)=M_3;

Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_4=1;
M_4(2:7,2:7)=L;
M_4(2:8,2:8)=M_4;
M_4(2:9,2:9)=M_4;
M_4(2:10,2:10)=M_4;

Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_5=1;
M_5(2:6,2:6)=L;
M_5(2:7,2:7)=M_5;
M_5(2:8,2:8)=M_5;
M_5(2:9,2:9)=M_5;
M_5(2:10,2:10)=M_5;

Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_6=1;
M_6(2:5,2:5)=L;
M_6(2:6,2:6)=M_6;
M_6(2:7,2:7)=M_6;
M_6(2:8,2:8)=M_6;
M_6(2:9,2:9)=M_6;
M_6(2:10,2:10)=M_6;

Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_7=1;
M_7(2:4,2:4)=L;
M_7(2:5,2:5)=M_7;
M_7(2:6,2:6)=M_7;
M_7(2:7,2:7)=M_7;
M_7(2:8,2:8)=M_7;
M_7(2:9,2:9)=M_7;
M_7(2:10,2:10)=M_7;

Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_8=1;
M_8(2:3,2:3)=L;
M_8(2:4,2:4)=M_8;
M_8(2:5,2:5)=M_8;
M_8(2:6,2:6)=M_8;
M_8(2:7,2:7)=M_8;
M_8(2:8,2:8)=M_8;
M_8(2:9,2:9)=M_8;
M_8(2:10,2:10)=M_8;

Aprime=L*Aprime*L';
L(1,:)=[];
L(:,1)=[];
Aprime(1,:)=[];
Aprime(:,1)=[];
L(2,1)=-Aprime(2,1)/Aprime(1,1);
M_9=1;
M_9(2:2,2:2)=L;
M_9(2:3,2:3)=M_9;
M_9(2:4,2:4)=M_9;
M_9(2:5,2:5)=M_9;
M_9(2:6,2:6)=M_9;
M_9(2:7,2:7)=M_9;
M_9(2:8,2:8)=M_9;
M_9(2:9,2:9)=M_9;
M_9(2:10,2:10)=M_9;

%Construct D
D = L_final*A*L_final';

%Check results
answer = A - L_final*D*L_final';

%Create function

%Run
[L,D] = symmetric_tridiagonal_LU(A);
%% Question 2
m = 250;
% Strategy 1: uniform
t1 = linspace(0, 1, m)';
%% Strategy 2: clustered close to 0 and to 1, separately.
delta = 0.1;
ta = linspace(0, 0+delta, m/2)';
tb = linspace(1-delta, 1, m/2)';
t2 = [ta ; tb];
%% My work

 plot(t1,t2);
 
 %Construct two versions of matrix A. Then, QR factorize.
 %To articulate error bound of x, use 3.2 in lecture notes.
 %Is delta equal to 0.1 in worst-case perturbation of b?
 
%Build A_1 and A_2 
t1_0=t1.^0;
t1_1=t1;
t1_2=t1.^2;
A_1 = [t1_0 t1_1 t1_2];

t2_0=t2.^0;
t2_1=t2;
t2_2=t2.^2;

A_2 = [t2_0 t2_1 t2_2];

%Modified Gram-Schmidt Twice
[Q1, R1] = modified_gram_schmidt(A_1);
[Q2, R2] = modified_gram_schmidt(Q1);
Q_1 = Q2;
R_1 = R2*R1;

[Q1, R1] = modified_gram_schmidt(A_2);
[Q2, R2] = modified_gram_schmidt(Q1);
Q_2 = Q2;
R_2 = R2*R1;

Z=zeros(3,247);
R_1=[R_1 Z];
R_2=[R_2 Z];
%Compute error bounds on x.

perturbation = 10^(-3);

QR_1= norm(Q_1*R_1);
QR_1inv=norm(inv(R_1*Q_1));
xbound_1 = QR_1*QR_1inv*perturbation;

QR_2= norm(Q_2*R_2);
QR_2inv=norm(inv(R_2*Q_2));
xbound_2 = QR_2*QR_2inv*perturbation;

%Bound 2 means a smaller deviation in relative perturbance to x than bound
%1. Thus, Method 2 is better. 