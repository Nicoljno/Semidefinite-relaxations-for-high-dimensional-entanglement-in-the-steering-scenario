clear all
clc


% dimension d, number of measurements N and Schmidt number to test k

d=3;
N=9;
k=2;
data_=[];
rho=MaxEntangled(d);
rho=rho*rho';

% definition of the sic POVMs, change the generator v when needed
v=[+.78867513459481288225457439025097872782e+0+.00000000000000000000000000000000000000e+0i;
+.28867513459481288225457439025097872782e+0-.50000000000000000000000000000000000000e+0i;
+.10566243270259355887271280487451063608e+0-.18301270189221932338186158537646809173e+0i;];


% v=[+.48571221409126403909152153176812197109e+0+.00000000000000000000000000000000000000e+0i;
% +.60043369656069688700611847041568366744e+0-.44989636690811813902417022753091663501e+0i;
% +.00000000000000000000000000000000000000e+0-.20118858648686589293456281596678826706e+0i;
% -.39924511007383099407155565444889540038e+0-.35815847183145900067351304237205336087e-1i;];
v=v*v';
for i=1:d
    for j=1:d
        tmp=GenPauli(i-1,j-1,d);
        A{1,(i-1)*d+j}=tmp*v*tmp';
        A{2,(i-1)*d+j}=eye(d)-A{1,(i-1)*d+j};
    end
end


%loading constraints obtained from steering_hierarchy_constraints_2_out_simple.ipynb
load('sic_d3_lvl_2.mat');
sdp_mat_count=size(sdp_matrices);
dim_gamma=sdp_mat_count(2);
sdp_mat_count=sdp_mat_count(1);
complex_mat_count=size(complex_matrices);
complex_mat_count=complex_mat_count(1);
coefs_count=size(coefs);
coefs_count=coefs_count(1);
complex_coefs_count=size(complex_coefs);
complex_coefs_count=complex_coefs_count(1);
complex_C_count=double(length(complex_C_pairs));
complex_M_count=double(length(complex_M_pairs));

%loop over al SN from 2 to d
for k=2:d

    %clear yalmip memory at each iteration
    yalmip('clear')

    %definition of the variables:
    %reduced density matrix rho_B; assemblage sigma; complex matrices not_sigma; real and complex scalars c, cc; visibility v.  

    rho_B=sdpvar(d,d, 'hermitian', 'complex');
    sigma=sdpvar(d,d, sdp_mat_count-1, 'hermitian', 'complex');
    not_sigma=sdpvar(d,d, complex_mat_count, 'full', 'complex');
    c=sdpvar(coefs_count, 1, 'full', 'real');
    cc=sdpvar(complex_coefs_count, 1, 'full', 'complex');
    v=sdpvar(1);
    
    % normalization is not needed if we fix the assemblage
    % sigma3=sdpvar(d,d, N, 'hermitian', 'complex');
    
    
    F=[];
    F=[F, v<=1];
    
    %fixing the reduced desity matrix to I/d
    F=[F, rho_B==PartialTrace(rho*v+(1-v)*eye(d^2)/d^2,1)];

    %variables related by complex conjugate operation
    for i=1:complex_C_count
        F=[F, cc(complex_C_pairs(i,1)+1)==cc(complex_C_pairs(i,2)+1)'];
    end
    %variables related by hermitian conjugate operation
    for i=1:complex_M_count
        F=[F, not_sigma(:,:,complex_M_pairs(i,1)+1)==not_sigma(:,:,complex_M_pairs(i,2)+1)'];
    end

    % not needed if we fix the assemblage
    % for i=1:N
    %     F=[F, sigma3(:,:,i)==rho_B - sigma(:,:,i)];
    %     F=[F, sigma3(:,:,i)>=0];
    % end

    % the assemblage is fixed as coming from MUB measurements on the maximally entangled state
    for j=1:N
        F=[F, sigma(:,:,j) == PartialTrace((rho*v+(1-v)*eye(d^2)/d^2)*Tensor(A{1,j}.', eye(d)),1)];
    end

    %variables equal under full trace
    trace_pair=size(tracial_pairs);
    trace_pair=trace_pair(1);
    for i=1:trace_pair
        if(tracial_pairs(i,1)==0 && tracial_pairs(i,2)==0)
            F=[F, trace(sigma(:,:,tracial_pairs(i,3)))==trace(sigma(:,:,tracial_pairs(i,4)))];
        end
        if(tracial_pairs(i,1)==0 && tracial_pairs(i,2)==1)
            F=[F, trace(sigma(:,:,tracial_pairs(i,3)))==trace(not_sigma(:,:,tracial_pairs(i,4)))];
        end
        if(tracial_pairs(i,1)==1 && tracial_pairs(i,1)==0)
            F=[F, trace(not_sigma(:,:,tracial_pairs(i,3)))==trace(sigma(:,:,tracial_pairs(i,4)))];
        end
        if(tracial_pairs(i,1)==1 && tracial_pairs(i,2)==1)
            F=[F, trace(not_sigma(:,:,tracial_pairs(i,3)))==trace(not_sigma(:,:,tracial_pairs(i,4)))];
        end
    end
    
    % redundant if rho_B is fixed
    F=[F, trace(rho_B)==1];
    F=[F, rho_B>=0];

    % definition of Gamma
    Gamma=zeros(dim_gamma*d, dim_gamma*d);
    M=sparse(1, 1, 1, dim_gamma, dim_gamma);
    Gamma=Gamma+Tensor(M, k*eye(d));
    
    %TRACE_A(A_axA_by ⊗ I)=c*I; TRACE_A(A_axA_byA_ax ⊗ I)=c*I; ...
    for i=1:coefs_count
        index=find(coefs(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma dim_gamma], index);
        for j=1:length(row_ind)
            M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
            Gamma=Gamma+Tensor(M, c(i)*eye(d));
        end
    end
    
    %TRACE_A(A1A2A3 ⊗ I)=complex*I ...
    for i=1:complex_coefs_count
        index=find(complex_coefs(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma dim_gamma], index);
        for j=1:length(row_ind)
            if(row_ind(j)>col_ind(j))
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                Gamma=Gamma+Tensor(M, cc(i)*eye(d));
            else
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                Gamma=Gamma+Tensor(M, cc(i)'*eye(d));
            end
        end
    end
    
    %UNKNOWN HERMITIAN SDP MATRICES
    index=find(sdp_matrices(1,:,:)==1);
    [row_ind, col_ind] = ind2sub([dim_gamma dim_gamma], index);
    for j=1:length(row_ind)
        M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
        Gamma=Gamma+Tensor(M, rho_B);
    end
    for i=2:sdp_mat_count
        index=find(sdp_matrices(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma dim_gamma], index);
        F=[F, sigma(:,:,i-1)>=0];
        for j=1:length(row_ind)
            M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
            Gamma=Gamma+Tensor(M, reshape(sigma(:,:,i-1), [d,d]));
        end
    end
    
    %COMPLEX VALUED MATRICES
    for i=1:complex_mat_count
        index=find(complex_matrices(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma dim_gamma], index);
        for j=1:length(row_ind)
            if(row_ind(j)>col_ind(j))
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,i), [d,d]));
            else
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,i)', [d,d]));
            end
        end
    end

    % when the assemblage is fixed we optimise over v, if not the witness has to be defined
    
    obj=v;

    % obj=trace(A{1,1}*sigma(:,:,1)+A{2,1}*sigma3(:,:,1)+...
    %           A{1,2}*sigma(:,:,2)+A{2,2}*sigma3(:,:,2)+...
    %           A{1,3}*sigma(:,:,3)+A{2,3}*sigma3(:,:,3)+...
    %           A{1,4}*sigma(:,:,4)+A{2,4}*sigma3(:,:,4)+...
    %           A{1,5}*sigma(:,:,5)+A{2,5}*sigma3(:,:,5)+...
    %           A{1,6}*sigma(:,:,6)+A{2,6}*sigma3(:,:,6)+...
    %           A{1,7}*sigma(:,:,7)+A{2,7}*sigma3(:,:,7)+...
    %           A{1,8}*sigma(:,:,8)+A{2,8}*sigma3(:,:,8)+...
    %           A{1,9}*sigma(:,:,9)+A{2,9}*sigma3(:,:,9)+...
    %           A{1,10}*sigma(:,:,10)+A{2,10}*sigma3(:,:,10)+...
    %           A{1,11}*sigma(:,:,11)+A{2,11}*sigma3(:,:,11)+...
    %           A{1,12}*sigma(:,:,12)+A{2,12}*sigma3(:,:,12)+...
    %           A{1,13}*sigma(:,:,13)+A{2,13}*sigma3(:,:,13)+...
    %           A{1,14}*sigma(:,:,14)+A{2,14}*sigma3(:,:,14)+...
    %           A{1,15}*sigma(:,:,15)+A{2,15}*sigma3(:,:,15)+...
    %           A{1,16}*sigma(:,:,16)+A{2,16}*sigma3(:,:,16))/16;
    
    F=[F, Gamma>=0];
    
    ops=sdpsettings('solver', 'mosek');%, 'mosek.MSK_IPAR_NUM_THREADS', 0); 
    
    diagnostic=solvesdp(F, -obj, ops)
    value(obj)
    data_=[data_, value(obj)];
end
