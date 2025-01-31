clear all
clc


% dimension d, number of measurements N and Schmidt number to test k
d=4;
N=2;
k=2;

rho=MaxEntangled(d);
rho=rho*rho';

%%%MUB DEFINITION
P=GenPauli(1,0,d);
[a,b_x]=eig(P);
for i=1:d
    A{i,1}=a(:,i)*a(:,i)';
end
P=GenPauli(0,1,d);
[a,b_x]=eig(P);
for i=1:d
    A{i,2}=a(:,i)*a(:,i)';
end
P=GenPauli(1,1,d);
[a,b_x]=eig(P);
for i=1:d
    A{i,3}=a(:,i)*a(:,i)';
end
for j=2:16
    P=GenPauli(1,j,d);
    [a,b_x]=eig(P);
    for i = 1:d
        A{i,j+2}=a(:,i)*a(:,i)';
    end
end



%loading constraints obtained from steering_hierarchy_constraints_2_out_POVM.ipynb
load('steering_POVM_d4_N2_lvl_2.mat');
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
complex_M_pairs=complex_M_pairs+1;
complex_C_pairs=complex_C_pairs+1;


%loop over al SN from 2 to d-1
for k=2:d-1
    %clear yalmip memory at each iteration
    yalmip('clear')

    %definition of the variables:
    %reduced density matrix rho_B; assemblage sigma; complex matrices not_sigma; real and complex scalars c, cc; visibility v.  
    rho_B=sdpvar(d,d, 'hermitian', 'complex');
    sigma=sdpvar(d,d, sdp_mat_count-1, 'hermitian', 'complex');
    not_sigma=sdpvar(d,d, complex_mat_count-size(complex_M_pairs,1), 'full', 'complex');    
    c=sdpvar(coefs_count, 1, 'full', 'real');
    cc=sdpvar(complex_coefs_count-size(complex_C_pairs,1), 1, 'full', 'complex');
    v=sdpvar(1);
    
    % normalization is not required if we fix the assemblage
    % sigma3=sdpvar(d,d, N, 'hermitian', 'complex');

    F=[];

    % not required if we fix the assemblage
    % for i=1:N
    %     F=[F, sigma3(:,:,i)==rho_B - sigma(:,:,i) - sigma(:,:,i+N)];% - sigma(:,:,i+3*N) - sigma(:,:,i+4*N) - sigma(:,:,i+5*N) - sigma(:,:,i+6*N)
    %     F=[F, sigma3(:,:,i)>=0];
    % end

    %fixing the reduced desity matrix to I/d
    F=[F, rho_B==PartialTrace(rho*v+(1-v)*eye(d^2)/d^2,1)];
    % the assemblage is fixed as the result of MUBs measurements on the noisy maximally entangled state
    for i=1:d-1
        for j=1:N
            F=[F, sigma(:,:,(i-1)*N+j) == PartialTrace((rho*v+(1-v)*eye(d^2)/d^2)*Tensor(A{i,j}.', eye(d)),1)];
        end
    end

    %variables equal under full trace
    trace_pair=size(tracial_pairs);
    trace_pair=trace_pair(1);
    for i=1:trace_pair
        if(tracial_pairs(i,1)==0 && tracial_pairs(i,2)==0)
            F=[F, trace(sigma(:,:,tracial_pairs(i,3)))==trace(sigma(:,:,tracial_pairs(i,4)))];
        elseif(tracial_pairs(i,1)==0 && tracial_pairs(i,2)==1)
            if ismember(tracial_pairs(i,4), complex_M_pairs(:,1))
                target = find(complex_M_pairs(:,1)==tracial_pairs(i,4));
                F=[F, trace(sigma(:,:,tracial_pairs(i,3)))==trace(not_sigma(:,:,target))];
            elseif ismember(tracial_pairs(i,4), complex_M_pairs(:,2))
                target = find(complex_M_pairs(:,2)==tracial_pairs(i,4));
                F=[F, trace(sigma(:,:,tracial_pairs(i,3)))==trace(not_sigma(:,:,target))'];
            else
                F=[F, trace(sigma(:,:,tracial_paris(i,3)))==trace(not_sigma(:,:,tracial_pairs(i,4)-complex_mat_count))];
            end
        elseif(tracial_pairs(i,1)==1 && tracial_pairs(i,2)==0)
            if ismember(tracial_pairs(i,3), complex_M_pairs(:,1))
                target = find(complex_M_pairs(:,1)==tracial_pairs(i,3));
                F=[F, trace(not_sigma(:,:,target))==trace(sigma(:,:,tracial_pairs(i,4)))];
            elseif ismember(tracial_pairs(i,3), complex_M_pairs(:,2))
                target = find(complex_M_pairs(:,2)==tracial_pairs(i,3));
                F=[F, trace(not_sigma(:,:,target))'==trace(sigma(:,:,tracial_pairs(i,4)))];
            else
                F=[F, trace(not_sigma(:,:,tracial_paris(i,3)-complex_mat_count))==trace(sigma(:,:,tracial_pairs(i,4)))];
            end
        elseif(tracial_pairs(i,1)==1 && tracial_pairs(i,2)==1)
            if ismember(tracial_pairs(i,3), complex_M_pairs(:,1)) 
                target1 = find(complex_M_pairs(:,1)==tracial_pairs(i,3));
                if ismember(tracial_pairs(i,4), complex_M_pairs(:,1)) 
                    target2 = find(complex_M_pairs(:,1)==tracial_pairs(i,4));
                    F=[F, trace(not_sigma(:,:,target1))==trace(not_sigma(:,:,target2))];
                elseif ismember(tracial_pairs(i,4), complex_M_pairs(:,2))
                    target2 = find(complex_M_pairs(:,2)==tracial_pairs(i,4));
                    F=[F, trace(not_sigma(:,:,target1))==trace(not_sigma(:,:,target2))'];
                else
                    F=[F, trace(not_sigma(:,:,target1))==trace(not_sigma(:,:,tracial_pairs(i,4)-complex_mat_count))];
                end
            elseif ismember(tracial_pairs(i,3), complex_M_pairs(:,2))
                target1 = find(complex_M_pairs(:,2)==tracial_pairs(i,3));
                if ismember(tracial_pairs(i,4), complex_M_pairs(:,1)) 
                    target2 = find(complex_M_pairs(:,1)==tracial_pairs(i,4));
                    F=[F, trace(not_sigma(:,:,target1))'==trace(not_sigma(:,:,target2))];
                elseif ismember(tracial_pairs(i,4), complex_M_pairs(:,2))
                    target2 = find(complex_M_pairs(:,2)==tracial_pairs(i,4));
                    F=[F, trace(not_sigma(:,:,target1))==trace(not_sigma(:,:,target2))];
                else
                    F=[F, trace(not_sigma(:,:,target1))==trace(not_sigma(:,:,tracial_pairs(i,4)-complex_mat_count))];
                end
            elseif not(ismember(tracial_pairs(i,4), complex_M_pairs(:,1))) && not(ismember(tracial_pairs(i,4), complex_M_pairs(:,2)))
                F=[F, trace(not_sigma(:,:,tracial_pairs(i,3)-complex_mat_count))==trace(not_sigma(:,:,tracial_pairs(i,4)-complex_mat_count))];
            end
        end
    end
    
    % redundant if rho_B is fixed
    F=[F, trace(rho_B)==1];
    F=[F, rho_B>=0];
    F=[F, c(N*(d-1)+1)==k];

    % definition of Gamma
    Gamma=zeros(dim_gamma*d, dim_gamma*d);
    M=sparse(1, 1, 1, dim_gamma, dim_gamma);
    Gamma=Gamma+Tensor(M, k*d*eye(d));
    
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
    %Here, to reduce the number fo variables we exploit the fact that 
    %the elements in complex_C_pairs are equal up to the complex conjugate
    tmp=0;
    for i=1:complex_coefs_count
        index=find(complex_coefs(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma dim_gamma], index);
        if ismember(i, complex_C_pairs(:,1))
            target = find(complex_C_pairs(:,1)==i);
            for j=1:length(row_ind)
                if(row_ind(j)>col_ind(j))
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, cc(target)*eye(d));
                else
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, cc(target)'*eye(d));
                end
            end
        elseif ismember(i, complex_C_pairs(:,2))
            target = find(complex_C_pairs(:,2)==i);
            for j=1:length(row_ind)
                if(row_ind(j)>col_ind(j))
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, cc(target)'*eye(d));
                else
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, cc(target)*eye(d));
                end
            end
        else
            tmp=tmp+1;
            target=complex_C_count+tmp;
            for j=1:length(row_ind)
                if(row_ind(j)>col_ind(j))
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, cc(target)*eye(d));
                else
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, cc(target)'*eye(d));
                end
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
    %Here, to reduce the number fo variables we exploit the fact that the
    %elements in complex_M_pairs are equal up to the hermitian conjugate
    tmp=0;
    for i=1:complex_mat_count
        index=find(complex_matrices(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma dim_gamma], index);
        if ismember(i, complex_M_pairs(:,1))
            target = find(complex_M_pairs(:,1)==i);
            for j=1:length(row_ind)
                if(row_ind(j)>col_ind(j))
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,target), [d,d]));
                else
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,target)', [d,d]));
                end
            end
        elseif ismember(i, complex_M_pairs(:,2))
            target = find(complex_M_pairs(:,2)==i);
            for j=1:length(row_ind)
                if(row_ind(j)>col_ind(j))
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,target)', [d,d]));
                else
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,target), [d,d]));
                end
            end
        else
            tmp=tmp+1;
            target=complex_M_count+tmp;
            for j=1:length(row_ind)
                if(row_ind(j)>col_ind(j))
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,target), [d,d]));
                else
                    M=sparse(row_ind(j), col_ind(j), 1, dim_gamma, dim_gamma);
                    Gamma=Gamma+Tensor(M, reshape(not_sigma(:,:,target)', [d,d]));
                end
            end            
        end
    end

    % when the assemblage is fixed we optimise over v, if not the witness has to be defined
    % obj=trace(A{1,1}*sigma(:,:,1)+A{2,1}*sigma(:,:,4)+A{3,1}*sigma3(:,:,1)+...
    %           A{1,2}*sigma(:,:,2)+A{2,2}*sigma(:,:,5)+A{3,2}*sigma3(:,:,2)+...
    %           A{1,3}*sigma(:,:,3)+A{2,3}*sigma(:,:,6)+A{3,3}*sigma3(:,:,3));
    
    obj=v;
    
    F=[F, Gamma>=0];
    
    ops=sdpsettings('solver', 'mosek', 'mosek.MSK_IPAR_NUM_THREADS', 0); 
    
    diagnostic=solvesdp(F, -obj, ops)
    value(obj)
end