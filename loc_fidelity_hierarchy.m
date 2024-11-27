clear all
clc


%define dimension d and number of measurements N, the initial fidelity is 1/d

d=3;
N=3;
f_0=1/d;


%define the steps of the for loop and the delta to increase the fidelity at each step

N_count=10;
delta_f=(1-1/d)/N_count;
data_=[];

%load the data saved from loc_fidelity_hierarchy_constraints_2.ipynb
load('loc_fidelity_d3_N3_lvl_2.mat');

for opt_count=1:N_count
    f=f_0+(delta_f*(opt_count-1));
    sdp_mat_count=double(var_num(2));
    dim_gamma_rho_y=size(sdp_matrices_y);
    dim_gamma_rho_y=dim_gamma_rho_y(2);
    complex_mat_count=double(var_num(4));
    coefs_count=double(var_num(1));
    complex_coefs_count=double(var_num(3));
    complex_C_count=double(length(complex_C_pairs));
    complex_M_count=double(length(complex_M_pairs));

    %variables definition:
    %reduced density matrix rho_B, assemblages sigma, normalization of the assemblages sigma3, complex matrix variables sigma_3;
    %real valued scalars c, complex valued scalars c3
    
    rho_B=sdpvar(d,d, 'hermitian', 'complex');
    sigma=sdpvar(d,d, sdp_mat_count-1, 'hermitian', 'complex');
    sigma3=sdpvar(d,d,N,'hermitian','complex');
    not_sigma=sdpvar(d,d, complex_mat_count, 'full', 'complex');
    c=sdpvar(coefs_count, 1, 'full', 'real');
    cc=sdpvar(complex_coefs_count, 1, 'full', 'complex');
     
    var_count=1;
    
    F=[];
    
    %RHO_B IS A DENSITY MATRIX
    F=[F, trace(rho_B)==1];
    F=[F, rho_B>=0];

    %variables related by the complex conjugate operation
    for i=1:complex_C_count
        F=[F, cc(complex_C_pairs(i,1)+1)==cc(complex_C_pairs(i,2)+1)'];
    end

    %variables related by the hermitian conjugate operation
    for i=1:complex_M_count
        F=[F, not_sigma(:,:,complex_M_pairs(i,1)+1)==not_sigma(:,:,complex_M_pairs(i,2)+1)'];
    end

    %variables equal under trace operation
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

    %NO SIGNALING ON SIGMA, THE CONDITION HAS TO BE CHANGED TO RESPECT THE NUMBER OF OUTCOMES
    for i=1:N
        F=[F, sigma3(:,:,i) == rho_B - sigma(:,:,i) - sigma(:,:,i+N)];%
        F=[F, sigma3(:,:,i) >= 0];
    end

    % %NO SIGNALING ON SIGMA
    % for i=1:N
    %     F=[F, sigma3(:,:,i) == rho_B - sigma(:,:,i)];
    %     F=[F, sigma3(:,:,i) >= 0];
    %     F=[F, c3(i)==f*d-c(i+1)];
    % end

    %definition of the block matrices
    Gamma_rho=zeros(dim_gamma_rho_y*d, dim_gamma_rho_y*d);
    Gamma_y=zeros(dim_gamma_rho_y*d, dim_gamma_rho_y*d);
    
    
    
    %LOCALISED MATRIX RHO, SDP MATRICES
    index=find(sdp_matrices_rho(1,:,:)==1);
    [row_ind, col_ind] = ind2sub([dim_gamma_rho_y dim_gamma_rho_y], index);
    for j=1:length(row_ind)
        M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
        Gamma_rho=Gamma_rho-Tensor(M, rho_B);
    end
    var_count=var_count+1;
    for i = 2:length(rho_indices)
        F=[F, sigma(:,:,i-1)>=0];
        index=find(sdp_matrices_rho(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma_rho_y dim_gamma_rho_y], index);
        for j=1:length(row_ind)
            M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
            Gamma_rho=Gamma_rho-Tensor(M, reshape(sigma(:,:,rho_indices(i)-1), [d,d]));
        end
        var_count=var_count+1;
    end
    
    %LOCALISED MATRIX RHO, COMPLEX MATRICES
    for i = 1:length(complex_rho_indices)
        index=find(complex_matrices_rho(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma_rho_y dim_gamma_rho_y], index);
        for j=1:length(row_ind)
            if(row_ind(j)>col_ind(j))
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
                Gamma_rho=Gamma_rho-Tensor(M, reshape(not_sigma(:,:,complex_rho_indices(i)), [d,d]));
            else
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
                Gamma_rho=Gamma_rho-Tensor(M, reshape(not_sigma(:,:,complex_rho_indices(i))', [d,d]));
            end
        end
    end
    
    
    %LOCALISED MATRIX Y, SDP MATRICES
    for i = 1:length(hermitian_y_indices)
        index=find(sdp_matrices_y(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma_rho_y dim_gamma_rho_y], index);
        for j=1:length(row_ind)
            M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
            Gamma_y=Gamma_y+Tensor(M, reshape(sigma(:,:,hermitian_y_indices(i)-1), [d,d]));
        end
    end
    
    %LOCALISED MATRIX Y, COMPLEX MATRICES
    for i = 1:length(complex_mat_y_indices)
        index=find(complex_matrices_y(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma_rho_y dim_gamma_rho_y], index);
        for j=1:length(row_ind)
            if(row_ind(j)>col_ind(j))
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
                Gamma_y=Gamma_y+Tensor(M, reshape(not_sigma(:,:,complex_mat_y_indices(i)), [d,d]));
            else
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
                Gamma_y=Gamma_y+Tensor(M, reshape(not_sigma(:,:,complex_mat_y_indices(i))', [d,d]));
            end
        end
    end
    
    %LOCALISED MATRIX Y, REAL COEFS
    for i=1:length(real_y_indices)
        index=find(coefs_y(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma_rho_y dim_gamma_rho_y], index);
        for j=1:length(row_ind)
            M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
            Gamma_y=Gamma_y+Tensor(M, c(real_y_indices(i))*eye(d));
        end
    end
    
    %LOCALISED MATRIX Y, COMPLEX COEFS
    for i=1:length(complex_y_indices)
        index=find(complex_coefs_y(i,:,:)==1);
        [row_ind, col_ind] = ind2sub([dim_gamma_rho_y dim_gamma_rho_y], index);
        for j=1:length(row_ind)
            if(row_ind(j)>col_ind(j))
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
                Gamma_y=Gamma_y+Tensor(M, cc(complex_y_indices(i))*eye(d));
            else
                M=sparse(row_ind(j), col_ind(j), 1, dim_gamma_rho_y, dim_gamma_rho_y);
                Gamma_y=Gamma_y+Tensor(M, cc(complex_y_indices(i))'*eye(d));
            end
        end
    end

    %fidelity constraint
    F=[F, c(1)==f*d];

    
    Gamma_rho_y=Gamma_y+Gamma_rho;
    
    %semidefinite positiveness constraint
    F=[F, Gamma_rho<=0];
    F=[F, Gamma_rho_y>=0];
    

    %definition of the objective
    P=GenPauli(1,0,d);
    [a,b_x]=eig(P);
    for i=1:d
        B{i,1}=a(:,i)*a(:,i)';
    end
    P=GenPauli(0,1,d);
    [a,b_x]=eig(P);
    for i=1:d
        B{i,2}=a(:,i)*a(:,i)';
    end
    P=GenPauli(1,1,d);
    [a,b_x]=eig(P);
    for i=1:d
        B{i,3}=a(:,i)*a(:,i)';
    end
    P=GenPauli(1,2,d);
    [a,b_x]=eig(P);
    for i=1:d
        B{i,4}=a(:,i)*a(:,i)';
    end

    % obj=trace(B{1,1}*sigma(:,:,1))+trace(B{2,1}*sigma(:,:,1+N))+trace(B{3,1}*sigma3(:,:,1))+ ...
    %     trace(B{1,2}*sigma(:,:,2))+trace(B{2,2}*sigma(:,:,2+N))+trace(B{3,2}*sigma3(:,:,2))+ ...
    %     trace(B{1,3}*sigma(:,:,3))+trace(B{2,3}*sigma(:,:,3+N))+trace(B{3,3}*sigma3(:,:,3))+ ...
    %     trace(B{1,4}*sigma(:,:,4))+trace(B{2,4}*sigma(:,:,4+N))+trace(B{3,4}*sigma3(:,:,4));

    obj=trace(B{1,1}*sigma(:,:,1))+trace(B{2,1}*sigma(:,:,1+N))+trace(B{3,1}*sigma3(:,:,1))+ ...
        trace(B{1,2}*sigma(:,:,2))+trace(B{2,2}*sigma(:,:,2+N))+trace(B{3,2}*sigma3(:,:,2))+ ...
        trace(B{1,3}*sigma(:,:,3))+trace(B{2,3}*sigma(:,:,3+N))+trace(B{3,3}*sigma3(:,:,3));

    % obj=trace(B{1,1}*sigma(:,:,1))+trace(B{2,1}*sigma3(:,:,1))+ ...
    %     trace(B{1,2}*sigma(:,:,2))+trace(B{2,2}*sigma3(:,:,2))+ ...
    %     trace(B{1,3}*sigma(:,:,3))+trace(B{2,3}*sigma3(:,:,3));

    % obj=trace(B{1,1}*sigma(:,:,1))+trace(B{2,1}*sigma(:,:,3))+trace(B{3,1}*sigma3(:,:,1))+ ...
    %     trace(B{1,2}*sigma(:,:,2))+trace(B{2,2}*sigma(:,:,4))+trace(B{3,2}*sigma3(:,:,2));

    ops=sdpsettings('solver', 'mosek');
    diagnostic=solvesdp(F, - obj, ops)
    value(obj)
    data_=[data_; f, obj];
end

