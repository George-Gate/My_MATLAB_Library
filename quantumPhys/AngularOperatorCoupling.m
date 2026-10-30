% GeorgeGate @ 2026-10-30
function ops = AngularOperatorCoupling(num_I, num_S, checkCommutator)
    %  return the addition of two spin operators
    %  [Inputs]
    %  num_I: spin 1 quantum number
    %  num_S: spin 2 quantum number
    %  checkCommutator: whether to check the commutation relations of angular operators in coupled representation
    %
    %  [Outputs]
    %  ops.rep_kron:    operators using the kron representation,    Hilbert space is H = H_I \otimes H_S
    %  ops.rep_coupled: operators using the coupled representation, Hilbert space is H = H_F1 \oplus H_F2 \oplus ...
    %  
    % [Eigen State Order (kron representation), I=3/2, S=1/2 for example]
    % | mI ,  mS >
    % |+3/2, +1/2>
    % |+3/2, -1/2>
    % |+1/2, +1/2>
    % |+1/2, -1/2>
    % |-1/2, +1/2>
    % |-1/2, -1/2>
    % |-3/2, +1/2>
    % |-3/2, -1/2>
    %
    % [Eigen State Order (coupled representation), I=3/2, S=1/2 for example]
    % |F, mF>
    % |2, +2>
    % |2, +1>
    % |2,  0>
    % |2, -1>
    % |2, -2>
    % |1, +1>
    % |1,  0>
    % |1, -1>

    if nargin<3
        checkCommutator = true;
    end
    
    %% construct operators in kron representation

    % angular operator in their own Hilbert space
    [Ix_single, Iy_single, Iz_single, Ip_single, Im_single] = AngularOperator(num_I);
    [Sx_single, Sy_single, Sz_single, Sp_single, Sm_single] = AngularOperator(num_S);
    IdenOp_I_single = eye(size(Ix_single));  % Identical operator in spin 1 space
    IdenOp_S_single = eye(size(Sx_single));  % Identical operator in spin 2 space

    % angular momentum operators in combined Hilbert space (kron representation)
    Ix = kron(Ix_single, IdenOp_S_single);
    Iy = kron(Iy_single, IdenOp_S_single);
    Iz = kron(Iz_single, IdenOp_S_single);
    Im = kron(Im_single, IdenOp_S_single);
    Ip = kron(Ip_single, IdenOp_S_single);
    
    Sx = kron(IdenOp_I_single, Sx_single);
    Sy = kron(IdenOp_I_single, Sy_single);
    Sz = kron(IdenOp_I_single, Sz_single);
    Sm = kron(IdenOp_I_single, Sm_single);
    Sp = kron(IdenOp_I_single, Sp_single);
    
    Fx = Sx + Ix;
    Fy = Sy + Iy;
    Fz = Sz + Iz;
    Fm = Sm + Im;
    Fp = Sp + Ip;
    
    I2 = Ix^2+Iy^2+Iz^2;
    S2 = Sx^2+Sy^2+Sz^2;
    F2 = Fx^2+Fy^2+Fz^2;
    
    % pack ops
    op_list = {'Fx', 'Fy', 'Fz', 'Fm', 'Fp', 'F2', 'Sx', 'Sy', 'Sz', 'Sm', 'Sp', 'S2', 'Ix', 'Iy', 'Iz', 'Ip', 'Im', 'I2'};

    rep_kron = struct();
    for ii = 1:length(op_list)
        rep_kron.(op_list{ii}) = eval(sprintf('%s', op_list{ii}));
    end


    %% transform to coupled representation
    % eigen states in coupled representation
    [V, D]  = eig(F2+Fz);
    % sort eigen values / eigen states
    [d,ind] = sort(diag(D), 'descend');
    D = D(ind,ind);
    V = V(:,ind);
    

    rep_coupled = struct();
    rep_coupled.eig_states_in_kron_basis = {};
    rep_coupled.eig_val_F2 = zeros(1, size(V,2));
    rep_coupled.eig_val_Fz = zeros(1, size(V,2));
    % Extract eigenvalues ​​and eigenvectors, and check accuracy
    for ii = 1:size(V,2)
        rep_coupled.eig_states_in_kron_basis{ii} = V(:,ii);
        rep_coupled.eig_val_F2(ii) = V(:,ii)' * F2 * V(:,ii);
        rep_coupled.eig_val_Fz(ii) = V(:,ii)' * Fz * V(:,ii);
        assert(norm( F2 * V(:,ii) - rep_coupled.eig_val_F2(ii)* V(:,ii) )<1e-10);
        assert(norm( Fz * V(:,ii) - rep_coupled.eig_val_Fz(ii)* V(:,ii) )<1e-10);
    end
    % basis transformation 
    for ii = 1:length(op_list)
        rep_coupled.(op_list{ii}) = eval(sprintf(' V'' * %s * V ', op_list{ii}));
    end

    % fix numerical residuals
    for ii = 1:length(op_list)
        tmp = rep_coupled.(op_list{ii});
        tmp(abs(tmp)<1e-12) = 0;
        rep_coupled.(op_list{ii}) = tmp;
    end
    

    %% check commutation relation of angular momentum operator
    if checkCommutator
        checkList = {rep_coupled.Fx, rep_coupled.Fy, 1i*rep_coupled.Fz, '[Fx,Fy]=i*Fz';
                     rep_coupled.Fy, rep_coupled.Fz, 1i*rep_coupled.Fx, '[Fy,Fz]=i*Fx';
                     rep_coupled.Fz, rep_coupled.Fx, 1i*rep_coupled.Fy, '[Fz,Fx]=i*Fy';
                     rep_coupled.F2, rep_coupled.Fz, 0, '[F^2,Fz]=0';
                     rep_coupled.Ix, rep_coupled.Iy, 1i*rep_coupled.Iz, '[Ix,Iy]=i*Iz';
                     rep_coupled.Iy, rep_coupled.Iz, 1i*rep_coupled.Ix, '[Iy,Iz]=i*Ix';
                     rep_coupled.Iz, rep_coupled.Ix, 1i*rep_coupled.Iy, '[Iz,Ix]=i*Iy';
                     rep_coupled.I2, rep_coupled.Iz, 0, '[I^2,Iz]=0';
                     rep_coupled.Sx, rep_coupled.Sy, 1i*rep_coupled.Sz, '[Sx,Sy]=i*Sz';
                     rep_coupled.Sy, rep_coupled.Sz, 1i*rep_coupled.Sx, '[Sy,Sz]=i*Sx';
                     rep_coupled.Sz, rep_coupled.Sx, 1i*rep_coupled.Sy, '[Sz,Sx]=i*Sy';
                     rep_coupled.S2, rep_coupled.Sz, 0, '[S^2,Sz]=0';};
        
        for ii = 1:size(checkList,1)
            A = checkList{ii,1};
            B = checkList{ii,2};
            C = checkList{ii,3};
            expName = checkList{ii,4};
            % err = max(abs( A*B-B*A - C ), [],'all');
            err = norm(A*B-B*A - C);
            if abs(err)>1e-12
                warning('The accuracy of relation "%s" in [.rep_coupled] is lower than the expected value. norm(err)=%g\n',expName, err);
            end
        end
    end

    ops.rep_kron = rep_kron;
    ops.rep_coupled = rep_coupled;
end