function [Lam_min,eigenvalues] = mPOD(N_s, filename, dt, Lam_min, M, M_max)

% mPOD Script
% This script calculates the proper orthogonal decomposition
% using the mPOD approach

% Inputs
% N_s - number of snapshots
% filename - Snapshot filenames are assumed to be of form filenameXXX.mat 
% where XXX is the snapshot number (from 1 to N_s).
% Snapshot files are assumed to contain a single variable named "data" 
% which is a 1D vector. 
% The inner product over the domain is assumed to be the dot product
% This can be modified for more complex data structures 
% dt - timestep of snapshots
% M - inital group size for mPOD decomposition
% M_max - maximum value of M before Lam_min gets increased
% Lam_min - cut-off below which modes are no longer kept


% Outputs
% Lam_min - final value of lambda_min used to respect M_max (see Rule 3)
% eigenvalues - final set of POD eigenvalues
% files named mode1.XX.txt where XX is the the mode number

% Internal Routines  (These can be modified for different data formats)
% These routines are at the bottom of this file
%
% phi = myload(filename) - function to load modes/snapshots 
%
% mysave(filename,phi) - function to save modes/snapshots 
%
% [result] = inner_product(snapshot1,snapshot2) - function to perform inner product 
%

% %%%%%%%%%%%%%%%%%%%%
% Some sanity checks %
%%%%%%%%%%%%%%%%%%%%%%
if (N_s <= 0)
    disp('Error N_s must be a positive integer');
    return;
end
if (dt < 0.0)
    disp('Error dt must be positive');
    return;
end
if (M <= 0) 
    disp('Error M must be a positive integer');
    return;
end
if (M_max < M)
    disp('Error M_max cannot be smaller than the initial group size');
    M_max = N_s;
    disp('Using M_max = N_s');
end
phi0 = myload([filename num2str(1)]);
max_eig_estimate = inner_product(phi0,phi0)*N_s*dt;
if (Lam_min < 1e-16*max_eig_estimate)
    disp('Warning: Lam_min is possibly requesting more than sixteen orders of magnitude in eigenvalues');
end
if (Lam_min > max_eig_estimate)
    disp('Warning: Lam_min is large, possibly no modes will be calculated');
end    

%%%%%%%%%%%%%
% BEGIN CODE
%%%%%%%%%%%%%
eigenvalues = dt*ones(N_s,1);   % 1st level weights snapshots by dt (eq. 5)
next_eigenvalues = zeros(size(eigenvalues)); % stores eigenvalues for next level
first = true;  % Boolean so that snapshots don't get deleted
N_g = ceil(N_s/M);  % Calculate initial number of sets based on M and N_s

% BEGIN LOOP FOR LEVELS OF MPOD
while (true)  % do while loop implemented with a break statement at end when N_g = 1
  disp(['N_g = ' num2str(N_g) ', M = ' num2str(M) ', Lam_min = ' num2str(Lam_min)]);
  modename = ['mode' num2str(N_g) '.'];  % output filename for this level's modes
  next_N_s = 0;     % counter of number of modes for next level
  
  % BEGIN LOOP OVER GROUPS
  for i=1:N_g       % Loop over each set (this loop could be parallelized...)
    I0 = (i-1)*M +1;    % index of beginning of set by eq. 15
    I1 = min(i*M,N_s);  % index of end of set avoiding integer errors
    Mi = I1-I0+1;       % M for this set, again avoiding integer errors
    
    % BEGIN LOOPS TO CREATE SNAPSHOT MATRIX FOR THIS GROUP
    C = zeros(Mi,Mi);   % this is the snapshot matrix (see eq. 17)
    for I = I0:I1                   % rows of snapshot matrix
      phiI = myload([filename num2str(I)]);  % load data
      row = I-I0+1;                 % local row number 
      C(row,row) = inner_product(phiI,phiI);  % diagonal entry of snapshot matrix
      for J = I+1:I1                  % off diagonal entries of snapshot matrix
        phiJ = myload([filename num2str(J)]); % load data for dot product
        col = J-I0+1;                 % local column number
        C(row,col) = inner_product(phiI,phiJ);      % this is a simple dot product 
        C(col,row) = C(row,col);      % C is symmetric 
      end
    end
    
    % SOLVE EIGENVALUE PROBLEM
    sqrtD = diag(sqrt(eigenvalues(I0:I1)));     % see eq. 18
    [V,Lambda] = eig(sqrtD*C*sqrtD');           % see eq. 18
    Lambda = diag(Lambda);                      % Extract diagonal of matrix
    Psi = sqrtD*V;                              % see eq. 18
    [Lambda,I] = sort(Lambda,'descend');        % sort eigenvalues 
    Psi = Psi(:,I);                             % reorder eigenvectors
    N_m = find(Lambda < Lam_min,1)-1;           % number of modes to keep
    
    % LOOP TO CALCULATE AND SAVE MODES FOR THIS GROUP
    for j = 1:N_m               % calculate new modes
      next_N_s = next_N_s+1;    % increment number of modes for next level
      phi = zeros(size(phiI));  % for summing snapshot/mode (eq. 4)
      for I = I0:I1
        phiI = myload([filename num2str(I)]);     % load snapshot/mode
        phi = phi +Psi(I-I0+1,j)*phiI;        % add to sum (eq. 4)
      end
      phi = phi / sqrt(inner_product(phi,phi));   % Renormalize mode
      mysave([modename num2str(next_N_s) '.txt'],phi); % mode for next level
      next_eigenvalues(next_N_s) = Lambda(j);       % eig for next level
    end
  end  % LOOP BACK FOR NEXT GROUP
  
  % SHIFTS FOR NEXT LEVEL OF mPOD
  N_s = next_N_s;     % shift N_s for next level
  eigenvalues = next_eigenvalues(1:N_s); % shift eigenvalues for next level
  if (~first) 
    delete([filename '*']);         % delete temporary mode files
  end
  filename = modename;              % change filename so next level loads modes
  first = false;
  
  % CHECK IF WE ARE DONE
  if (N_g == 1)  % if only one group we are done
    break;       % exit loop
  end
  
  % SET PARAMETERS FOR NEXT LEVEL
  N_g = ceil(min([N_g/2,N_s/M]));   % rules 1 & 2 in paper
  M = ceil(N_s/N_g);                % avoid noninteger values
  if (M > M_max)                    % rule 3
    M = M_max;                      % hit max size
    N_s = N_g*M_max;                % maximum allowed for next level
    [sorted_eigenvalues,I] = sort(eigenvalues,'descend');
    Lam_min = sorted_eigenvalues(N_s+1);  % new lambda_min
    I = I(1:N_s);                   % indices of modes to keep
    I = sort(I);                    % indices in ascending order
    for i=1:N_s                     % sequentially rename files and eigenvalues to keep 
       if (I(i) ~= i)
         movefile([filename num2str(I(i)) '.txt'],[filename num2str(i) '.txt']);
         eigenvalues(i) = eigenvalues(I(i));
       end
    end
  end
end  % LOOP FOR NEXT LEVEL OF mPOD

function [snapshot] = myload(filename)
    % function to load modes/snapshots (faster than matlab's load command)
    % Customize to load different formats here
    fileID = fopen([filename +'.txt'],'r');
    snapshot = fscanf(fileID,'%f');
    fclose(fileID);
return

function mysave(filename,phi)
    % function to save modes/snapshots 
    % Customize to save different formats here
    save(filename,'phi','-ascii'); % mode for next level
return

function [result] = inner_product(snapshot1,snapshot2)
    % function to perform inner product 
    % Customize for different snapshot data here
    result = snapshot1'*snapshot2;
return
    

    
