function X = rbm_2_coo_3d(RBM)

ndof = size(RBM,1);
nnod = ndof/3;
X = zeros(nnod,3);

% Extract rotational modes
Ry_z =  RBM(1:3:end,5);
Rz_y = -RBM(1:3:end,6);
Rx_z = -RBM(2:3:end,4);
Rz_x =  RBM(2:3:end,6);
Rx_y =  RBM(3:3:end,4);
Ry_x = -RBM(3:3:end,5);

% Determine max rotational displacements
max_Ry_z = max(Ry_z);
max_Rz_y = max(Rz_y);
max_Rx_z = max(Rx_z);
max_Rz_x = max(Rz_x);
max_Rx_y = max(Rx_y);
max_Ry_x = max(Ry_x);

% Determine theta_x, theta_y and theta_z (up to a scaling factor)
fac_Rx_Rz = max_Rx_y / max_Rz_y;
fac_Rz_Ry = max_Rz_x / max_Ry_x;
fac_Rx_Ry = max_Rx_z / max_Ry_z;

% Arbitarily set Rx and find the other two
Rx = 1;
Ry = Rx / fac_Rx_Ry;
Rz = Rx / fac_Rx_Rz;

% Compute new coordinates
X(:,1) = Rz_x / Rz + Ry_x / Ry;
X(:,2) = Rz_y / Rz + Rx_y / Rx;
X(:,3) = Rx_z / Rx + Ry_z / Ry;

end
