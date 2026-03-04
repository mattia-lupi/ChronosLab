% Script to prepare binaries for aspAMG

fileIN = fopen('Prepare_Bin.fnames','r');
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_MATRIX  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_TV0     = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_RHS     = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_BIN     = D{1};
fclose(fileIN);

% Read the system matrix
A = read_csr(file_MATRIX);

% Read the initial test space
TV0 = dlmread(file_TV0,'',1,0);

rhs = load(file_RHS);

% Write the binary file
save(file_BIN)
