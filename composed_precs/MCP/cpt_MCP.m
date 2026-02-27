function MCP_prec = cpt_MCP(param,AA,TV0,BB,CC)

AMG_prec = cpt_aspAMG(param,AA,TV0);

nthread = 8;
nstep = 30;
step_size = 1;
epsilon = 1.e-3;
G11 = afsai_cpp(AA,nthread,nstep,step_size,epsilon);

HH = G11*BB;
SS = CC - HH'*HH;

%ZZ = CC - BB'*inv(AA)*BB;
%save('SS');
%CC

G22 = afsai_cpp(-SS,nthread,nstep,step_size,epsilon);

MCP_prec.AMG_prec = AMG_prec;
MCP_prec.G22 = G22;

end
