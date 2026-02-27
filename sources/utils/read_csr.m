function A = read_csr(filename)

A = spconvert(dlmread(filename,'',1,0));

end
