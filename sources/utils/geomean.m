function gmean = geomean(v)
  
gmean = exp(mean(log(abs(v))));

end
