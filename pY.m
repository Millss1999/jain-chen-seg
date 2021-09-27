function pY = pY(yi, y_hat, sigma)
  pY = (1/(sqrt(2*pi)*sigma)) * (1/exp((yi-y_hat)^2/sigma^2));