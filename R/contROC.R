"contROC" <-
function(gldstd, test1, test2 = NULL)
{
N = length(gldstd) 
vcomp1 = rep(0, length = N) 

for(i in 1:N)
{
vcomp1[i] = sum(gldstd[i] == gldstd[-i] | test1[i] == test1[-i])/2 +
sum(gldstd[i] < gldstd[-i] & test1[i] < test1[-i]) +
sum(gldstd[i] > gldstd[-i] & test1[i] > test1[-i])
}

vcomp1 = vcomp1/(N - 1)
acc1 = sum(vcomp1)/N
cov11 = sum((vcomp1 - acc1)^2)/(N/2 * (N/2 - 1))

if(is.null(test2) == TRUE)
{
acc1 = max(acc1, 1 - acc1)
return(list("Accuracy" = data.frame("Estimate" = acc1, "Standard Error" = sqrt(cov11))))
}
else
{
vcomp2 = rep(0, length = N) 

for(i in 1:N)
{
vcomp2[i] = sum(gldstd[i] == gldstd[-i] | test2[i] == test2[-i])/2 + 
sum(gldstd[i] < gldstd[-i] & test2[i] < test2[-i]) +
sum(gldstd[i] > gldstd[-i] & test2[i] > test2[-i])

}

vcomp2 = vcomp2/(N - 1)
acc2 = sum(vcomp2)/N
cov22 = sum((vcomp2 - acc2)^2)/(N/2 * (N/2 - 1))
cov12 = sum((vcomp1 - acc1) * (vcomp2 - acc2))/(N/2 * (N/2 - 1))

acc1 = max(acc1, 1 - acc1)
acc2 = max(acc2, 1 - acc2)
cov12 = abs(cov12)

z = (acc1 - acc2)/sqrt(cov11 + cov22 - 2 * cov12)
p = 2 * pnorm(abs(z), lower.tail = FALSE)

return(list("Accuracy" = data.frame("Estimates" = c(acc1, acc2), "Standard Errors" = 
sqrt(c(cov11, cov22))), "Covariance" = cov12, "Two-Sided Hypothesis Test" = 
data.frame("Z" = z, "p-value" = p)))
}
}

