"ordROC" <-
function(gldstd, test1, test2 = NULL, penalty = NULL)
{
num.st = length(levels(factor(gldstd)))
grp.st1 = split(test1, factor(gldstd))
num.pt = as.numeric(sapply(grp.st1, length))
tot.pt = length(gldstd)
p.pt = num.pt/tot.pt

if(is.null(penalty) == TRUE)
{
penalty = matrix(0, nrow = num.st, ncol = num.st)
penalty[col(penalty) > row(penalty)] = 1
}

acc1 = matrix(0, nrow = num.st, ncol = num.st)
var1 = matrix(0, nrow = num.st, ncol = num.st)
cov1 = array(0, c(num.st, num.st, num.st, num.st))
vcomp1 = array(0, c(num.st, num.st, tot.pt))    
    
for(k in 1:(num.st - 1))
{
        for(l in (k + 1):num.st)
        {
for(m in 1:num.pt[k])
{
vcomp1[k,l,m] = (sum(as.numeric(grp.st1[[k]][m]) == as.numeric(grp.st1[[l]]))/2 +
sum(as.numeric(grp.st1[[k]][m]) < as.numeric(grp.st1[[l]])))/num.pt[l]
}

for(n in 1:num.pt[l])
{
vcomp1[l,k,n] = (sum(as.numeric(grp.st1[[k]]) == as.numeric(grp.st1[[l]][n]))/2 +
sum(as.numeric(grp.st1[[k]]) < as.numeric(grp.st1[[l]][n])))/num.pt[k]
}

acc1[k,l] = sum(vcomp1[k,l,])/num.pt[k]
acc1[l,k] = acc1[k,l]

var1[k,l] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[k,l])^2)/(num.pt[k] * (num.pt[k] - 1)) +
sum((vcomp1[l,k,1:num.pt[l]] - acc1[l,k])^2)/(num.pt[l] * (num.pt[l] - 1))
var1[l,k] = var1[k,l]
}
}

for(k in 1:(num.st - 2))
{
for(l in (k + 1):(num.st - 1))
{
for(t in (l + 1):num.st)
{
cov1[k,l,k,t] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[k,l]) * (vcomp1[k,t,1:num.pt[k]] - acc1[k,t]))/
(num.pt[k] * (num.pt[k] - 1))
}               
}
}

if(num.st >= 3)
{
for(k in 3:num.st)
{
for(l in 1:(k - 2))
{
for(t in (l + 1):(k - 1))
{
                cov1[l,k,t,k] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[l,k]) * (vcomp1[k,t,1:num.pt[k]] - acc1[t,k]))/
(num.pt[k] * (num.pt[k] - 1))
}
}               
}
}

for(k in 2:(num.st - 1))
{
for(l in 1:(num.st - 2))
{
if(l < k)
{
for(t in (k + 1):num.st)
{
cov1[l,k,k,t] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[l,k]) * (vcomp1[k,t,1:num.pt[k]] - acc1[k,t]))/
(num.pt[k] * (num.pt[k] - 1))
}
}
}
}

denom = 0
for(o in 1:(num.st - 1))
{
for(s in (o + 1):num.st)
{
denom = denom + p.pt[o] * p.pt[s]
}
}

over.acc1 = 0
over.var1 = 0

for(o in 1:(num.st - 1))
{
for(s in (o + 1):num.st)
{
over.acc1 = over.acc1 + (p.pt[o] * p.pt[s]/denom) * penalty[o,s] * (1 - acc1[o,s])
over.var1 = over.var1 + (p.pt[o] * p.pt[s]/denom)^2 * penalty[o,s]^2 * var1[o,s]
}
}
over.acc1 = 1 - over.acc1

over.covar1 = 0

for(o in 1:num.st)
{
for(s in 1:num.st)
{
for(t in 1:num.st)
{
for(u in 1:num.st)
{
over.covar1 = over.covar1 + (p.pt[o] * p.pt[s]/denom) * (p.pt[t] * 
p.pt[u]/denom) * penalty[o,s] * penalty[t,u] * cov1[o,s,t,u] * 2
}    
}
}
}

over.se1 = sqrt(over.var1 + over.covar1)

acc11 = NULL
var11 = NULL
name.pair = NULL

rownames(penalty) = 1:num.st
colnames(penalty) = 1:num.st

for(i in 1:(num.st - 1))
{
acc11 = c(acc11, acc1[i, (i + 1):num.st])
var11 = c(var11, var1[i, (i + 1):num.st])
name.pair = c(name.pair, paste(i, " vs ", (i + 1):num.st, sep = ""))
}

if(is.null(test2) == TRUE)
{
acc11 = pmax(acc11, 1 - acc11)
over.acc1 = max(over.acc1, 1 - over.acc1)

return(list("Pairwise Accuracy" = data.frame("Pair" = name.pair, "Estimate" = acc11, 
"Standard Error" = sqrt(var11)), "Penalty Matrix" = penalty, 
"Overall Accuracy" = data.frame("Estimate" = over.acc1, "Standard Error" = over.se1)))
}

else
{
grp.st2 = split(test2, factor(gldstd))
acc2 = matrix(0, nrow = num.st, ncol = num.st)
var2 = matrix(0, nrow = num.st, ncol = num.st)
cov2 = array(0, c(num.st, num.st, num.st, num.st))
vcomp2 = array(0, c(num.st, num.st, tot.pt))
    
for(k in 1:(num.st - 1))
{
        for(l in (k + 1):num.st)
        {
for(m in 1:num.pt[k])
{
vcomp2[k,l,m] = (sum(as.numeric(grp.st2[[k]][m]) == as.numeric(grp.st2[[l]]))/2 +
sum(as.numeric(grp.st2[[k]][m]) < as.numeric(grp.st2[[l]])))/num.pt[l]
}

for(n in 1:num.pt[l])
{
vcomp2[l,k,n] = (sum(as.numeric(grp.st2[[k]]) == as.numeric(grp.st2[[l]][n]))/2 +
sum(as.numeric(grp.st2[[k]]) < as.numeric(grp.st2[[l]][n])))/num.pt[k]
}

acc2[k,l] = sum(vcomp2[k,l,])/num.pt[k]
acc2[l,k] = acc2[k,l]



var2[k,l] = sum((vcomp2[k,l,1:num.pt[k]] - acc2[k,l])^2)/(num.pt[k] * (num.pt[k] - 1)) +
sum((vcomp2[l,k,1:num.pt[l]] - acc2[l,k])^2)/(num.pt[l] * (num.pt[l] - 1))
var2[l,k] = var2[k,l]
}
}

for(k in 1:(num.st - 2))
{
for(l in (k + 1):(num.st - 1))
{
for(t in (l + 1):num.st)
{
cov2[k,l,k,t] = sum((vcomp2[k,l,1:num.pt[k]] - acc2[k,l]) * (vcomp2[k,t,1:num.pt[k]] - acc2[k,t]))/
(num.pt[k] * (num.pt[k] - 1))
}               
}
}

if(num.st >= 3)
{
for(k in 3:num.st)
{
for(l in 1:(k - 2))
{
for(t in (l + 1):(k - 1))
{
                cov2[l,k,t,k] = sum((vcomp2[k,l,1:num.pt[k]] - acc2[l,k]) * (vcomp2[k,t,1:num.pt[k]] - acc2[t,k]))/
(num.pt[k] * (num.pt[k] - 1))
}
}               
}
}

for(k in 2:(num.st - 1))
{
for(l in 1:(num.st - 2))
{
if(l < k)
{
for(t in (k + 1):num.st)
{
cov2[l,k,k,t] = sum((vcomp2[k,l,1:num.pt[k]] - acc2[l,k]) * (vcomp2[k,t,1:num.pt[k]] - acc2[k,t]))/
(num.pt[k] * (num.pt[k] - 1))
}
}
}
}


over.acc2 = 0
over.var2 = 0

for(o in 1:(num.st - 1))
{
for(s in (o + 1):num.st)
{
over.acc2 = over.acc2 + (p.pt[o] * p.pt[s]/denom) * penalty[o,s] * (1 - acc2[o,s])
over.var2 = over.var2 + (p.pt[o] * p.pt[s]/denom)^2 * penalty[o,s]^2 * var2[o,s]
}
}
over.acc2 = 1 - over.acc2

over.covar2 = 0

for(o in 1:num.st)
{
for(s in 1:num.st)
{
for(t in 1:num.st)
{
for(u in 1:num.st)
{
over.covar2 = over.covar2 + (p.pt[o] * p.pt[s]/denom) * (p.pt[t] * 
p.pt[u]/denom) * penalty[o,s] * penalty[t,u] * cov2[o,s,t,u] * 2
}    
}
}
}

over.se2 = sqrt(over.var2 + over.covar2)

cov12 = array(0, c(num.st, num.st, num.st, num.st))
for(k in 1:(num.st - 1))
{
        for(l in (k + 1):num.st)
        {
cov12[k,l,k,l] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[k,l]) * (vcomp2[k,l,1:num.pt[k]] - acc2[k,l]))/
(num.pt[k] * (num.pt[k] - 1)) + sum((vcomp1[l,k,1:num.pt[l]] - acc1[l,k]) * 
(vcomp2[l,k,1:num.pt[l]] - acc2[l,k]))/(num.pt[l] * (num.pt[l] - 1))
}
}

for(k in 1:(num.st - 2))
{
for(l in (k + 1):num.st)
{
for(t in (k + 1):num.st)
{
if(t != l)
{
                cov12[k,l,k,t] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[k,l]) * 
(vcomp2[k,t,1:num.pt[k]] - acc2[k,t]))/(num.pt[k] * (num.pt[k] - 1))
}               
}
}
}

if(num.st >= 3)
{
for(k in 3:num.st)
{
for(l in 1:(k - 1))
{
for(t in 1:(k - 1))
{
if(l != t)
{
cov12[l,k,t,k] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[l,k]) * 
(vcomp2[k,t,1:num.pt[k]] - acc2[t,k]))/(num.pt[k] * (num.pt[k] - 1))
}               
}
}               
}
}

for(k in 2:(num.st - 1))
{
for(l in 1:(num.st - 2))
{
if(l < k)
{
for(t in (k + 1):num.st)
{
cov12[l,k,k,t] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[l,k]) * 
(vcomp2[k,t,1:num.pt[k]] - acc2[k,t]))/(num.pt[k] * (num.pt[k] - 1))
}
}
}
}

for(k in 2:(num.st - 1))
{
for(t in 1:(num.st - 2))
{
if(t < k)
{
for(l in (k + 1):num.st)
{
cov12[k,l,t,k] = sum((vcomp1[k,l,1:num.pt[k]] - acc1[l,k]) * 
(vcomp2[k,t,1:num.pt[k]] - acc2[k,t]))/(num.pt[k] * (num.pt[k] - 1))
}
}
}
}

over.cov12 = 0

for(o in 1:num.st)
{
for(s in 1:num.st)
{
for(t in 1:num.st)
{
for(u in 1:num.st)
{
over.cov12 = over.cov12 + (p.pt[o] * p.pt[s]/denom) * (p.pt[t] * p.pt[u]/denom) * 
penalty[o,s] * penalty[t,u] * cov12[o,s,t,u] 
}    
}
}
}
over.acc1 = max(over.acc1, 1 - over.acc1)
over.acc2 = max(over.acc2, 1 - over.acc2)
over.cov12 = abs(over.cov12)

z = (over.acc1 - over.acc2)/sqrt(over.se1^2 + over.se2^2 - 2 * over.cov12)
p = 2 * pnorm(abs(z), lower.tail = FALSE)

acc22 = NULL
var22 = NULL

for(i in 1:(num.st - 1))
{
acc22 = c(acc22, acc2[i, (i + 1):num.st])
var22 = c(var22, var2[i, (i + 1):num.st])
}

acc11 = pmax(acc11, 1 - acc11)
acc22 = pmax(acc22, 1 - acc22)

return(list("Pairwise Accuracy for Test 1" = data.frame("Pair" = name.pair, "Estimate" = acc11, 
"Standard Error" = sqrt(var11)), "Pairwise Accuracy for Test 2" = data.frame("Pair" = name.pair, 
"Estimate" = acc22, "Standard Error" = sqrt(var22)), "Penalty Matrix" = penalty, "Overall Accuracy" = 
data.frame("Estimate" = c(over.acc1, over.acc2), "Standard Error" = c(over.se1, over.se2)), 
"Covariance" = over.cov12, "Two-Sided Hypothesis Test" = data.frame("Z" = z, "p-value" = p)))
}
}

