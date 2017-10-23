f=dir(pattern="*.txt")
null = vector(mode="numeric")

for(i in 1:length(f)) {
    d = read.table(f[[i]])
    null = c(null, d$V1)
}
print(length(null))
print(quantile(null,probs = c(0,.025, .5, .975, 1)))
