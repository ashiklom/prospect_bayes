a <- rectGrob(gp = gpar(fill="red"))
b <- circleGrob()
cc <- linesGrob()
d <- nullGrob()

roww <- matrix(list(a,b,cc,d), nrow=1)
coll <- matrix(list(a,b,cc,d), ncol=1)
mat <- matrix(list(a,b,cc,d), nrow=2)

gtable_matrix("demo", mat, widths=c(1,2), heights=c(2,1))

