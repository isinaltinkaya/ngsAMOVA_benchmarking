
# read csv file into a data frame
# from first argument
args <- commandArgs(trailingOnly = TRUE)
d1 <- read.csv(args[1], header = FALSE, sep = "")
df2 <- read.csv(args[2], header = FALSE, sep = ",")
d2<-df2[4:12]


# get difference between d1 and d2
d3 <- d1 - d2

# print summary statistics to output file
# output is 3rd argument
sink(args[3])


print(d3)
summary(d3)

sink()
