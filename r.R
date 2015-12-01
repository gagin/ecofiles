setwd("/Users/User/Desktop/ecofiles")

library(lubridate)

library(digest)
entries <- list()

splitUnlist <- function(s, p) {
        res <- unlist(strsplit(s, p))
        a <- 1
        if(res[[1]]=="") a <- 2
        res[a:length(res)]
}

countCharOccurrences <- function(s, char) {
        length(strsplit(s, char, fixed=TRUE)[[1]])-1
}

src <- file("2015-11-25 EROD pr .txt")
src.lines <- readLines(src)
dat <- list()
for(i in 1:length(src.lines)) {
        line.split <- splitUnlist(src.lines[i], "[[:blank:]]+")
        num.elements <- length(line.split)
        if(num.elements > 20) {
                entry <- digest(src.lines[i])
                entries[[length(entries)+1]] <- c(entry, src.lines[i])
        }
        if(num.elements == 14 & !grepl("T", line.split[1], fixed=TRUE)) {
                timer <- line.split[1]
                temp <- line.split[2]
                dat[[length(dat)+1]] <- c(entry, line.split)
        }
        if(num.elements == 12) {
                dat[[length(dat)+1]] <- c(entry, timer, temp, line.split)
        }
}
close(src)
data <- data.frame(matrix(unlist(dat), nrow=length(dat), byrow = TRUE),
                      stringsAsFactors = FALSE)
colnames(data) <- c("entry", "time", "temp", sapply(1:12,function(x)paste0("m",x)))

### Come back to this if data will actually  include hours
#data$time <- unlist(lapply(data$time, function(x) {
#        semicolumns.plusone <- length(strsplit(x, ":", fixed=TRUE)[[1]])
#        if(semicolumns.plusone == 2) return(as.POSIXct(x,format="%M:%S"))
#        if(semicolumns.plusone == 3) return(as.POSIXct(x,format="%H:%M:%S"))
#}))

data$time <- as.POSIXct(data$time, format="%M:%S")
for(i in 3:15) data[,i] <- as.numeric(data[,i])
