setwd("/Users/User/Desktop/helga")

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
                if(countCharOccurrences(line.split[1], ":")==1)
                        timer <- ms(line.split[1])
                if(countCharOccurrences(line.split[1], ":")==2)
                        timer <- hms(line.split[1])
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
#data$time <- sapply(data$time, function(x) strsplit(x, ":", fixed=TRUE)[[1]]
        