setwd("/Users/User/Desktop/ecofiles")

filename <- "2015-11-25 EROD pr .txt"
src <- file(filename)
src.lines <- readLines(src)
data.list <- list()
proteins.list <- list()
src.split <- strsplit(src.lines,"\\t")
for(i in 1:length(src.lines)) {
        if(src.lines[i] == "~End") {
                timer <- NA
                temp <- NA
                sample.number <- NA
        }
        if(length(src.split[[i]]) != 15) next
        first.split <- strsplit(src.split[[i]][1], ":0", fixed=TRUE)[[1]]
        if(length(first.split) > 1) {
                timer <- first.split[1]
                temp <- src.split[[i]][2]
                sample.number <- 1
        }
        if(!is.na(timer)) {
                if(src.split[[i]][1] == "") sample.number <- sample.number + 1
                # repopulate missing fields
                src.split[[i]][1] <- timer
                src.split[[i]][2] <- temp
                data.list[[length(data.list)+1]] <- c(sample.number,src.split[[i]])
                next
        }
        # Handle last group of lines
        if(grepl("Temp", src.split[[i]][2], fixed=TRUE)) next
        if(src.split[[i]][2] != "") temp <- src.split[[i]][2]
        if(src.split[[i]][2] == "") src.split[[i]][2] <- temp
        proteins.list[[length(proteins.list)+1]] <- src.split[[i]]
}
close(src)
data <- data.frame(matrix(as.numeric(unlist(data.list)),
                          nrow=length(dat),
                          byrow = TRUE))
# drop last empty column
data <- data[, 1:15]
colnames(data) <- c("sample", "time", "temp", sapply(1:12,function(x)paste0("m",x)))

proteins <- data.frame(matrix(as.numeric(unlist(proteins.list)),
                          nrow=length(proteins.list),
                          byrow = TRUE))
# drop last empty column
proteins <- proteins[, 2:14]

