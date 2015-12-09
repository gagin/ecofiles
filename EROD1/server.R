library(shiny)

library(ggplot2)
#library(dplyr)
#library(tidyr)
####### Helper functions for calibration charts

get.formula <- function(model) {
        if(length(model$coefficients) == 2) {
                sign.char <- "+"
                if(model$coefficients[1] < 0) sign.char <- "-"
                res <- paste0("y = ",
                              sprintf("%.4f", model$coefficients[2],4),
                              "x ",
                              sign.char,
                              " ",
                              sprintf("%.4f", abs(coefficients(model)[1]))
                )
        }
        if(length(model$coefficients) == 1) {
                res <- paste0("y = ",
                              sprintf("%.4f", model$coefficients[1],4),
                              "x"
                )    
        }
        res
}

get.r.squared <- function(model) {
        sprintf("R squared = %.4f", summary(model)$r.squared)
}

point.two <- function(x, shift=0.2) shift*(range(x)[2]-range(x)[1])


shinyServer(function(input, output) {
        #output$chart <- renderPlot({
         react<-reactive({
                 
                inFile1 <- input$file1
        
        inFile2 <- input$file2
        
        if (is.null(inFile1))
               return(NULL)
        
        if (is.null(inFile2))
               return(NULL)
        
        
        #filename <- "2015-11-26 EROD pr .txt"
        filename <- inFile1$datapath
        #namebase <- sub(".txt$", "", filename, fixed=FALSE)
        #pattern.name <- paste0(namebase, ".csv")
        # TODO make it read Excel files directly
        #if(!file.exists(pattern.name)) pattern.name <- "default.csv"
        pattern.name <- inFile2$datapath
        pattern <- read.csv(pattern.name, header = FALSE, na.strings="",
                            stringsAsFactors = FALSE)
        block.length <- dim(pattern)[1]
        block.width  <- dim(pattern)[2]
        measurements.count <- sum(!is.na(pattern))
        pattern.vector <- as.vector(as.matrix(pattern))
        
        timers <- seq(0, 30, 2)
        drop.lines <- 5
        # TODO Is BLOCK= 3 comment meaningful?
        src <- file(filename)
        src.lines <- readLines(src)
        close(src)
        
        data.length <- length(timers) * measurements.count
        datata <- data.frame(timer=numeric(data.length),
                           probe=character(data.length),
                           fluor=numeric(data.length),
                           stringsAsFactors = FALSE)
        
        for(block in 0:(length(timers)-1)) {
                timer <- timers[block+1]
                data.lines <- measurements.count*block + 1:measurements.count
                datata[data.lines, "timer"] <- timer
                datata[data.lines, "probe"] <- pattern.vector[!is.na(pattern.vector)]
                datata[data.lines, "fluor"] <-
                        as.numeric(
                                as.vector(
                                        matrix(
                                                unlist(
                                                        strsplit(
                                                                src.lines[drop.lines + block*(block.length+1) + 1:block.length],
                                                                "\\t"
                                                        )
                                                ),
                                                ncol = block.width + 3,
                                                byrow = TRUE
                                        )[, 3:(block.width+2)]
                                )[!is.na(pattern.vector)]
                        )
        }
        
        ### Protein calibration
        src.split <- strsplit(src.lines, "\\t")
        ## TODO get rid of lists here
        proteins.list <- list()
        for( i in 153:160) {
                if(src.split[[i]][2] != "") temp <- src.split[[i]][2]
                if(src.split[[i]][2] == "") src.split[[i]][2] <- temp
                proteins.list[[length(proteins.list)+1]] <- src.split[[i]]
        }
        proteins <- data.frame(matrix(as.numeric(unlist(proteins.list)),
                                      nrow=length(proteins.list),
                                      byrow = TRUE))
        proteins.x <- c(proteins[1:6, 13],proteins[1:6 ,14])
        proteins.y <- rep(c(0, 0.1, 0.3, 0.75, 1.1, 1.4), times=2)
        
        
        
        ##### Without averaging
        protein1.calib <- data.frame(fluor=proteins.x,
                                     conc=proteins.y)
        protein1.model <- lm(conc ~ fluor, data=protein1.calib)
        
        print(
                ggplot(aes(fluor, conc),
                       data=protein1.calib) +
                        geom_point() +
                        stat_smooth(method="lm", slope=1, show_guide=TRUE) +
                        ggtitle(paste(
                                "Protein Calibration Curve w/o averaging\n",
                                filename)) +
                        annotate("text",
                                 x=min(proteins.x)+point.two(proteins.x),
                                 y=max(proteins.y)-point.two(proteins.y),
                                 label=paste(
                                         get.formula(protein1.model),
                                         get.r.squared(protein1.model),
                                         sep="\n")
                        )+
                        xlab("Fluorescence, unitless") +
                        ylab("Protein conc., mg/ml")
        )
        
        #####
        
        # Replicate Excel instead - overwrite previous two lines
        proteins.x <- rowMeans(cbind(proteins[1:6, 13],proteins[1:6 ,14]))
        proteins.y <- c(0, 0.1, 0.3, 0.75, 1.1, 1.4)
        
        protein.calib <- data.frame(fluor=proteins.x,
                                    conc=proteins.y)
        protein.model <- lm(conc ~ fluor, data=protein.calib)
        
        ### Resorufin calibration
        
        rs.list <- list()
        for( i in 141:146) {
                if(grepl("Temp", src.split[[i]][2], fixed=TRUE)) next
                if(src.split[[i]][2] != "") temp <- src.split[[i]][2]
                if(src.split[[i]][2] == "") src.split[[i]][2] <- temp
                rs.list[[length(rs.list)+1]] <- src.split[[i]]
        }
        rs <- data.frame(matrix(as.numeric(unlist(rs.list)),
                                nrow=length(rs.list),
                                byrow = TRUE))
        
        rs.x <- c(rs[1:6, 11],rs[1:6 ,12])
        rs.y <- rep(c(0, 0.1, 0.2, 0.3, 0.4, 0.5), times=2)
        
        rs.calib <- data.frame(fluor=rs.x, conc=rs.y)
        rs.model <- lm(conc ~ fluor-1, data=rs.calib)
        
        ######## Protein Calibration Chart
        
#         print(
                p1<-ggplot(aes(fluor, conc),
                       data=protein.calib) +
                        geom_point() +
                        stat_smooth(method="lm", slope=1, show_guide=TRUE) +
                        ggtitle("Protein Calibration Curve") +
                        annotate("text",
                                 x=min(proteins.x)+point.two(proteins.x),
                                 y=max(proteins.y)-point.two(proteins.y),
                                 label=paste(
                                         get.formula(protein.model),
                                         get.r.squared(protein.model),
                                         sep="\n")
                        )+
                        xlab("Fluorescence, unitless") +
                        ylab("Protein conc., mg/ml")
#         )
        
        ####### Resorufin Calibration Chart
        
#         print(
                p2<-ggplot(aes(fluor, conc),
                       data=rs.calib) +
                        geom_point() +
                        stat_smooth(method="lm", slope=1, show_guide=TRUE) +
                        ggtitle("Resorufin Calibration Curve") +
                        annotate("text",
                                 x=min(rs.x)+point.two(rs.x),
                                 y=max(rs.y)-point.two(rs.y),
                                 label=paste(
                                         get.formula(rs.model),
                                         get.r.squared(rs.model),
                                         sep="\n")
                        )+
                        xlab("Fluorescence, unitless") +
                        ylab("Resorufin conc, uM")
#         )
        
        ###### Calculating slopes
        
        datata$protein.abs <- rep(as.vector(as.matrix(proteins[,3:14]))[!is.na(pattern.vector)], times = length(timers))
        datata$rs.conc <- predict(rs.model, datata[, "fluor", drop=FALSE])
        datata$protein.conc <- predict(protein.model, data.frame(fluor=datata$protein.abs))
        datata$resor.prot <- datata$rs.conc/datata$protein.conc
        
        # Note what rounding did in Excel
        #> 5159*0.0000962-0.1656
        #[1] 0.3306958
        #> 5159*0.0001-0.16
        #[1] 0.3559
        
        samples <- unique(pattern.vector)
        samples <- samples[!is.na(samples)]
        
        
        models <- list()
        submodels <- list()
        legendary.i <- character(length(samples))
        legendary <- character()
        
        per.probe <- nrow(datata[datata$probe==samples[1] & datata$timer==0,])
        submask <- rep(FALSE, per.probe)
        
        for(i in 1:length(samples)) {
                templist <- list()
                for(j in 1:per.probe) {
                        
                        mask <- submask
                        mask[j] <- TRUE
                        mask <- rep(mask, length(timers))
                        subdata <- datata[datata$probe==samples[i], ]
                        subdata <- subdata[mask,]
                        templist[[j]] <- lm(resor.prot ~ timer, subdata)
                }  
                submodels[[i]] <- templist
                
                models[[i]] <- lm(resor.prot ~ timer, datata[datata$probe==samples[i], ])
                legendary.i[i] <- get.formula(models[[i]])
                legendary <- paste0(legendary, 
                                    "\n",
                                    samples[i],
                                    ": ",
                                    legendary.i[i]
                )
        }
        
        slopes <- data.frame(probe = samples)
        
        for(j in 1:per.probe) {
                subslope <- numeric(length(samples))
                for(i in 1:length(samples)) {
                        subslope[i] <- coefficients(submodels[[i]][[j]])["timer"]
                }
                slopes[, paste("Fluor", j)] <- subslope
        }
        total.slope <- numeric(length(samples))
        for(i in 1:length(samples)) {
                total.slope[i] <- coefficients(models[[i]])["timer"]
        }
        slopes[, "All"] <- total.slope
        
        #write.csv(x = slopes, file = paste0(namebase, " - Slopes.csv"), row.names = FALSE)
        
        ## Common scale
        
        #print(
               p3<-ggplot(aes(timer, resor.prot, color=probe), data=datata) +
                       ggtitle("Resorufin/Protein per time, single scale") +
                       geom_point() +
                       guides(colour=FALSE) +
                       facet_wrap(~probe) +
                       geom_text(data=data.frame(probe=samples,l=legendary.i),
                                 size=3,
                                 aes(x=10,#min(datata$timer)+point.two(datata$timer, 0.5),
                                     y=0.5,#max(datata$resor.prot)-point.two(datata$resor.prot, 0.2),
                                     label = l),
                                 inherit.aes=FALSE, parse=FALSE)
         
               p4<- ggplot(aes(timer, resor.prot, color=probe), data=datata) +
                       ggtitle("Resorufin/Protein per time, separate scales") +                
                       geom_point() +
                       guides(colour=FALSE) +
                       facet_wrap(~probe, scales="free", ncol=4) +
                       geom_text(data=data.frame(probe=samples,l=legendary.i),
                                 size=3,
                                 aes(x=10,#min(datata$timer)+point.two(datata$timer, 0.35),
                                     y=0.5,#min(datata$resor.prot)+point.two(datata$resor.prot, 0.1),
                                     label = l),
                                 inherit.aes=FALSE, parse=FALSE)
                
        #))
        #print(p)
        list(p1=p1,p2=p2,p3=p3,p4=p4,slopes=slopes)
         })
         
        output$contents <- renderTable({react()$slopes}, digits=5)
        output$p1 <- renderPlot({react()$p1})
        output$p2 <- renderPlot({react()$p2})
        output$p3 <- renderPlot({react()$p3}, height = 1000)
        output$p4 <- renderPlot({react()$p4}, height = 1000)
        #output$contents <- renderTable(slopes)
})