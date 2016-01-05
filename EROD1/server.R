library(shiny)

library(ggplot2)
####### Legend and formatting functions for calibration charts

# We want to draw slope and intercept as a nice formula, Excel-like
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

# Keep 4 digits after dot
get.r.squared <- function(model) {
        sprintf("R squared = %.4f", summary(model)$r.squared)
}

point.two <- function(x, shift=0.2) shift*(range(x)[2]-range(x)[1])



# Console debug
# setwd("/Users/User/Desktop/ecofiles"); filename <- "2015-11-26 EROD pr .txt"; pattern.name <- "default.csv"
# setwd("/Users/User/Dropbox/eco"); filename <- "2016-01-04-EROD-Pr-a.txt"; pattern.name <- "Mappnig-2016-01-04-a.csv"
# input <- list(); input$proteins <- "0, 0.1, 0.3, 0.75, 1.1, 1.4"; input$rs <- "0, 0.1, 0.2, 0.3, 0.4, 0.5"; input$proteins.scale <- "0.172"; input$rs.scale <- "150"; input$inGroup <- "2"
# Deployment
# library(checkpoint); checkpoint("2015-12-01"); library(rsconnect); deployApp('EROD1',appName="EROD1")

shinyServer(function(input, output) {
        #output$chart <- renderPlot({
         react<-reactive({
                 
                inFile1 <- input$file1
        
        inFile2 <- input$file2
        
        # Don't do anything until both input files uploaded
        if (is.null(inFile1))
               return(NULL)
        
        if (is.null(inFile2))
               return(NULL)
        
        
        filename <- inFile1$datapath
        pattern.name <- inFile2$datapath
        #### Run manual debung from below this line
        pattern <- read.csv(pattern.name, header = FALSE, na.strings="",
                            stringsAsFactors = FALSE)
        block.length <- dim(pattern)[1]
        block.width  <- dim(pattern)[2]
        measurements.count <- sum(!is.na(pattern))
        pattern.vector <- as.vector(as.matrix(pattern))
        
        # These parameters are considered fixed for the machine used as data
        # source and experiment setting. Still, where possible we keep
        # calculation dynamic, so future adjustments can potentially be done
        # via change in this single place
        timers <- seq(0, 30, 2)
        drop.lines <- 5
        src <- file(filename)
        src.lines <- readLines(src)
        close(src)

        # Initialize new structure for the data.
        # Because data isn't just squarish, but also some cells are
        # arbitrary have different meaning (used for calibration),
        # we don't just use tidyr, but cycle and remember cell positions
        # in separate two columns
        data.length <- length(timers) * measurements.count
        dat <- data.frame(rnum=integer(data.length),
                             cnum=integer(data.length),
                             timer=numeric(data.length),
                           probe=character(data.length),
                           fluor=numeric(data.length),
                           stringsAsFactors = FALSE)
        
        # We extract all probes with non-empty names in mapping file
        for(block in 0:(length(timers)-1)) {
                timer <- timers[block+1]
                nonNA.counter <- 1
                zone <- matrix(unlist(strsplit(
                        src.lines[drop.lines + block*(block.length+1) +
                                          1:block.length],
                        "\\t")),
                        ncol = block.width + 3,
                        byrow = TRUE)[, 3:(block.width+2)]
                data.lines <- measurements.count*block + 1:measurements.count
                dat[data.lines, "timer"] <- timer
                for(co in 1:block.width) {
                        for (ro in 1:block.length) {
                                if(!is.na(pattern[ro, co])) {
                                        nonNA.line <- data.lines[nonNA.counter]
                                        dat[nonNA.line, "rnum"] <- ro
                                        dat[nonNA.line, "cnum"] <- co
                                        dat[nonNA.line, "probe"] <- pattern[ro, co]
                                        dat[nonNA.line, "fluor"] <- as.numeric(zone[ro, co])
                                        nonNA.counter <- nonNA.counter + 1
                                }
                        }
                }
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
        
        #####
        #proteins.x <- rowMeans(cbind(proteins[1:6, 13],proteins[1:6 ,14]))
        
        proteins.y <- eval(parse(text=paste0("c(",input$proteins,")")))
        proteins.y.length <- length(proteins.y)
        proteins.y <- rep(proteins.y, times=2)
        proteins.y <- proteins.y * as.numeric(input$proteins.scale)
        
        proteins.x <- c(proteins[1:proteins.y.length, 13],
                        proteins[1:proteins.y.length ,14])
        
        protein.calib <- data.frame(fluor=proteins.x,
                                    conc=proteins.y)
        protein.model <- lm(conc ~ fluor, data=protein.calib)
        
        ### Resorufin calibration
        # This copying of temperature to every row isn't actually needed,
        # it's an artifact from a previous version where table was made unified
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
        
        
        rs.y <- eval(parse(text=paste0("c(",input$rs,")")))
        rs.y.length <- length(rs.y)
        rs.y <- rep(rs.y, times=2)
        rs.y <- rs.y * as.numeric(input$rs.scale)
        
        rs.x <- c(rs[1:rs.y.length , 11],
                  rs[1:rs.y.length  ,12])
        
        rs.calib <- data.frame(fluor=rs.x, conc=rs.y)
        # Notice that -1 means it's regression through the origin
        rs.model <- lm(conc ~ fluor, data=rs.calib)
        
        ######## Protein Calibration Chart
        
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

        ####### Resorufin Calibration Chart
        
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

        ###### Calculating slopes
        
        dat$protein.abs <- rep(as.vector(as.matrix(proteins[, 3:14])
                                            )[!is.na(pattern.vector)],
                                  times = length(timers))
        dat$rs.conc <- predict(rs.model, dat[, "fluor", drop=FALSE])
        dat$protein.conc <- predict(protein.model,
                                       data.frame(fluor=dat$protein.abs))
        dat$resor.prot <- dat$rs.conc/dat$protein.conc
        
        
        samples <- unique(pattern.vector)
        samples <- samples[!is.na(samples)]
        

        # Calculate slopes per probe name, based on all points.
        # It's the most meaningful way, although isn't what was set as the task.
        # The version that does it for each pair separately, as was asked,
        # is farther below.
        models <- list()
        legendary.i <- character(length(samples))
        legendary <- character()
        slopes <- data.frame(Probe = samples)
        for(i in 1:length(samples)) {
                models[[i]] <- lm(resor.prot ~ timer,
                                  dat[dat$probe==samples[i], ])
                legendary.i[i] <- get.formula(models[[i]])
                legendary <- paste0(legendary, 
                                    "\n",
                                    samples[i],
                                    ": ",
                                    legendary.i[i]
                )
                slopes[i, "All points"] <- coefficients(models[[i]])["timer"]
        }
        
        
        # Calculate slopes for every cell in pattern separately
        
        slopes.grid <- matrix(ncol=block.width, nrow=block.length)
        for(co in 1:block.width)
                for (ro in 1:block.length) {
                        if(!is.na(pattern[ro, co]))
                                slopes.grid[ro, co] <- coefficients(
                                        lm(resor.prot ~ timer,
                                           dat[dat$rnum==ro &
                                                          dat$cnum==co,]
                                        ))["timer"]
                }
                                
        # Now let's make it narrow                                                          
        columns.per.probe <- as.integer(input$inGroup)
        slopes.lines <- measurements.count/columns.per.probe
        slopes.narrow <- data.frame(Probe=character(slopes.lines),
                             stringsAsFactors = FALSE)
        for(i in 1:columns.per.probe)
                slopes.narrow[, paste("Fluor", i)] <- numeric(slopes.lines)
        slopes.set <- 1
        for(co in seq(1, block.width, columns.per.probe))
                for(ro in 1:block.length)
                        if(!is.na(pattern[ro, co])) {
                                slopes.narrow[slopes.set,
                                              "Probe"] <- pattern[ro, co]
                                for(i in 1:columns.per.probe)
                                        slopes.narrow[slopes.set, paste("Fluor",
                                                                        i)] <- 
                                                slopes.grid[ro, co + i - 1]
                                slopes.set <- slopes.set + 1      
                        }
        
        ## Common scale
        
               p3<-ggplot(aes(timer, resor.prot, color=probe), data=dat) +
                       ggtitle("Resorufin/Protein per time, single scale") +
                       geom_point() +
                       guides(colour=FALSE) +
                       facet_wrap(~probe) +
                       geom_text(data=data.frame(probe=samples,l=legendary.i),
                                 size=3,
                                 # Shiny breaks on out-of-reactive call for some
                                 # reason, so we just use approximate position
                                 aes(x=10,
                                     y=0.5,
                                     label = l),
                                 inherit.aes=FALSE, parse=FALSE)
         
               p4<- ggplot(aes(timer, resor.prot, color=probe), data=dat) +
                       ggtitle("Resorufin/Protein per time, separate scales") +                
                       geom_point() +
                       guides(colour=FALSE) +
                       facet_wrap(~probe, scales="free", ncol=4) +
                       geom_text(data=data.frame(probe=samples,l=legendary.i),
                                 size=3,
                                 aes(x=10,
                                     y=0.5,
                                     label = l),
                                 inherit.aes=FALSE, parse=FALSE)
                

        # Return results
               list(p1=p1,
                    p2=p2,
                    p3=p3,
                    p4=p4,
                    slopes.grid=slopes.grid,
                    slopes.narrow=slopes.narrow,
                    slopes=slopes)
         })

# Calculations finished, let's pass resulting data and charts back to ui.R
         
        output$slopes <- renderTable({react()$slopes},
                                     digits=5, include.rownames=FALSE)
        output$slopes.narrow <- renderTable({react()$slopes.narrow},
                                            digits=5, include.rownames=FALSE)
        output$slopes.grid <- renderTable({react()$slopes.grid},
                                          digits=5, include.rownames=FALSE)
        output$p1 <- renderPlot({react()$p1})
        output$p2 <- renderPlot({react()$p2})
        output$p3 <- renderPlot({react()$p3}, height = 1000)
        output$p4 <- renderPlot({react()$p4}, height = 1000)
})