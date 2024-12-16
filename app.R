library(shiny)
#library(bslib)
library(shinyWidgets)
library(ggplot2)
library(anomalous)

# Workaround for Chromium Issue 468227
# Need this to properly download the csv file
# this bug and workaround is only for shinylive, you do not need it in your regular app
downloadButton <- function(...) {
  tag <- shiny::downloadButton(...)
  tag$attribs$download <- NULL
  tag
}

ui <- navbarPage(
    "CAPAView",
    tabPanel("Introduction",
             p("First tab content.")),
    tabPanel(title = "Methodology",
             p("Second tab content.")),
    tabPanel(title = "Changing Classes",
             p("Second tab content.")),
    tabPanel(title = "Profiles",
             sidebarLayout(
                 sidebarPanel(
                     selectInput(
                         inputId = "profile_row",
                         label = h3("Profile:"),
                         choices =  "",
                         selectize = FALSE,
                         selected = ""
                     ),
                     sliderInput(
                         inputId = "profile_train",
                         label = h3("Training Window:"),
                         min=3,max=20,value=14,
                         width="100%"
                     ),
                     h3("Anomaly Selection:"),
                     sliderInput("profile_betaP", 
                                 label = "Point Anomaly Threshold:",
                                 min = 2, max = 6, step=0.1,value = 3,
                                 width="100%"),
                     sliderTextInput("profile_nc", 
                                     label = "Number of Collective Anomalies:",
                                     "",
                                     width="100%"),
                     helpText("Select anomalies by specifying the number of deviations away from the mean for a point anomaly. Then select the number of collective anomalies from the options"),
                     h3("Data:"),
                     helpText("The provided data set can be download using the link below. You own data in the same format can be uploaded for analysis."),
                     downloadButton("profile_downloadData", "Download Example Data",width="50%"),
                     fileInput("profile_fileInput", "Upload your data", multiple = FALSE, accept = "csv", width = "50%", buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                     width=3
                 ),
                 mainPanel(
                     fluidRow({ plotOutput(outputId = "profile_plot") }),
                     fluidRow(
                         h3("Anomalies"),
                         tableOutput('profile_summaryTable')
                     ),
                     width=9
                 )
             )
             ),
    tabPanel(title = "Count Profiles",
             p("Second tab content."))
)

server <- function(input, output) {

    ## set reactive values
    DATA <- reactiveValues(categorical = NULL,
                           profile = as.matrix( read.csv("example_data.csv",row.names=1) ),
                           profile_anom = NULL,
                           count = NULL)


    ## ###########################################
    ## for profile example

    ## load the profile data
    profile_get_data <- observe({
        DATA$profile <- as.matrix( read.csv( input$profile_fileInput$datapath,row.names=1 ) )
    })
    bindEvent(profile_get_data,input$profile_fileInput)

    profile_change_options <- observe({
        updateSelectInput(session = getDefaultReactiveDomain(),
                          inputId="profile_row",
                          choices= rownames(DATA$profile)[(input$profile_train+1):nrow(DATA$profile)]
                          )
    })
    bindEvent(profile_change_options,input$profile_train,DATA$profile)

    
    profile_run_anom <- observe({
        print("running detection")
        if( is.null(DATA$profile) ){ return(NULL) }

        D <- DATA$profile
        betaP <- input$profile_betaP

        ii <- which(rownames(D)==input$profile_row)
        if( length(ii)==0){return(NULL)}
        idx <- ii-(1:input$profile_train)
        if( any(idx<=0) ){return(NULL)}
        x <- D[ii,]
        X <- D[idx,]
        m <- apply(X,2,median,na.rm=T)
        s <- pmax(1e-6,apply(X,2,mad,na.rm=T))
        nx <- length(x)
        
        betaP <- (betaP^2) - log(betaP^2) - 1 ## convert z-score to betaP
        suppressWarnings({
            anm <- crops(log(nx),
                         10*log(nx)+1,
                         gaussMean$new(x,m=m,s=s),
                         alg=capa,betaP=betaP)
        })
        str <- paste(sort(unique(anm$mRec)))
        if(length(str)==1){
            str <- c(str,str)
            opt <- str[1]
        }else{
            opt <- paste(anm$mRec[ which.min( anm$betaRec + (anm$betaRec < 3*log(nx))*2e300 ) ])
        }
            
        updateSliderTextInput(
            session = getDefaultReactiveDomain(),
            "profile_nc",
            choices = str,
            selected = opt
        )
        DATA$profile_anom <- list(x=x,m=m,s=s,crops=anm)  
    })
    bindEvent(profile_run_anom,DATA$profile,input$profile_betaP,
              input$profile_row,input$profile_train)

    
    output$profile_plot <- renderPlot({
        anm <- DATA$profile_anom
        if( is.null(anm) ){ return(NULL) }
        msg <- ""
        
        if( sum(is.finite(anm$x))==0){msg <- "No observations available"}
        if( sum(is.finite(anm$m+anm$s))==0){msg <- paste(msg,"Unable to compute expected distribution",sep=" \n ")}
        if( (nchar(msg)==0)  & (input$profile_nc %in% names(anm$crops$outRec)) ){
            pd <- data.frame(x = names(anm$x),value <- as.numeric(anm$x))
            pd$m <- anm$m
            
            pa <- point_anomalies( anm$crops$outRec[[input$profile_nc]] )
            ca <- collective_anomalies( anm$crops$outRec[[input$profile_nc]] )
            ca$start <- ca$start - 0.5
            ca$end <- ca$end + 0.5
            clr <- rep("grey",nrow(pd))
            clr[ pa$loc ] <- "red"
            ggplot(pd) +
                geom_rect(aes(xmin=start, xmax=end,ymin=-Inf,ymax=Inf),alpha=0.2,data=ca,fill="red") + #,fill=clr) +
                geom_bar(stat="identity",aes(x=x,y=value),col=clr,fill=clr,na.rm=TRUE) +
                geom_point(aes(x=x,y=m),pch=19,size=5,col="blue") + xlab(NULL)+
                ylab("kWh") +
                theme_bw(base_size = 20)
            
        }else{
            ggplot() + 
                annotate("text", x = 4, y = 25, size=8,
                         label = msg) +
                theme_void()
        }
    })
    
    output$profile_summaryTable <- renderTable({
        anm <- DATA$profile_anom
        if(is.null(anm) | !(input$profile_nc %in% names(anm$crops$outRec))){ NULL }
        else{
            d <- data.frame(x = as.numeric(anm$x), m = as.numeric(anm$m),s = as.numeric(anm$s))
            tmp <- summary( anm$crops$outRec[[input$profile_nc]] )
            tmp[["Change"]] <- numeric(nrow(tmp))
            tmp[["Total Change"]] <- numeric(nrow(tmp))
            for(ii in seq_len(nrow(tmp))){
                if(tmp$type[ii] == "background"){ next }
                idx <- tmp$start[ii]:tmp$end[ii]
                e <- d$x[idx] - d$m[idx]
                tmp[["Change"]][ii] <- mean( e ,na.rm=TRUE )
                tmp[["Total Change"]][ii] <- mean( e ,na.rm=TRUE ) * length(idx)
            }
            tmp$Start <- rownames(DATA$profile)[tmp$start]
            tmp$End <- rownames(DATA$profile)[tmp$end]
            tmp$start <- tmp$end <- NULL
            cap <- function(s){
                paste0(toupper(substring(s, 1, 1)),substring(s, 2))
            }
            tmp$Type <- sapply(tmp$type,cap)
            tmp[,c("Type","Start","End","Change","Total Change")]
        }
    })



}

shinyApp(ui = ui, server = server)

