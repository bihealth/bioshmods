if(interactive()) {
    cntr1 <- data.frame(log2FoldChange=rnorm(5000),
                        pvalue=runif(5000),
                        PrimaryID=paste0("ID", 1:5000))
    cntr2 <- data.frame(log2FoldChange=cntr1$log2FoldChange + 
                                       rnorm(5000),
                        pvalue=runif(5000) * cntr1$pvalue,
                        PrimaryID=paste0("ID", 1:5000))
    cntr2$SecondaryID <- sample(cntr2$PrimaryID)
    rownames(cntr2) <- paste0("ID", 1:5000)
    cntr <- list("Contrast 1"=cntr1, "Contrast 2"=cntr2)
    shinyApp(ui=fluidPage(discoUI("disco", names(cntr))),
             server=function(input, output, session) {
                discoServer("disco", cntr)
             })
}

