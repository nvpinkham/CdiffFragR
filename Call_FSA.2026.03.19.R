getwd()

source("R/Functions.2026.03.19.R")

unzip("data/F-Ribotyping_Reps.lite.15.zip", exdir = "data")

rox.ladder <- c(50,  75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475, 500,
                550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)

# biologically informative peaks: 
peaks.pick <- c(200,  250,  300,  350,  400,  450,  475,  500,  550,  600)

# Restore the object
cd.db <- readRDS(file = "data/Cdiff_DB_lite.15.RDS")
cd.db.dist <- readRDS(file = "data/Cdiff_DB_lite.15.dist.RDS")

#posible.ribotypes <- CAF::get.element(names(cd.db$database), sep = "/", 2)

p1 <- list.files("Files_to_analyze",
                 recursive = T, 
                 full.names = T,
                 pattern = ".fsa")

ids <-strsplit(p1, split = "[/]")
ids <- sapply(ids, "[[", 2)

jpegs <- gsub("fsa", "jpeg", ids)

current_time <- Sys.time()
hour_of_day <- format(current_time, "%Y.%m.%d-%H")
results.dir <- paste0("Results_", hour_of_day)

dir.create(results.dir, showWarnings = F)


dir.create("Files_analyzed", showWarnings = F)

results <-  as.data.frame(matrix(ncol = 10, nrow = length(p1)))
colnames(results) <- c("query_file",
                       "Confidence",
                       "max_ladder_intensity",
                       "max_query_intensity",
                       "Hit_1", 
                       "ribo_1",  
                       "Dist_1",
                       "hit_1_jpeg",
                       "chormatogram_jpeg", 
                       "warnings")


querys.db <- list()

for(i in 1 : nrow(results)){

  print(i)
  
  chrom.file <- paste0(results.dir, "/chrom_",  jpegs[i])
  
  jpeg(  chrom.file, height = 1200, width = 2400, res = 300)
  
  res.i <- find.match(file.query = p1[i], 
                      ladder = rox.ladder, 
                      peaks.needed = peaks.pick,
                      database = cd.db, 
                      look4triplet = T, 
                      plot.fsa = T)
  
  dev.off()
  
  results$chormatogram_jpeg[i]    <- chrom.file

  if(!identical(res.i,  "no match found")){
    
    results$max_ladder_intensity[i] <- res.i$metrics$max_ladder_intensity
    results$max_query_intensity[i]  <- res.i$metrics$max_query_intensity
    
    ribotypes.i <- get.element(names(res.i$matches), sep = "/", 3)
    
    plot(res.i$matches[1:10], 
         bg = color.groups(ribotypes.i[1:10]),
         pch = 21,
         cex = 2, 
         main = p1[i], 
         ylab = "BC distance",
         xlab = "top 10 matches")
    
    text(res.i$matches[1:10], labels = ribotypes.i[1:10], pos = 4)
    
    res.i.ordered <- res.i$matches[match(colnames(cd.db.dist), names(res.i$matches))]

    all.dists <- rbind(res.i.ordered, cd.db.dist)
    all.dists <- cbind(c(0, res.i.ordered), all.dists)
    
    row.names(all.dists)[1] <- p1[i]
    colnames(all.dists)[1] <- p1[i]
    
    ribotypes.db <- get.element(row.names(all.dists),
                                sep = "/", 
                                3)
    
    
    p <- ribotypes.db %in% unique(ribotypes.i)[1:5]
    p[1] <- T
    
    clust <- hclust(as.dist(all.dists[p,p]),
                    method = "single")
    
    plot(clust, 
         labels = ribotypes.db[p], 
         cex = .65, 
         main = "Custer Dendrogram of 5 most similar ribotypes")
    
################################################################################
    
    hit.file <- paste0(results.dir, "/hit_",  jpegs[i])
    
    jpeg(hit.file, width = 10, height = 6, units = "in", res = 300,
         quality = 100, pointsize = 14)
    
    res1 <- compare.frags(query.file = p1[i] ,
                          hit.inDB = names(res.i$matches[1]))
    
    dev.off()
    
    results$Hit_1[i] <- res1$hit.file
    results$hit_1_jpeg[i] <- hit.file
    
    results$ribo_1[i] <- res1$ribotype
    results$Dist_1[i] <- res1$BCdist
    
    
  }else{
    
    results$Hit_1[i] <- res.i[1]

  }

  file.copy(from = p1[i], to = gsub("Files_to_analyze", "Files_analyzed", p1[i]))

}

results$query_file <- p1
results$Confidence[results$Dist_1 < 0.1] <- "good match"
results$Confidence[results$Dist_1 >=  0.1] <- "questionable match"
results$Confidence[results$Dist_1 >=  0.2] <- "poor match"

write.csv(results, paste0(results.dir, "/SUMMARY.csv"))
print("Success!")

zip(zipfile = paste0("Results", Sys.Date(), ".zip"), files = results.dir)



