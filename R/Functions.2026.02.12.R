

packages.needed <- c("bayestestR", "vegan", "zip")

packages2install <- !packages.needed %in% rownames(installed.packages())
packages.needed <- packages.needed[packages2install]

packages.needed <-c(packages.needed)

install.packages(packages.needed)

ladder <- c(50,  75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475, 500,
            550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)

rox.ladder <- ladder

tightest_triplet <- function(x) {
  x_sorted <- sort(x)
  
  spans <- x_sorted[3:length(x_sorted)] - 
    x_sorted[1:(length(x_sorted) - 2)]
  
  i <- which.min(spans)
  
  x_sorted[i:(i + 2)]
}

scalar1 <- function(x) {x / sqrt(sum(x^2))}

Closest <- function (x, a, which = FALSE, na.rm = FALSE){
  if (na.rm)
    x <- x[!is.na(x)]
  mdist <- min(abs(x - a))
  if (is.na(mdist))
    res <- NA
  else {
    idx <- DescTools::IsZero(abs(x - a) - mdist)
    if (which == TRUE)
      res <- which(idx)
    else res <- x[idx]
  }
  return(res)
}

local.max <- function(a, b){
  return(a[b - 21 + find_peaks(a[(b - 20) : (b + 20)], m = 20, n = 20)])
}

find_peaks <- function (x, m=50, n=50){ # default m is set to %; look at the points +/- % in the vector;
  # I addin "n" so you can chnage the range you look both behind and in front of the point
  shape <- diff(sign(diff(x, na.pad = FALSE))) # use sign to show the dirrection of change
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1 #
    z <- ifelse(z > 0, z, 1)
    w <- i + n + 1 #
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

find_area <- function(x, peaks,  m=50, n=50){

  if(length(peaks) > 0){

    #incase oeaks are at begining or end of vwctor
    peaks <- peaks + 51
    x <- c(rep(0, 50), x, rep(0, 50))

    area <- NULL


    for(j in 1:length(peaks)){

      # plot(x[(peaks[j] - m) : (peaks[j] + n)])

      area[j] <- bayestestR::auc(1 : 101 , x[(peaks[j] - m) : (peaks[j] + n)], "trapezoid")

    }
    return(area)
  }
}

#' extrapolate
#'
#' This function extrapolate the size of each peak in the query given a ladder of know size
#' @param revise.cutoff default = TRUE
#' @keywords lm
#' @export
#' @examples
#' interpolate()

extrapolate <- function(stand.bp, stand.time, sample){

  stand.lm <- lm(stand.bp ~ stand.time)

  res <- predict(stand.lm,
                 newdata = data.frame(stand.time = sample))


  return(res)

}

#' interpolate
#'
#' This function interpolates the size of each peak in the query given a ladder of know size
#' @keywords lm
#' @export
#' @examples
#' interpolate()
interpolate <- function(stand.bp, stand.time, sample){

  if(length(sample) > 0){

    res <- NULL

    for(i in 1:length(sample)){

      a <- sample[i] - stand.time

      if(length(unique(a > 0)) == 2){

        pick <- which(a == min(a[a >= 0]))

        bp <- stand.bp[pick : (pick + 1)]
        time <- stand.time[pick : (pick + 1)]

        stand.lm <- lm(bp ~ time)

        res.i <- predict(stand.lm,
                         newdata = data.frame(time = sample[i]))

        #res[i] <- round(res.i)
        res[i] <- res.i

      }else{

        #print("can't interpolate")
        res[i] <- NA
      }

    }

    return(res)
  }
}

#' filter.bp
#'
#' This function removes DNA fragments outside a given size range
#' @keywords fsa
#' @export
#' @examples
#' filter.bp()
filter.bp <- function(res.bp, bp.min = 200,  bp.max = 600){

  if(class(res.bp) == "data.frame"){

    res <- res.bp[res.bp$bp >= bp.min & res.bp$bp <= bp.max & !is.na(res.bp$bp) , ]
    return(res)
  }else{
    warning("Filter.bp need a dataframe input")
    return("no peaks found")
  }
}

#' Plot fsa Function
#'
#' This function allows you to visualize the chromatogram.Only plots ladder and query channel
#' @param revise.cutoff Would you like to change the cutoff for finding peaks? Defaults to TRUE.
#' @keywords fsa
#' @export
#' @examples
#' plot.fsa()
plot.fsa <- function(file_path,
                     ladder = rox.ladder,
                     cutoff = 250,
                     revise.cutoff = T,
                     sd = 3,
                     only.up = T,
                     channel.ladder = 4,
                     channel.query = 1){
  # the cutoff applies to all channels
  # Revise cutoff only applies to the ladder channel
  # In the ladder channel fragment sizes are assigned to peaks from the largest to smallest fragments of DNA.
  # if revise.cutoff = TRUE, then all peaks in the ladder channel are called,
  # then the median and sd are calculated from the tallest peaks.
  # If the ladder has 16 fragments than the 16 largest peaks will be used for this calculation.
  # anything under the number of specified standard deviations from the median is removed when the
  # revised cutoff is implemented.


  ladder.id <- deparse(substitute(ladder))

  cdiff <- try(read.fsa(file_path))
  if(class(cdiff) != "fsa"){
    return("fsa not able to be parsed")
  }

  plot(cdiff[, 1], type = "n",
       ylim = range(cdiff) * 1.05,
       ylab = "", xlab = "")

  abline(h = cutoff, lty = 3)

  num.channels <- length(  attr(cdiff, "colors"))

  channels <- c(channel.query, channel.ladder)
  dyes <- attr(cdiff, "colors")
  dyes <- dyes[channels]

  res <- list()# save output

  for(i in channels){

    channel.i <- cdiff[, i]
    channel.i.2 <- cdiff[, i]
    lines(channel.i,
          col = attr(cdiff, "colors")[i],
          type = "l")

    channel.i[channel.i <  cutoff] <-  cutoff
    time <- find_peaks(channel.i, 50, 50)
    height <- channel.i.2[time]

    if(revise.cutoff & i == channel.ladder){

      time <- find_peaks(channel.i.2, 50, 50)
      height <- channel.i.2[time]

      tallest <- rev(sort(height))[1:length(ladder)]
      new.cutoff <-  median(  tallest ) - (sd(  tallest ) * sd)

      if(only.up){
        new.cutoff<- max(c(new.cutoff, cutoff))
      }

      channel.i.2[channel.i.2 <  new.cutoff] <-  new.cutoff
      time <- find_peaks(channel.i.2, 50, 50)

      int <- diff(time)
      time <- time[c(T, int < median(int) * 2)]

      height <- channel.i.2[time]

      abline(h = new.cutoff, lty = 3)
      text(1000, new.cutoff,
           paste("     ladder cutoff", round(new.cutoff)), pos = 3)


    }

    points(time,   height,
           bg = attr(cdiff, "colors")[i],
           pch = 21,
           col = i)

    area <- find_area(x = channel.i.2, peaks = time)
    res[[i]] <- as.data.frame(cbind(time, height, area))

  }

  names(res) <- names(attr(cdiff, "colors"))[channels]

  legend("topright",
         pch = 21,
         col = channels,
         pt.bg = dyes,
         bg = "snow1",
         legend = paste(names(dyes),
                        "  wavelength =",
                        attr(cdiff, "wavelengths")[channels],
                        "NM"))
  legend("bottomright", paste("file path =", file_path))

  l.r <- rev(ladder)
  h.r <- rev(res[[channel.ladder]]$height)
  t.r <- rev(res[[channel.ladder]]$time)

  height <- rev(h.r[1 : min(c(length(l.r), length(h.r)))])
  time <- rev(t.r[1 : min(c(length(l.r), length(t.r)))])
  ladder <- rev(l.r[1 : min(c(length(l.r), length(h.r)))])

  ######
  
  triplet <- which(time %in% tightest_triplet(time[5:(length(time) - 5)]))[1]
  adj  <-  triplet - which(ladder == 450) + 1
  
  if(length(adj) > 0){
    time <- time[adj : length(time)]
    height <- height[adj : length(time)]
    ladder <- ladder[1:(length(ladder) - adj + 1)]
  }else{
    message("Warning: ladder looks off")
  }
  ######
  
  text(time,
       height,
       paste("        ", ladder),
       srt = 45, cex = .8, col = "grey")

  
  text(time[length(time)],
       height[length(height)],
       ladder.id,
       col = attr(cdiff, "colors")[channel.ladder],
       pos = 4)

  bp <- interpolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)

  # extrapolated bp:
  #  ext <- extrapolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)
  #  ext[!is.na(int)] <- NA

  res.bp <- cbind(res[[channel.query]], bp)
  # res.bp <- res.bp[res.bp$bp >= 200 & !is.na(res.bp$bp),]

  if(nrow(res.bp) > 0){
    text(res.bp$time, res.bp$height,
         round(res.bp$bp), font = 4, pos = 3)

    return(res.bp)
  }else{
    return("no peaks found")
  }
}

#' Plot fsa all Function
#'
#' This function allows you to visualize the chromatogram. Plotsall channels
#' @param revise.cutoff Would you like to change the cutoff for finding peaks? Defaults to TRUE.
#' @keywords fsa
#' @export
#' @examples
#' plot.fsa.all()
plot.fsa.all <- function(file_path,
                          ladder,
                          cutoff = 250,
                          revise.cutoff = T,
                          sd = 3,
                          channel.ladder = 4,
                          channel.query = 1){
  # all channels are plotted

  # the cutoff applies to all channels
  # Revise cutoff only applies to the ladder channel
  # In the ladder channel fragment sizes are assigned to peaks from the largest to smallest fragments of DNA.
  # if revise.cutoff = TRUE, then all peaks in the ladder channel are called,
  # then the median and sd are calculated from the tallest peaks.
  # If the ladder has 16 fragments than the 16 largest peaks will be used for this calculation.
  # anything under the number of specified standard deviations from the median is removed when the
  # revised cutoff is implemented.


  ladder.id <- deparse(substitute(ladder))
  cdiff <- read.fsa(file_path)

  plot(cdiff[, 1], type = "n",
       ylim = range(cdiff) * 1.05,
       ylab = "", xlab = "")

  abline(h = cutoff, lty = 3)

  num.channels <- length(  attr(cdiff, "colors"))

  res <- list()# save output

  for(i in 1 : num.channels){

    channel.i <- cdiff[, i]
    channel.i.2 <- cdiff[, i]
    lines(channel.i,
          col = attr(cdiff, "colors")[i],
          type = "l")

    channel.i[channel.i <  cutoff] <-  cutoff
    time <- find_peaks(channel.i, 50, 50)
    height <- channel.i.2[time]

    if(revise.cutoff & i == channel.ladder){

      time <- find_peaks(channel.i.2, 50, 50)
      height <- channel.i.2[time]

      tallest <- rev(sort(height))[1:length(ladder)]
      new.cutoff <-  median(  tallest ) - (sd(  tallest ) * sd)

      channel.i.2[channel.i.2 <  new.cutoff] <-  new.cutoff
      time <- find_peaks(channel.i.2, 150, 150)

      int <- diff(time)
      time <- time[c(T, int < median(int) * 2)]

      height <- channel.i.2[time]

      abline(h = new.cutoff, lty = 3)
      text(1000, new.cutoff,
           paste("revised ladder cutoff", round(new.cutoff)), pos = 3)


    }

    points(time,   height,
           bg = attr(cdiff, "colors")[i],
           pch = 21,
           col = i)

    area <- find_area(x = channel.i.2, peaks = time)
    res[[i]] <- as.data.frame(cbind(time, height, area))

  }

  names(res) <- names(attr(cdiff, "colors"))

  legend("topright",
         pch = 21,
         col = 1:   num.channels,
         pt.bg = attr(cdiff, "colors"),
         legend = paste(names(attr(cdiff, "colors")),
                        "  wavelength =",
                        attr(cdiff, "wavelengths"),
                        "NM"))
  legend("bottomright", paste("file path =", file_path))

  l.r <- rev(ladder)
  h.r <- rev(res[[channel.ladder]]$height)
  t.r <- rev(res[[channel.ladder]]$time)

  height <- rev(h.r[1 : min(c(length(l.r), length(h.r)))])
  time <- rev(t.r[1 : min(c(length(l.r), length(t.r)))])
  ladder <- rev(l.r[1 : min(c(length(l.r), length(h.r)))])

  text(time, height,
       paste("        ", ladder)[1:length(ladder)],
       srt = 45, cex = .8, col = "grey")

  text(time[length(time)],
       height[length(height)],
       ladder.id,
       col = attr(cdiff, "colors")[channel.ladder],
       pos = 4)

  bp <- interpolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)

  # extrapolated bp:
  #  ext <- extrapolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)
  #  ext[!is.na(int)] <- NA

  res.bp <- cbind(res[[channel.query]], bp)
 # res.bp <- res.bp[res.bp$bp >= 200 & !is.na(res.bp$bp),]

  text(res.bp$time, res.bp$height,
       round(res.bp$bp), font = 4, pos = 3)

  return(res.bp)
}

#' Summerize fsa Function
#'
#' This function allows you to visualize the chromatogram.Only plots ladder and query channel
#' the cutoff applies to all channels
#' Revise cutoff only applies to the ladder channel
#' In the ladder channel fragment sizes are assigned to peaks from the largest to smallest fragments of DNA.
#' if revise.cutoff = TRUE, then all peaks in the ladder channel are called,
#' then the median and sd are calculated from the tallest peaks.
#' If the ladder has 16 fragments than the 16 largest peaks will be used for this calculation.
#' anything under the number of specified standard deviations from the median is removed when the
#' revised cutoff is implemented.
#' @param revise.cutoff Would you like to change the cutoff for finding peaks? Defaults to TRUE.
#' @keywords fsa summary
#' @export
#' @examples
#' sum.fsa()
sum.fsa <- function(file_path,
                    ladder = rox.ladder,
                    cutoff = 250,
                    sd = 3,
                    revise.cutoff = T,
                    only.up = T,
                    look4triplet = T, 
                    channel.ladder = 4,
                    channel.query = 1){

    ladder.id <- deparse(substitute(ladder))

    cdiff <- try(read.fsa(file_path))
    if(class(cdiff) != "fsa"){
      return("fsa not able to be parsed")
    }

    num.channels <- length(  attr(cdiff, "colors"))

    channels <- c(channel.query, channel.ladder)

    res <- list()# save output

    for(i in channels){

      channel.i <- cdiff[, i]
      channel.i.2 <- cdiff[, i]

      channel.i[channel.i <  cutoff] <-  cutoff
      time <- find_peaks(channel.i, 50, 50)
      height <- channel.i.2[time]

      if(revise.cutoff & i == channel.ladder){

        time <- find_peaks(channel.i.2, 50, 50)
        height <- channel.i.2[time]

        tallest <- rev(sort(height))[1:length(ladder)]
        new.cutoff <-  median(  tallest ) - (sd(  tallest ) * sd)

        if(only.up){
          new.cutoff<- max(c(new.cutoff, cutoff))
        }

        channel.i.2[channel.i.2 <  new.cutoff] <-  new.cutoff
        time <- find_peaks(channel.i.2, 150, 150)

        int <- diff(time)
        time <- time[c(T, int < median(int) * 2)]

        height <- channel.i.2[time]
      }

      area <- find_area(x = channel.i.2, peaks = time)
      res[[i]] <- as.data.frame(cbind(time, height, area))

    }

    l.r <- rev(ladder)
    h.r <- rev(res[[channel.ladder]]$height)
    t.r <- rev(res[[channel.ladder]]$time)

    height <- rev(h.r[1 : min(c(length(l.r), length(h.r)))])
    time <- rev(t.r[1 : min(c(length(l.r), length(t.r)))])
    ladder <- rev(l.r[1 : min(c(length(l.r), length(h.r)))])

    #############
    
    if(look4triplet){
      if(length(time) >= 12){
        if(all(!ladder %in% c(450, 475, 500))){
          message("Warning: peaks at 450, 475, 500 to found, is this a Rox Ladder?")
        }
      triplet <- which(time %in% tightest_triplet(time[4:(length(time) - 4)]))[1]
      
      adj  <-  triplet - which(ladder == 450) + 1
      
      time <- time[adj : length(time)]
      height <- height[adj : length(time)]
      ladder <- ladder[1:(length(ladder) - adj + 1)]
      }else{
        message("Warning: Ladder looks off")
      }
    }
    #############
    
        
    
    bp <- interpolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)

    # extrapolated bp:
    #  ext <- extrapolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)
    #  ext[!is.na(int)] <- NA

    res.bp <- cbind(res[[channel.query]], bp)
    # res.bp <- res.bp[res.bp$bp >= 200 & !is.na(res.bp$bp),]

    return(res.bp)
}

#' Pairwise frag summary
#'
#' This function allows you calculate the Bray Curtis (BC) or Jaccard distance between two
#' fragment summaries
#' @param  window the largest difference in fragment size (bp) between two
#' fragement summaries at which peaks/fragments are considered an observation
#' of the same fragment default is 2.5bp.
#' @param  method Default distacne is Bray-Curtis,"jd" Jaccard
#' @keywords fsa summary Bray-Curtis BC
#' @export
#' @examples
#' pairwise.dist()


pairwise.dist <- function(str.1 = database[[1]],  str.2, window = 2.5, method = "bc"){

  same.1 <- same.2 <- NULL

  for(i in 1:nrow(str.1)){

    a <- abs(str.1$bp[i] - str.2$bp)

    if(min(a) <= window){

      same.1 <- c(same.1, i)
      same.2 <- c(same.2, which(a == min(a))[1])
      # same.2 <- c(same.2, which(a <= 5))
    }
  }

  str.1$bp[same.1] <- str.2$bp[same.2]

  h.1 <- str.1$height
  names(h.1) <- str.1$bp

  h.2 <- str.2$height
  names(h.2) <- str.2$bp

  ######################## deduplicate names

  h.11 <- aggregate(h.1, list(names(h.1)), sum)
  h.1 <- h.11$x
  names(h.1) <- h.11$Group.1


  h.21 <- aggregate(h.2, list(names(h.2)), sum)
  h.2 <- h.21$x
  names(h.2) <- h.21$Group.1

  #######################

  aa <- names(h.2)[!names(h.2) %in% names(h.1) ]
  bb <- vector(length = length(aa))
  names(bb) <- aa
  h.1 <- c(h.1, bb)
  h.1 <- h.1[order(names(h.1))]

  aa <- names(h.1)[!names(h.1) %in% names(h.2) ]
  bb <- vector(length = length(aa))
  names(bb) <- aa
  h.2 <- c(h.2, bb)
  h.2 <- h.2[order(names(h.2))]

  hs <- rbind(h.1, h.2)
  hs <- hs * 25# needed for the log transformation, rarefy needs integers
  hs <- round(hs)

  set.seed(42)
  hs <- vegan::rrarefy(hs, min(rowSums(hs)))

  jd <- vegan::vegdist(hs, method = "jaccard")
  bc <- vegan::vegdist(hs)

  if(method == "bc"){

    return(bc)
  }
  if(method == "jd"){

    return(jd)
  }
}



compare.frags <- function(query.file =  "Promega_001_96-well-plate_ready-to-go/001-G09.fsa",
                          query.ladder = rox.ladder,
                          query.channel= 1,
                          query.channel.ladder = 4,
                          hit.inDB,
                          hit.ladder = rox.ladder,
                          hit.channel= 1,
                          hit.channel.ladder = 4,
                          db = cd.db,
                          meth = "bc",
                          look4triplet = T,
                          cutoff = 200){

  str.1 <- sum.fsa(file_path = query.file,
                   ladder = query.ladder,
                   channel.query = query.channel,
                   channel.ladder = query.channel.ladder,
                   look4triplet = look4triplet,
                   cutoff = cutoff)
  str.1 <- filter.bp(str.1)

  str.2 <- db$database[[which(names(db$database) == hit.inDB)]]
  dist <- pairwise.dist(str.1, str.2, method = meth)
  
  cdiff1 <- read.fsa(query.file)
  cdiff2 <- read.fsa(hit.inDB)

  cdiff1 <-   cdiff1[, query.channel]
  cdiff2 <-   cdiff2[, hit.channel]

  cdiff1[cdiff1 < 0] <- 0
  cdiff2[cdiff2 < 0] <- 0

  bp.1 <- str.1$bp
  time.1 <- str.1$time

  lm.1 <- lm(bp.1 ~ time.1)
  str1.x <- predict(lm.1, newdata = data.frame(time.1 = 1 : length(cdiff1)))

  cdiff1 <- (cdiff1 / max(cdiff1[  str1.x > 150 &   str1.x < 650])) * 100
  # normaliz to peak value in the 150 to 650 bp range

  plot(str1.x, cdiff1, 
       type = "l",
       xlim = c(150, 650), 
       ylim = c(-20,100),
       xlab = "bp", 
       ylab = "peak intensity",
       lwd = 2)

  text( str1.x[str.1$time],
        cdiff1[str.1$time],
        round(str.1$bp, 0), pos = 2)

  bp.2 <- str.2$bp
  time.2 <- str.2$time

  lm.2 <- lm(bp.2 ~ time.2)
  str2.x <- predict(lm.2, newdata = data.frame(time.2 = 1 : length(cdiff2)))

  cdiff2 <- (cdiff2 / max(cdiff2[  str2.x > 150 &   str2.x < 650])) * 100

  lines(str2.x, cdiff2, type = "l", col = 2, lwd = 2, lty = 1)

  points( str2.x[str.2$time],
          cdiff2[str.2$time],
          col = 2)

  text( str2.x[str.2$time],
        cdiff2[str.2$time],
        round(str.2$bp, 0), pos = 4, col = 2)

  ################ make hit plot ################
  a <- strsplit(query.file, "/")[[1]]
  b <- strsplit(hit.file, "/")[[1]]

  ribo <- b[length(b) - 1]# Check to make sure each folder is a ribotype

  title(paste0(round(dist, 3), "BC distance from\n",
               a[length(a)], "\nto refference ",
               b[length(b)], " (",
               ribo, ")"))

  legend("bottom", lty =1, col = c(1,2),
         legend = c( a[length(a)],
                     b[length(b)]), 
         bg = "snow2")

  res <- c(query.file, hit.file, dist, ribo)
  names(res) <- c("query.file", "hit.file", "BCdist", "ribotype")
  res <- as.matrix(res)
  res <- as.data.frame(t(res))
  res$BCdist <- as.numeric(res$BCdist)

  return(res)
}

find.match <- function(file.query, 
                       ladder = rox.ladder,
                       database = cd.db,
                       meth= "bc", 
                       look4triplet = T, 
                       log.trans = F){

  # the most common arguements used to build the database are used to find the match
  db.settings <- database$settings
  cut.m <- median(sapply(db.settings, "[", "cutoff"))
  sd.m <- median(sapply(db.settings, "[", "sd"))

  print(file.query)

  bp.query <- sum.fsa(file_path = file.query, 
                      ladder = ladder, 
                      cutoff = cut.m, 
                      sd = sd.m, 
                      look4triplet = look4triplet) 
  
  bp.query <- filter.bp(bp.query)

  if(class(  bp.query) != "data.frame"){
    return("no match found")
  }

  if(log.trans){
    bp.query$height <- log(bp.query$height)
    bp.query$area <- log(abs(bp.query$area))
  }

  if( nrow(bp.query) > 0 ){

    db <- database$database
    res <- lapply(db,  function(db) pairwise.dist(str.1 = db, str.2 = bp.query, method = meth))
    res <- unlist(res)

    return(sort(res))

  }else{
    return("no match found")
  }
}


make.database <- function(folder = "F-RibotypingFiles/F001/", ladder,  revise.cutoff= T,
                          cutoff = 100, sd = 3, bp.min = 200, bp.max = 600){

  print("this may take a while please be patient")

  db.fsa <- list.files(folder, full.names = T, recursive = T, pattern = ".fsa")
  str.bp <- sum.fsa(db.fsa[1], ladder, cutoff, sd)
  bp.query <- filter.bp(str.bp, bp.min, bp.max)

  summary.list <- lapply(db.fsa, sum.fsa, ladder, cutoff, sd)

  db <- lapply(db.fsa, function(db.fsa) sum.fsa(file_path = db.fsa,
                                                cutoff = cutoff,
                                                ladder = rox.ladder,
                                                revise.cutoff =  revise.cutoff,
                                                sd = sd))

  db.filtered <- lapply(db, function(db) filter.bp(res.bp = db,
                                                   bp.min = bp.min,
                                                   bp.max = bp.max))

  names(db.filtered) <- db.fsa

  a <- sapply(db.filtered, lengths) # remove empty entries

  no.peaks <- colnames(a)[a[4,] == 0]

  if(length(no.peaks) > 0){
    print(paste("No useable peaks found in",  no.peaks))
  }


  db.filtered <- db.filtered[a[4,] > 0]

  db.args <- c(cutoff, sd, bp.min, bp.max, length(  db.fsa))
  names(db.args) <-  c("cutoff", "sd", "bp.min", "bp.max", "n")
  settings.individual <- rep( list(db.args) , length(db.filtered))
  names(settings.individual) <- names(db.filtered)

  db.final <- list(db.filtered,   settings.individual, ladder)
  names(  db.final) <- c("database", "settings", "ladder")

  t <- gsub("-|:", ".", Sys.time())
  t <- gsub(" ", "_", t)

  rds <- paste0("FSA_database_", t, ".RDS")

  saveRDS(object = db.final, file = rds)
  print(paste("saved database as", rds))

  return(db.final)
}


eval.bin <- function(ribo.pick, db,  bc.cut = .25, plot = T){

  #  Performs Complete-linkage clustering hierarchical of each bin

  ids <- names(db$database)
  ids <- gsub("//", "/", ids)
  ribo <- sapply(strsplit(ids, "/"), "[[", 2)

  db.pick <- db$database[  ribo == ribo.pick]

  n <- length(db.pick)

  if(n == 0){
    stop("ERROR ribotype/bin not found")
  }
  if(n == 1){

    rep <- names(  db.pick)[1]
    names(rep) <- paste0(ribo.pick, ".1")

    res <- list(NA, NA, rep)
    names(res) <- c("db.dist", "groups", "representative")

    print("only on FSA found, bin can not be clustered")


    return(res)

  }else{
    db.dist <- matrix(ncol = n, nrow = n)
    row.names(db.dist) <- ids[  ribo == ribo.pick]
    colnames(db.dist) <- ids[  ribo == ribo.pick]


    for(i in 1 : n){
      for(j in 1 : i){
        if(i != j){

          db.dist[row.names(db.dist) == names(db.pick[i]),
                  colnames(db.dist) == names(db.pick[j])] <- pairwise.dist(db.pick[[i]],
                                                                           db.pick[[j]])

        }
      }
    }

    db.dist <- as.dist(db.dist)
    db.clust <- hclust(db.dist, method = "complete")

    if(plot & n > 2){

      plot(db.clust , hang = -1, main = paste(ribo.pick, "\nn =", n), ylim = c(0,1), col = 4)
      abline(h = bc.cut, lty = 2, col = 2)
    }

    groups  <- cutree(db.clust , h =  bc.cut)
    db.dist <- as.matrix(db.dist )
    cents <- NULL

    for(i in unique(groups)){

      if(length(groups[groups == i]) > 1){

        g.dist <- db.dist[row.names(db.dist) %in% names(groups[groups == i]),
                          row.names(db.dist) %in% names(groups[groups == i])]

        g.mean <- rowMeans(g.dist)
        cents[i] <- names(g.mean[g.mean == min(g.mean)[1]])[1]

      }else{

        cents[i] <- names(groups[groups == i])
      }
    }

    names(cents) <- paste0(ribo.pick, ".", unique(groups))

    res <- list(db.dist, groups, cents)
    names(res) <- c("db.dist", "groups", "representative")

    return(res)
  }
}

fsa2df <- function(file_path){

  cdiff <- read.fsa(file_path)

  num.channels <- length(  attr(cdiff, "names"))

  channels <- cdiff[, 1:num.channels]
  colnames(channels) <- names(attr(cdiff, "colors"))
  return(channels)
}




find.missclass <- function(db){

  entries <- names(db$database)
  ribos <- sapply(strsplit(entries, "/"), "[[", 2)
  res <- as.data.frame(matrix( nrow = length(entries), ncol = 10))

  colnames(res) <- c("ribotype",
                     "hit.1", "hit.2", "hit.3",
                     "ribo.1", "ribo.2", "ribo.3",
                     "dist.1", "dist.2", "dist.3", "reps")

  row.names(res) <- entries
  res$ribotype <- ribos


  for(i in 1 : length(entries)){

    print(paste("evaluating entry", entries[i]))

    ribotype <- strsplit(entries[i], "/")[[1]][2]

    hits3 <- find.match(file.query = entries[i],database =  db)[2:4]

    hit3.r <- sapply(strsplit( names(hits3), "/"), "[[", 2)

    res[ i , 2:4 ] <-    names(hits3)
    res[ i , 5:7 ] <-    hit3.r
    res[ i , 8:10] <-    hits3

    res[i , 11] <- length(ribos[ribotype  == ribos])

  }
  return(res)
}


remove.entry <- function(db, db.entry){

  db.pick <- db
  rm <- which(names(db$database) == db.entry)

  if(length(rm) == 0){

    stop("ERROR entry not found")

  }

  db.pick$database <- db.pick$database[-rm]
  db.pick$settings <- db.pick$settings[-rm]


  return(db.pick)

}

###################
read.abif <- function(filename, max.bytes.in.file = file.info(filename)$size, 
                      pied.de.pilote = 1.2, verbose = FALSE){
  #
  # Suppress warnings when reading strings with internal nul character:
  #
  RTC <- function(x, ...) suppressWarnings(rawToChar(x, ...))
  #
  # Define some shortcuts:
  #
  SInt32 <- function(f, ...) readBin(f, what = "integer", signed = TRUE, endian = "big", size = 4, ...)
  SInt16 <- function(f, ...) readBin(f, what = "integer", signed = TRUE, endian = "big", size = 2, ...)
  SInt8 <- function(f, ...) readBin(f, what = "integer", signed = TRUE, endian = "big", size = 1, ...)
  UInt32 <- function(f, ...) readBin(f, what = "integer", signed = FALSE, endian = "big", size = 4, ...)
  UInt16 <- function(f, ...) readBin(f, what = "integer", signed = FALSE, endian = "big", size = 2, ...)
  UInt8 <- function(f, ...) readBin(f, what = "integer", signed = FALSE, endian = "big", size = 1, ...)
  f32 <- function(f, ...) readBin(f, what = "numeric", size = 4, ...)
  f64 <- function(f, ...) readBin(f, what = "numeric", size = 8, ...)
  #
  # Load raw data in memory:
  #
  fc <- file(filename, open = "rb")
  rawdata <- readBin(fc, what = "raw", n = pied.de.pilote*max.bytes.in.file)
  if(verbose) print(paste("number of bytes in file", filename, "is", length(rawdata)))
  close(fc)
  #
  # Make a list to store results:
  #
  res <- list(Header = NULL, Directory = NA, Data = NA)
  #
  # Header section is 128 bytes long, located at a fixed position at the
  # beginning of the file. We essentially need the number of item and dataoffset
  #
  res$Header$abif <- RTC(rawdata[1:4])
  if(res$Header$abif != "ABIF") stop("file not in ABIF format")
  if(verbose) print("OK: File is in ABIF format")
  
  res$Header$version <- SInt16(rawdata[5:6])
  if(verbose) print(paste("File in ABIF version", res$Header$version/100))
  
  res$Header$DirEntry.name <- rawdata[7:10]
  if(verbose) print(paste("DirEntry name: ", RTC(res$Header$DirEntry.name)))
  
  res$Header$DirEntry.number <- SInt32(rawdata[11:14])
  if(verbose) print(paste("DirEntry number: ", res$Header$DirEntry.number))
  
  res$Header$DirEntry.elementtype <- SInt16(rawdata[15:16])
  if(verbose) print(paste("DirEntry elementtype: ", res$Header$DirEntry.elementtype))
  
  res$Header$DirEntry.elementsize <- SInt16(rawdata[17:18])
  if(verbose) print(paste("DirEntry elementsize: ", res$Header$DirEntry.elementsize))
  
  # This one is important:
  res$Header$numelements <- SInt32(rawdata[19:22])
  if(verbose) print(paste("DirEntry numelements: ", res$Header$numelements))
  
  # This one is important too:
  res$Header$dataoffset <- SInt32(rawdata[27:30])
  if(verbose) print(paste("DirEntry dataoffset: ", res$Header$dataoffset))
  dataoffset <- res$Header$dataoffset + 1 # start position is 1 in R vectors
  
  res$Header$datahandle <- SInt32(rawdata[31:34])
  if(verbose) print(paste("DirEntry datahandle: ", res$Header$datahandle))
  
  res$Header$unused <- SInt16(rawdata[35:128], n = 47)
  # Should be ingnored and set to zero
  res$Header$unused[1:length(res$Header$unused)] <- 0
  if(verbose) print(paste("DirEntry unused: ", length(res$Header$unused), "2-byte integers"))
  
  #
  # The directory is located at the offset specified in the header,
  # and consist of an array of directory entries.
  # We scan the directory to put values in a data.frame:
  #
  dirdf <- data.frame(list(name = character(0)))
  dirdf$name <- as.character(dirdf$name) # force to characters
  
  for(i in seq_len(res$Header$numelements)){
    deb <- (i-1)*res$Header$DirEntry.elementsize + dataoffset
    direntry <- rawdata[deb:(deb + res$Header$DirEntry.elementsize)]
    dirdf[i, "name"] <- RTC(direntry[1:4])
    dirdf[i, "tagnumber"] <- SInt32(direntry[5:8])
    dirdf[i, "elementtype"] <- SInt16(direntry[9:10])
    dirdf[i, "elementsize"] <- SInt16(direntry[11:12])
    dirdf[i, "numelements"] <- SInt32(direntry[13:16])
    dirdf[i, "datasize"] <- SInt32(direntry[17:20])
    dirdf[i, "dataoffset"] <- SInt32(direntry[21:24])
  }
  if(verbose){
    print("Element found:")
    print(dirdf$name)
  }
  #
  # Save Directory and make a list to store data:
  #
  res$Directory <- dirdf
  res$Data <- vector("list", nrow(dirdf))
  names(res$Data) <- paste(dirdf$name, dirdf$tagnumber, sep = ".")
  #
  # Data extraction:
  #
  for(i in seq_len(res$Header$numelements)){
    deb <- (i-1)*res$Header$DirEntry.elementsize + dataoffset
    # Short data are stored in dataoffset directly:
    if(dirdf[i, "datasize"] > 4){
      debinraw <- dirdf[i, "dataoffset"] + 1
    } else {
      debinraw <- deb + 20
    }
    elementtype <- dirdf[i, "elementtype"]
    numelements <- dirdf[i, "numelements"]
    elementsize <- dirdf[i, "elementsize"]
    data <- rawdata[debinraw:(debinraw + numelements*elementsize)]
    # unsigned 8 bits integer:
    if(elementtype == 1) res$Data[[i]] <- UInt8(data, n = numelements)
    # char or signed 8 bits integer
    if(elementtype == 2){
      res$Data[[i]] <- tryCatch(RTC(data),finally=paste(rawToChar(data,multiple=TRUE),collapse=""),error=function(er){cat(paste("an error was detected with the following  message:",er," but this error was fixed\n",sep=" "))})
    }
    
    # unsigned 16 bits integer:
    if(elementtype == 3) res$Data[[i]] <- UInt16(data, n = numelements)
    # short:
    if(elementtype == 4) res$Data[[i]] <- SInt16(data, n = numelements)
    # long:
    if(elementtype == 5) res$Data[[i]] <- SInt32(data, n = numelements)
    # float:
    if(elementtype == 7) res$Data[[i]] <- f32(data, n = numelements)
    # double:
    if(elementtype == 8) res$Data[[i]] <- f64(data, n = numelements)
    # date:
    if(elementtype == 10)
      res$Data[[i]] <- list(year = SInt16(data, n = 1), 
                            month = UInt8(data[-(1:2)], n = 1), day = UInt8(data[-(1:3)], n = 1))
    # time:
    if(elementtype == 11)
      res$Data[[i]] <- list(hour = UInt8(data, n = 1), 
                            minute = UInt8(data[-1], n = 1), second = UInt8(data[-(1:2)], n = 1),
                            hsecond = UInt8(data[-(1:3)], n = 1))
    # bool:
    if(elementtype == 13) res$Data[[i]] <- as.logical(UInt8(data))
    # pString:
    if(elementtype == 18){
      n <- SInt8(rawdata[debinraw])
      pString <- RTC(rawdata[(debinraw+1):(debinraw+n)])
      res$Data[[i]] <- pString
    }
    # cString:
    if(elementtype == 19) res$Data[[i]] <- RTC(data[1:(length(data) - 1) ])
    # user:
    if(elementtype >= 1024) res$Data[[i]] <- data
    # legacy:
    if(elementtype %in% c(12)) 
      warning("unimplemented legacy type found in file")
    if(elementtype %in% c(6, 9, 14, 15, 16, 17, 20, 128, 256, 384))
      warning("unsupported legacy type found in file")
  }
  return(res)
}

### Verbatim copy from MLPA package source 1.10.0
# Imports a .fsa file from Applied Biosystems, using the provided converter
# Author : Sylvain Mareschal <maressyl@gmail.com>
read.fsa <- function(
  file,
  applyLowess = TRUE,
  processed = FALSE,
  meta.extra = NULL,
  ...
) {
  # Parse ABIF
  fsa <- read.abif(file, ...)
  
  # Scan dimensions
  channelCount <- fsa$Data$`Dye#.1`
  if("SCAN.1" %in% names(fsa$Data)) { scanCount <- fsa$Data$SCAN.1
  } else                            { scanCount <- fsa$Data$Scan.1 ### Legacy
  }
  if(scanCount == 0)    warning("No scan stored in this file")
  if(channelCount == 0) warning("No channel stored in this file")
  
  # Extract data
  dyeNames <- character(0)
  dyeWavelengths <- integer(0)
  dyeColors <- character(0)
  
  # DATA indexes
  if(is.na(processed))  { processed <- any(sprintf("DATA.%i", c(9:12, 205:299)) %in% names(fsa$Data))
  }
  if(isTRUE(processed)) { dataTracks <- c(9:12, 205:299)
  } else                { dataTracks <- c(1:4, 105:199)
  }
  
  channelValues <- list()
  for(i in 1:channelCount) {
    # Get raw DATA
    values <- fsa$Data[[ sprintf("DATA.%i", dataTracks[i]) ]]
    
    # Apply lowess to reduce time bias
    if(isTRUE(applyLowess) && scanCount > 0) values <- values - lowess(x=1:length(values), y=values)$y
    
    # Update value matrix
    channelValues[[i]] <- values
    
    # Get dye name
    if(sprintf("DyeN.%i", i) %in% names(fsa$Data))                      dyeNames[i] <- fsa$Data[[ sprintf("DyeN.%i", i) ]]
    if(dyeNames[i] == "" && sprintf("DyeS.%i", i) %in% names(fsa$Data)) dyeNames[i] <- fsa$Data[[ sprintf("DyeS.%i", i) ]]
    
    # Get dye wavelength
    if(sprintf("DyeW.%i", i) %in% names(fsa$Data)) {
      dyeWavelengths[i] <- fsa$Data[[ sprintf("DyeW.%i", i) ]]
      dyeColors[i] <- wav2RGB(dyeWavelengths[i])
    } else {
      dyeWavelengths[i] <- NA
      dyeColors[i] <- NA
    }
  }
  
  # Channel size consistency
  channelLengths <- sapply(channelValues, length)
  if(isTRUE(processed)) scanCount <- channelLengths[1]
  if(any(channelLengths != scanCount)) stop("Data length inconsistency: ", paste(channelLengths, collapse=", "))
  
  # Store channels
  x <- matrix(as.double(NA), ncol=channelCount, nrow=scanCount)
  for(i in 1:channelCount) x[,i] <- channelValues[[i]]
  
  # Automatic colors
  if(all(is.na(dyeColors))) {
    dyeColors <- rainbow(channelCount)
    warning("No dye wavelength found, using arbitrary colors")
  }
  
  # Use dye names
  names(dyeWavelengths) <- dyeNames
  names(dyeColors) <- dyeNames
  colnames(x) <- dyeNames
  
  # Attributes
  attr(x, "lowess") <- isTRUE(applyLowess)
  attr(x, "wavelengths") <- dyeWavelengths
  attr(x, "colors") <- dyeColors
  
  # Metadata to collect
  collect <- c(
    user = "User",
    machine = "MCHN",
    runModule.name = "RMdN",
    runModule.version = "RMdV",
    runProtocole.name = "RPrN",
    runProtocole.version = "RPrV",
    runDate = "RUND",
    runTime = "RUNT",
    meta.extra
  )
  
  # Start collection
  meta <- list()
  for(metaName in names(collect)) {
    values <- fsa$Data[ grep(sprintf("^%s\\.[0-9]+$", collect[ metaName ]), names(fsa$Data)) ]
    if(length(values) > 0L) {
      if(all(sapply(values, is.atomic))) {
        meta[[ metaName ]] <- unlist(values)
        if(length(values) == 1L) names(meta[[ metaName ]]) <- NULL
      } else {
        meta[[ metaName ]] <- as.matrix(as.data.frame(lapply(values, unlist)))
      }
    }
  }
  
  # Reshape dates
  if("runDate" %in% names(meta) && "runTime" %in% names(meta)) {
    dates <- sprintf("%04i-%02i-%02i %02i:%02i:%02i", meta$runDate["year",], meta$runDate["month",], meta$runDate["day",], meta$runTime["hour",], meta$runTime["minute",], meta$runTime["second",])
    names(dates) <- c("runStart", "runStop", "collectionStart", "collectionStop")[ 1:length(dates) ]
    meta$runDate <- NULL
    meta$runTime <- NULL
    meta <- c(meta, dates)
  }
  
  # Injection time (not portable)
  regex <- "^.+<Token>DC_Injection_Time</Token><Value>(.+?)</Value>.+$"
  if("RMdX.1" %in% names(fsa$Data) && grepl(regex, fsa$Data$RMdX.1)) meta$injectionTime <- sub(regex, "\\1", fsa$Data$RMdX.1)
  
  # Store metadata
  attr(x, "runMetaData") <- meta
  
  # Off scale values (if any)
  attr(x, "offScale") <- as.integer(fsa$Data$OfSc.1) + 1L
  
  # S3 class
  class(x) <- "fsa"
  
  return(x)
}

# Converts a light wavelength (in nm) into an RGB color name
# Author	 : RL <http://codingmess.blogspot.fr/2009/05/conversion-of-wavelength-in-nanometers.html>
# Adapted by : Sylvain Mareschal <maressyl@gmail.com>
wav2RGB <- function(wav) {
  # Initialize
  SSS <- R <- G <- B <- integer(length(wav))
  
  # Color
  zone <- wav >= 380 & wav < 440
  R[zone] <- - (wav[zone] - 440) / (440 - 350)
  B[zone] <- 1
  
  zone <- wav >= 440 & wav < 490
  G[zone] <- (wav[zone] - 440) / (490 - 440)
  B[zone] <- 1
  
  zone <- wav >= 490 & wav < 510
  G[zone] <- 1
  B[zone] <- - (wav[zone] - 510) / (510 - 490)
  
  zone <- wav >= 510 & wav < 580
  R[zone] <- (wav[zone] - 510) / (580 - 510)
  G[zone] <- 1
  
  zone <- wav >= 580 & wav < 645
  R[zone] <- 1
  G[zone] <- - (wav[zone] - 645) / (645 - 580)
  
  zone <- wav >= 645 & wav <= 780
  R[zone] <- 1
  
  # Intensity correction
  zone <- wav >= 380 & wav < 420
  SSS[zone] <- 0.3 + 0.7*(wav[zone] - 350) / (420 - 350)
  
  zone <- wav >= 420 & wav <= 700
  SSS[zone] <- 1
  
  zone <- wav > 700 & wav <= 780
  SSS[zone] <- 0.3 + 0.7*(780 - wav[zone]) / (780 - 700)
  
  # Color name
  out <- rgb(SSS*R, SSS*G, SSS*B)
  
  return(out)
}

