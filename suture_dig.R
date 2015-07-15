#READ FILES, PERFORM DIGITIZATION AND CALCULATION
calc.si <- function(file_name, dig_method="dijkstra", smooth=0.4, window_edge=50, window_inc=50, filled_line=F) {
  
  #read file from disk and convert to data frame with x, y coords and pixel intensity
  suture_image <- readbitmap::read.bitmap(paste(file_name, ".png", sep=""))
  suture <- data.frame(x=sort(rep(seq(ncol(suture_image)),nrow(suture_image))), y=rep(seq(nrow(suture_image)), ncol(suture_image)), pixel=c(suture_image))
  
  if (filled_line==T) {
    blackpts <- subset(suture, suture$pixel==0) #extract black pixels
    whitepts <- subset(suture, suture$pixel==1) #extract white pixels
    
    #calculate distance from each black pixel to the nearest white pixel
    blackpts$dists <- apply(blackpts, 1, function(x) {
      sqrt(min((x[1]-whitepts[,1])^2 + (x[2]-whitepts[,2])^2))
    })
    
    #selects only black pixel that touch a white pixel
    suture <- subset(blackpts, blackpts$dists==1)
  }
  
  if (dig_method=="spline") {

    suture_sub <- subset(suture[1:2], suture$pixel == 0) #selects only darkest pixels
    
    plot(suture_sub$x, suture_sub$y, pch=15, cex=0.5, ylim=rev(range(suture_sub$y))) #plots raw suture pixels  
    
    suture.ma <- smooth.spline(suture_sub$x, suture_sub$y, spar = smooth) #cubic spline smoothing to generate line
    
    suture_pts <- data.frame(x=suture.ma$x, y=suture.ma$y) #extracts x/y points from spline fit
    
    points(suture_pts, col=rainbow(nrow(suture_pts)))
    
  } else if (dig_method=="dijkstra") {
    
    path_x <- numeric(0) #creates empty variable to store suture trace
    path_y <- numeric(0)
    
    names(suture) <- c("x", "y", "pixel")  
    
    suture_sub <- subset(suture[1:2], suture$pixel == 0) #selects only darkest pixels
    
    plot(suture_sub$x, suture_sub$y, pch=15, cex=0.5, ylim=rev(range(suture_sub$y))) #plots raw suture pixels  
    
    startx <- min(suture_sub$x) #
    starty <- min(suture_sub$y[which(suture_sub$x==startx)])
    
    while (window_edge < max(suture_sub$x)+(window_inc-1)) {
      
      if (window_edge > max(suture_sub$x)) {
        window_edge <- max(suture_sub$x)+1
      }
      
      suture_part <- subset(suture_sub, suture_sub$x < window_edge)
      
      start_pt <- which(suture_part$x==startx & suture_part$y==starty)
      
      pt_pairs <- data.frame(t(combn(seq(nrow(suture_part)), 2))) #creates list of all possible point pairs (vertex points of graph)
      
      names(pt_pairs) <- c("v1", "v2")
      
      pt_pairs$weight <- as.numeric(dist(suture_part)) #calculates Euclidean distance btwn pairs
      pt_pairs$weight[pt_pairs$weight>1.5] <- 100000
      
      ptgraph <- igraph::graph.data.frame(pt_pairs, directed=F) #creates graph network
      
      for (i in which(suture_part$x == (window_edge-1))) {
        #uses Dijkstra's algorithm with weights to find shortest path between vertices
        sh_path <- igraph::get.shortest.paths(ptgraph, from=start_pt, to=i, weights=igraph::E(ptgraph)$weight)
        
        pt_order <- unlist(sh_path$vpath)
        
        if (length(pt_order) > 2) {
          
          path_x <- c(path_x, suture_part$x[pt_order]) #selects original points corresponding to shortest path
          path_y <- c(path_y, suture_part$y[pt_order])
          
          startx <- path_x[length(path_x)]
          starty <- path_y[length(path_y)]
          
          #delete extra pixels adjacent to line
          visited_points <- data.frame(x=path_x[-length(path_x)], y=path_y[-length(path_y)])
          
          suture_part$dists <- apply(suture_part, 1, function(x) {
            sqrt(min((x[1]-visited_points[,1])^2 + (x[2]-visited_points[,2])^2))
          })
          
          suture_part$final_dists <- apply(suture_part, 1, function(x) {
            sqrt(min((x[1]-path_x[length(path_x)])^2 + (x[2]-path_y[length(path_y)])^2))
          })
          
          adjacent_pts <- paste(subset(suture_part$x, suture_part$dists < 1.5 & suture_part$final_dists > 1.5), subset(suture_part$y, suture_part$dists < 1.5 & suture_part$final_dists > 1.5))
          
          suture_sub <- subset(suture_sub, !paste(suture_sub$x, suture_sub$y) %in% adjacent_pts) #removes visited points (except last visited)
          
          suture_sub <- subset(suture_sub, suture_sub$x > max(path_x) - 0.4*max(suture_sub$x))
          
          #delete 
          suture_sub <- subset(suture_sub, !paste(suture_sub$x, suture_sub$y) %in% paste(path_x, path_y)) #removes visited points (except last visited)
          suture_sub <- rbind(suture_sub, c(startx, starty))
          
          break
        } else if (i == max(which(suture_part$x == (window_edge-1)))) {
          stop("There is a gap in the suture line")
        }
      }
      
      points(path_x, path_y, col=rainbow(length(path_x)), cex=0.5)
      
      window_edge <- window_edge + window_inc
      
    }
    
    
    suture_pts <- data.frame(x=path_x, y=path_y)
    
  } else {stop("Invalid method name")}
  
  #calculate the suture width (in pixels) and enter in SW column
  SW <- max(suture_pts$x) - min(suture_pts$x) 
  
  #calculate the suture index and enter in SI column
  x.final <- suture_pts[nrow(suture_pts), 1]
  x.start <- suture_pts[1, 1]
  y.final <- suture_pts[nrow(suture_pts), 2]
  y.start <- suture_pts[1, 2]
  
  dist.straight <- sqrt((x.final - x.start)^2 + (y.final - y.start)^2) #calculates straight-line distance from start to end
  
  dist.curved <- sum(sqrt(apply(data.frame(diff(suture_pts[, 1])^2, diff(suture_pts[, 2])^2), 1, sum))) #calculates curved distance along suture perimeter
  
  SI <- dist.curved/dist.straight #calculates suture index (SI)
  
  #save the suture points as a file
  write.csv(suture_pts, paste(file_name, ".csv", sep=""))  

  suture_results <- list("SI"=SI, "SW"=SW)
  
  suture_results
}






