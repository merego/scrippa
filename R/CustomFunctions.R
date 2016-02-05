# Find the element in a vector x which is the closest to some value 
find.closest <- function(x,value) {
  
  index <- which( min( abs( x - value ) ) == ( abs( x - value ) ) )
  val <- x[index]
  return(list(index,val))
}

# Find indeces corresponding to max on a matrix
which.max.matrix <- function(mat) {
 indexes <- which(mat == max(mat,na.rm=TRUE), arr.ind = TRUE)
 return(indexes)
}

# Find indeces corresponding to min on a matrix
which.min.matrix <- function(mat) {
 indexes <- which(mat == min(mat,na.rm=TRUE), arr.ind = TRUE)
 return(indexes)
}



# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, n=1, centered=FALSE) {

    if (centered) {
        before <- floor  ((n-1)/2)
        after  <- ceiling((n-1)/2)
    } else {
        before <- n-1
        after  <- 0
    }

    # Track the sum and count of number of non-NA items
    s     <- rep(0, length(x))
    count <- rep(0, length(x))

    # Add the centered data 
    new <- x
    # Add to count list wherever there isn't a 
    count <- count + !is.na(new)
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new

    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new   <- c(rep(NA, i), x[1:(length(x)-i)])

        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new

        i <- i+1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new   <- c(x[(i+1):length(x)], rep(NA, i))

        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new

        i <- i+1
    }

    # return sum divided by count
    s/count
}


# Left alignement plots using ggplot and gridExtra
grid.arrange.aligned <- function(plots) {
  grobs <- list()
  widths <- list()

  for (i in 1:length(plots)){
      grobs[[i]] <- ggplotGrob(plots[[i]])
      widths[[i]] <- grobs[[i]]$widths[2:5]
  }

  maxwidth <- do.call(grid::unit.pmax, widths)

  for (i in 1:length(grobs)){
       grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }

  do.call("grid.arrange", c(grobs, ncol = 1))
}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}


lappend <- function (lst, ...){
lst <- c(lst, list(...))
  return(lst)
}
