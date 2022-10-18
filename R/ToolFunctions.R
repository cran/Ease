# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                TOOL FUNCTIONS                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


## Matrix tests ----


#' Test if a matrix is of probability
#'
#' @param x a matrix.
#'
#' @return A logical corresponding to whether \code{x} is a probability
#' matrix (sum of rows equal to 1).
#'
#' @author Ehouarn Le Faou
#'
is.probability.matrix <- function(x) {
  rowSum <- apply(x, 1, sum)
  if (any(rowSum != 1)) {
    return(FALSE)
  }
  return(TRUE)
}

#' Test if a matrix is a default matrix
#'
#' @param x a matrix.
#'
#' @return A logical corresponding to whether \code{x} is a default matrix
#' (matrix of dimension 0x0).
#'
#' @author Ehouarn Le Faou
#'
is.default.matrix <- function(x) {
  return(nrow(x) == 0 & ncol(x) == 0)
}

#' Test if a matrix is a correct transition matrix
#'
#' @param x a matrix.
#' @param type type of the matrice (mutation matrix ? recombination matrix ?)
#' @param name the name of the matrix.
#'
#' @return A logical corresponding to whether \code{x} is a correct transition
#' matrix, i.e. a square matrix with dimensions greater than 0 and whose rows
#' sum to 1.
#'
#' @author Ehouarn Le Faou
#'
is.correct.transition.matrix <- function(x, type, name) {
  if (!is.default.matrix(x)) { # Not default matrix ?
    if (nrow(x) == ncol(x)) { # Square matrix ?
      if (!is.probability.matrix(x)) { # Is it a probability matrix ?
        stop(paste0(
          "The (or one of the) ", type, " matrix(ces) in '", name,
          "' given as input is not a probability matrix."
        ))
      }
    } else {
      stop(paste0(
        "The (or one of the) ", type, " matrix(ces) in '", name,
        "' given as input is not square."
      ))
    }
  }
}

## Diverse ----

#' Concatenate, print and line break
#'
#' Object output in the same way as the function \link[base]{cat} but adding
#' a line break at the end.
#'
#' See \link[base]{cat}.
#'
#' @param ... see \link[base]{cat}.
#' @param file see \link[base]{cat}.
#' @param sep see \link[base]{cat}.
#' @param fill see \link[base]{cat}.
#' @param labels see \link[base]{cat}.
#' @param append see \link[base]{cat}.
#'
#' @return None (invisible \code{NULL}).
#'
#' @author Ehouarn Le Faou
#'
catn <- function(..., file = "", sep = " ", fill = FALSE, labels = NULL,
                 append = FALSE) {
  cat(..., "\n",
    file = file, sep = sep, fill = fill, labels = labels,
    append = append
  )
}

#' Listing for display
#'
#' Listing from the elements of a vector by producing a string, with comma
#' separation between each element and the word "and" between the last two
#' elements.
#'
#' @param vect a vector of any class.
#'
#' @return A listing of the elements of the input vector as a string.
#'
#' @author Ehouarn Le Faou
#'
listing <- function(vect) {
  if (length(vect) == 1) {
    return(as.character(vect))
  } else {
    le <- length(vect)
    return(
      paste0(
        Reduce(function(x, y) {
          paste(x, ", ", y, sep = "")
        }, vect[-le]),
        " and ",
        vect[le]
      )
    )
  }
}
