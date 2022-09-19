# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                TOOL FUNCTIONS                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#' Concatenate, print and line break
#'
#' Object output in the same way as the function \link[base]{cat} but adding
#' a line break at the end.
#'
#' See \link[base]{cat}.
#'
#' @param ... See \link[base]{cat}.
#' @param file See \link[base]{cat}.
#' @param sep See \link[base]{cat}.
#' @param fill See \link[base]{cat}.
#' @param labels See \link[base]{cat}.
#' @param append See \link[base]{cat}.
#'
#' @return None (invisible NULL).
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
#' @param vect a vector of any class
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
