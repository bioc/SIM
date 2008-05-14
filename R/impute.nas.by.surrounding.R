`impute.nas.by.surrounding` <-
function(dataset, window.size=5) {

  # Create a matrix to hold the results.
  #retval = as.data.frame( matrix(0, nrow(dataset), ncol(dataset)) )
  #rownames(retval) = rownames(dataset)
  #colnames(retval) = colnames(dataset)


  # The following schematic shows the meaning of the variables used in the
  # rest of the function. A (part of a) vector is depicted. "start" indicates
  # the start of the window under scrutiny. "end" indicates the "end" of the
  # window. "offset" indicates the center.
  #
  #
  # |  X  |  X  |  X  |  X  |  X  |  NA  |  X  |  X  |  X  |  X  |
  # |     |     |     |start|     |offset|     | end |     |     |
  #
  # window.radius = (window.size - 1) / 2
  # start  = offset - window.radius
  # end    = offset + window.radius
  #
  method=median
  window.radius = (window.size - 1) / 2


  # Iterate over the columns (samples) in the dataset.
  for (col.idx in 1:ncol(dataset)) {

    # If the current column contains NAs, process it ...
    NAs = is.na(dataset[, col.idx])

    if ( sum(NAs) > 0 ) {
      print( sprintf("Processing column: %i", col.idx) )

      for (offset in 1:nrow(dataset)) {
        # We'll only have to impute things if the entry at "offset" happens to
        # be an NA.
        if ( is.na(dataset[offset, col.idx]) ) {

          start = offset - window.radius
          end   = offset + window.radius

          # If we're at the beginning of the stream, add a couple of items from
          # the other end until we're at "window.size".
          if (start < 0) {
            end   = end + abs(start) + 1
            start = 0
          }

          # If we're at the end of the stream, add a couple of items from the
          # other end until we're at "window.size".
          if (end > nrow(dataset)) {
            start = start - (end - nrow(dataset)) - 1
            end   = nrow(dataset)
          }

          window = dataset[start:end, col.idx]
          NAs.in.window = is.na(window)

          # If there's only a single NA in the window, we can impute it.
          # If not, we'll let it be and remove the row later on (just before
          # we return the dataset).
          if ( sum(NAs.in.window) == 1 ) {
            # Before computing anything, we need to remove the NA from the
            # window. Otherwise, it'll screw up computations.
            indices = !is.na(window)
            valid.window = window[indices]

            imputed.value = method(valid.window)

            dataset[offset, col.idx] = imputed.value

            if (FALSE) {
              print("window:")
              print(window)

              print("indices:")
              print(indices)

              print("valid.window:")
              print(valid.window)

              print("imputed value:")
              print(imputed.value)
              print("")
            }
          }

        }
      }
    }
  }

  return( na.omit(dataset) )
}

