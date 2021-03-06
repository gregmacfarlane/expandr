#' ipfr: A package to perform iterative proportional fitting
#' 
#' There are currently three functions 
#' 
#' \code{\link{ipf}} (deprecated)
#' 
#' \code{\link{ipu}}
#' 
#' \code{\link{ipu_nr}}
#' 
#' Both \code{ipu} and \code{ipu_nr} implement list balancing. \code{ipu_nr}
#' implements a newton-raphson approach to balance primary and secondary
#' targets. \code{ipu} balances targets directly, which can mean faster
#' convergence.
#' 
#' @docType package
#' 
#' @name ipfr
NULL
#> NULL

#' Iterative Proportional Updating
#' 
#' @description A general case of iterative proportional fitting. It can satisfy
#'   two, disparate sets of marginals that do not agree on a single total. A
#'   common example is balancing population data using household- and person-level
#'   marginal controls. This could be for survey expansion or synthetic
#'   population creation. The second set of marginal/seed data is optional, meaning
#'   it can also be used for more basic IPF tasks.
#'   
#'   Vignette: \url{http://pbsag.github.io/ipfr/}
#' 
#' @references \url{http://www.scag.ca.gov/Documents/PopulationSynthesizerPaper_TRB.pdf}
#' 
#' @param primary_seed In population synthesis or household survey expansion, 
#'   this would be the household seed table (each record would represent a 
#'   household). It could also be a trip table, where each row represents an 
#'   origin-destination pair. Must contain a \code{pid} ("primary ID") field
#'   that is unique for each row. Must also contain a geography field that
#'   starts with "geo_".
#' 
#' @param primary_targets A \code{named list} of data frames.  Each name in the 
#'   list defines a marginal dimension and must match a column from the 
#'   \code{primary_seed} table. The data frame associated with each named list
#'   element must contain a geography field (starts with "geo_"). Each row in
#'   the target table defines a new geography (these could be TAZs, tracts,
#'   clusters, etc.). The other column names define the marginal categories that
#'   targets are provided for. The vignette provides more detail.
#' 
#' @param secondary_seed Most commonly, if the primary_seed describes households, the 
#'   secondary seed table would describe a unique person with each row. Must
#'   also contain the \code{pid} column that links each person to their 
#'   respective household in \code{primary_seed}. Must not contain any geography
#'   fields (starting with "geo_").
#' 
#' @param secondary_targets Same format as \code{primary_targets}, but they constrain 
#'   the \code{secondary_seed} table.
#'   
#' @param secondary_importance A \code{real} between 0 and 1 signifying the 
#'   importance of the secondary targets. At an importance of 1, the function
#'   will try to match the secondary tarets exactly. At 0, only the percentage
#'   distributions are used (see the vignette section "Target Agreement".)
#' 
#' @param relative_gap After each iteration, the weights are compared to the
#' previous weights and the %RMSE is calculated. If the %RMSE is less than
#' the \code{relative_gap} threshold, then the process terminates.
#' 
#' @param max_iterations maximimum number of iterations to perform, even if 
#'    \code{relative_gap} is not reached.
#'    
#' @param absolute_diff Upon completion, the \code{ipu()} function will report
#'   the worst-performing marginal category and geography based on the percent
#'   difference from the target. \code{absolute_diff} is a threshold below which
#'   percent differences don't matter.
#'   
#'   For example, if if a target value was 2, and the expanded weights equaled
#'   1, that's a 100% difference, but is not important because the absolute value
#'   is only 1.
#'   
#'   Defaults to 10.
#'   
#' @param weight_floor Minimum weight to allow in any cell to prevent zero
#'   weights. Set to .0001 by default.  Should be arbitrarily small compared to
#'   your seed table weights.
#'   
#' @param verbose Print iteration details and worst marginal stats upon 
#'   completion? Default \code{FALSE}.
#'   
#' @param max_ratio \code{real} number. The average weight per seed record is
#' calculated by dividing the total of the targets by the number of records.
#' The max_scale caps the maximum weight at a multiple of that average. Defaults
#' to \code{10000} (basically turned off).
#' 
#' @param min_ratio \code{real} number. The average weight per seed record is
#' calculated by dividing the total of the targets by the number of records.
#' The min_scale caps the minimum weight at a multiple of that average. Defaults
#' to \code{0.0001} (basically turned off).
#' 
#' @return a \code{named list} with the \code{primary_seed} with weight, a 
#'   histogram of the weight distribution, and two comparison tables to aid in
#'   reporting.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' hh_seed <- data.frame(
#'   pid = c(1, 2, 3, 4),
#'   siz = c(1, 2, 2, 1),
#'   weight = c(1, 1, 1, 1),
#'   geo_cluster = c(1, 1, 2, 2)
#' )
#' 
#' hh_targets <- list()
#' hh_targets$siz <- data.frame(
#'   geo_cluster = c(1, 2),
#'   `1` = c(75, 100),
#'   `2` = c(25, 150)
#' )
#' 
#' result <- ipu(hh_seed, hh_targets, max_iterations = 10)
#' }
#' 
#' @importFrom magrittr "%>%"

ipu <- function(primary_seed, primary_targets, 
                secondary_seed = NULL, secondary_targets = NULL,
                secondary_importance = 1,
                relative_gap = 0.01, max_iterations = 100, absolute_diff = 10,
                weight_floor = .00001, verbose = FALSE,
                max_ratio = 10000, min_ratio = .0001){

  # If person data is provided, both seed and targets must be
  if (xor(!is.null(secondary_seed), !is.null(secondary_targets))) {
    stop("You provided either secondary_seed or secondary_targets, but not both.")
  }
  
  # Check for valid values of secondary_importance.
  if (secondary_importance > 1 | secondary_importance < 0) {
    stop("`secondary_importance` argument must be between 0 and 1")
  }
  
  # Check hh and person tables
  if (!is.null(secondary_seed)) {
    check_tables(primary_seed, primary_targets, secondary_seed, secondary_targets)
  } else {
    check_tables(primary_seed, primary_targets)
  }
  
  # Scale target tables. 
  # All tables in the list will match the totals of the first table.
  primary_targets <- scale_targets(primary_targets, verbose)
  if (!is.null(secondary_seed)) {
    secondary_targets <- scale_targets(secondary_targets, verbose) 
  }
  
  # Balance secondary targets to primary.
  if (secondary_importance != 1 & !is.null(secondary_seed)){
    if (verbose) {message("Balancing secondary targets to primary")}
    secondary_targets_mod <- balance_secondary_targets(
      primary_targets, primary_seed, secondary_targets, secondary_seed,
      secondary_importance
    )
  } else {
    secondary_targets_mod <- secondary_targets
  }
  
  # Pull off the geo information into a separate equivalency table
  # to be used as needed.
  geo_equiv <- primary_seed %>%
    dplyr::select(dplyr::starts_with("geo_"), "pid")
  primary_seed_mod <- primary_seed %>%
    dplyr::select(-dplyr::starts_with("geo_"))
  
  # Remove any fields that aren't in the target list and change the ones
  # that are to factors.
  col_names <- names(primary_targets)
  primary_seed_mod <- primary_seed_mod %>%
    # Keep only the fields of interest (marginal columns and pid)
    dplyr::select(dplyr::one_of(c(col_names, "pid"))) %>%
    # Convert to factors and then to dummy columns if the column has more
    # than one category.
    dplyr::mutate_at(
      .vars = col_names,
      # .funs = dplyr::funs(ifelse(length(unique(.)) > 1, as.factor(.), .))
      .funs = dplyr::funs(as.factor(.))
    )
  # If one of the columns has only one value, it cannot be a factor. The name
  # must also be changed to match what the rest will be after one-hot encoding.
  for (name in col_names){
    if (length(unique(primary_seed_mod[[name]])) == 1) {
      # unfactor
      primary_seed_mod[[name]] <- type.convert(as.character(primary_seed_mod[[name]]))
      # change name
      value = primary_seed_mod[[name]][1]
      new_name <- paste0(name, ".", value)
      names(primary_seed_mod)[names(primary_seed_mod) == name] <- new_name
    }
  }
  # Use one-hot encoding to convert the remaining factor fields to dummies
  primary_seed_mod <- primary_seed_mod %>%
    mlr::createDummyFeatures()
  
  if (!is.null(secondary_seed)) {
    # Modify the person seed table the same way, but sum by primary ID
    col_names <- names(secondary_targets_mod)
    secondary_seed_mod <- secondary_seed %>%
      # Keep only the fields of interest
      dplyr::select(dplyr::one_of(c(col_names, "pid"))) %>%
      dplyr::mutate_at(
        .vars = col_names,
        .funs = dplyr::funs(as.factor(.))
      ) %>%
      mlr::createDummyFeatures() %>%
      dplyr::group_by(pid) %>%
      dplyr::summarize_all(
        .funs = sum
      )
    
    # combine the hh and per seed tables into a single table
    seed <- primary_seed_mod %>%
      dplyr::left_join(secondary_seed_mod, by = "pid")
  } else {
    seed <- primary_seed_mod
  }
  
  # Add the geo information back.
  seed <- seed %>%
    dplyr::mutate(weight = 1)  %>%
    dplyr::left_join(geo_equiv, by = "pid")
  
  # store a vector of attribute column names to loop over later.
  # don't include 'pid' or 'weight' in the vector.
  geo_pos <- grep("geo_", colnames(seed))
  pid_pos <- grep("pid", colnames(seed))
  weight_pos <- grep("weight", colnames(seed))
  seed_attribute_cols <- colnames(seed)[-c(geo_pos, pid_pos, weight_pos)]
  
  # modify the targets to match the new seed column names and
  # join them to the seed table
  if (!is.null(secondary_seed)) {
    targets <- c(primary_targets, secondary_targets_mod)
  } else {
    targets <- primary_targets
  }
  for (name in names(targets)) {
    # targets[[name]] <- targets[[name]] %>%
    temp <- targets[[name]] %>%
      tidyr::gather(key = "key", value = "target", -dplyr::starts_with("geo_")) %>%
      dplyr::mutate(key = paste0(!!name, ".", key, ".target")) %>%
      tidyr::spread(key = key, value = target)
    
    # Get the name of the geo column
    pos <- grep("geo_", colnames(temp))
    geo_colname <- colnames(temp)[pos]
    
    seed <- seed %>%
      dplyr::left_join(temp, by = geo_colname)
  }
  
  # Calculate average, min, and max weights and join to seed. If there are
  # multiple geographies in the first primary target table, then min and max
  # weights will vary by geography.
  pos <- grep("geo_", colnames(targets[[1]]))
  geo_colname <- colnames(targets[[1]])[pos]
  recs_by_geo <- seed %>%
    dplyr::group_by(!!as.name(geo_colname)) %>%
    dplyr::summarize(count = n())
  weight_scale <- targets[[1]] %>%
    tidyr::gather(key = category, value = total, -!!as.name(geo_colname)) %>%
    dplyr::group_by(!!as.name(geo_colname)) %>%
    dplyr::summarize(total = sum(total)) %>% 
    dplyr::left_join(recs_by_geo, by = geo_colname) %>%
    dplyr::mutate(
      avg_weight = total / count,
      min_weight = (!!min_ratio) * avg_weight,
      max_weight = (!!max_ratio) * avg_weight
    ) 
  seed <- seed %>%
    dplyr::left_join(weight_scale, by = geo_colname)
  
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations) {
    # Loop over each target and upate weights
    for (seed_attribute in seed_attribute_cols) {
      
      # Create lookups for targets list
      target_tbl_name <- strsplit(seed_attribute, ".", fixed = TRUE)[[1]][1]
      target_name <- paste0(seed_attribute, ".", "target")
      
      # Get the name of the geo column
      target_tbl <- targets[[target_tbl_name]]
      pos <- grep("geo_", colnames(target_tbl))
      geo_colname <- colnames(target_tbl)[pos]

      # Adjust weights
      seed <- seed %>%
        dplyr::mutate(
          geo = !!as.name(geo_colname),
          attr = !!as.name(seed_attribute),
          target = !!as.name(target_name)
        ) %>%
        dplyr::group_by(geo) %>%
        dplyr::mutate(
          total_weight = sum(attr * weight),
          factor = ifelse(attr > 0, target / total_weight, 1),
          weight = weight * factor,
          # Implement the floor on zero weights
          weight = pmax(weight, weight_floor),
          # Cap weights to to multiples of the average weight.
          # Not applicable if target is 0.
          weight = ifelse(attr > 0 & target > 0, pmax(min_weight, weight), weight),
          weight = ifelse(attr > 0 & target > 0, pmin(max_weight, weight), weight)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-geo, -attr, -target, -factor)
    }
    
    # Determine percent differences (by geo field)
    saved_diff_tbl <- NULL
    pct_diff <- 0
    for (seed_attribute in seed_attribute_cols) {
      # create lookups for targets list
      target_tbl_name <- strsplit(seed_attribute, ".", fixed = TRUE)[[1]][1]
      target_name <- paste0(seed_attribute, ".", "target")
      target_tbl <- targets[[target_tbl_name]]
      
      # Get the name of the geo column
      pos <- grep("geo_", colnames(target_tbl))
      geo_colname <- colnames(target_tbl)[pos]
      
      diff_tbl <- seed %>%
        dplyr::filter((!!as.name(seed_attribute)) > 0) %>%
        dplyr::select(
          geo = !!geo_colname, pid, attr = !!seed_attribute, weight,
          target = !!target_name
        ) %>%
        dplyr::group_by(geo) %>%
        dplyr::mutate(
          total_weight = sum(attr * weight),
          diff = total_weight - target,
          abs_diff = abs(diff),
          pct_diff = diff / (target + .0000001) # avoid dividing by zero
        ) %>%
        # Removes rows where the absolute gap is smaller than 'absolute_diff'
        dplyr::filter(abs_diff > absolute_diff) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
      
      # If any records are left in the diff_tbl, record worst percent difference 
      # and save that percent difference table for reporting.
      if (nrow(diff_tbl) > 0) {
        if (max(abs(diff_tbl$pct_diff)) > pct_diff) {
          pct_diff <- max(abs(diff_tbl$pct_diff))
          saved_diff_tbl <- diff_tbl
          saved_category <- seed_attribute
          saved_geo <- geo_colname
        }
      }
      
    }
    
    # Test for convergence
    if (iter > 1) {
      rmse <- mlr::measureRMSE(prev_weights, seed$weight)
      pct_rmse <- rmse / mean(prev_weights) * 100
      converged <- ifelse(pct_rmse <= relative_gap, TRUE, FALSE)
      if(verbose){
        cat("\r Finished iteration ", iter, ". %RMSE = ", pct_rmse)
      }
    }
    prev_weights <- seed$weight
    iter <- iter + 1
  }

  if (verbose) {
    message(ifelse(converged, "IPU converged", "IPU did not converge"))
    if (is.null(saved_diff_tbl)) {
      message("All targets matched within the absolute_diff of ", absolute_diff)
    } else {
      message("Worst marginal stats:")
      position <- which(abs(saved_diff_tbl$pct_diff) == pct_diff)[1]
      message("Category: ", saved_category)
      message(saved_geo, ": ", saved_diff_tbl$geo[position])
      message("Worst % Diff: ", round(
        saved_diff_tbl$pct_diff[position] * 100, 2), "%"
      )
      message("Difference: ", round(saved_diff_tbl$diff[position], 2))
    }
    utils::flush.console()
  }
  
  # Set final weights into primary seed table. Also include average weight
  # and distribution info.
  primary_seed$weight <- seed$weight
  primary_seed$avg_weight <- seed$avg_weight
  primary_seed$weight_factor <- primary_seed$weight / primary_seed$avg_weight
  
  # If the average weight is 0 (meaning the target was 0) set weight
  # and weight factor to 0.
  primary_seed <- primary_seed %>%
    mutate(
      weight = ifelse(avg_weight == 0, 0, weight),
      weight_factor = ifelse(avg_weight == 0, 0, weight_factor)
    )
  
  # Create the result list (what will be returned). Add the seed table and a
  # histogram of weight distribution.
  result <- list()
  result$weight_tbl <- primary_seed
  result$weight_dist <- ggplot2::ggplot(
    data = primary_seed, ggplot2::aes(primary_seed$weight_factor)
  ) +
    ggplot2::geom_histogram(bins = 10, fill = "darkblue", color = "gray") +
    ggplot2::labs(
      x = "Weight Ratio = Weight / Average Weight", y = "Count of Seed Records"
    )
  
  # Compare resulting weights to initial targets
  primary_comp <- compare_results(primary_seed, primary_targets)
  result$primary_comp <- primary_comp
  if (!is.null(secondary_seed)) {
    # Add geo fields to secondary seed
    pos <- grep("geo_", colnames(primary_seed))
    geo_cols <- colnames(primary_seed)[pos]
    seed <- secondary_seed %>%
      dplyr::left_join(
        primary_seed %>% dplyr::select(dplyr::one_of(geo_cols), pid, weight),
        by = "pid"
      )
    
    # Run the comparison against the original, unscaled targets 
    # and store in 'result'
    secondary_comp <- compare_results(
      seed, 
      secondary_targets
    )
    result$secondary_comp <- secondary_comp
  }
  
  return(result)
}

#' Check seed and target tables for completeness
#' 
#' @description Given seed and targets, checks to make sure that at least one
#'   observation of each marginal category exists in the seed table.  Otherwise,
#'   ipf/ipu would produce wrong answers without throwing errors.
#'
#' @inheritParams ipu

check_tables <- function(primary_seed, primary_targets, secondary_seed = NULL, secondary_targets = NULL){
  
  # If person data is provided, both seed and targets must be
  if (xor(!is.null(secondary_seed), !is.null(secondary_targets))) {
    stop("You provided either secondary_seed or secondary_targets, but not both.")
  }
  
  ## Primary checks ##
  
  # Check that there are no NA values in seed or targets
  if (any(is.na(unlist(primary_seed)))) {
    stop("primary_seed table contains NAs")
  }
  if (any(is.na(unlist(primary_targets)))) {
    stop("primary_targets table contains NAs")
  }
  
  # Check that primary_seed table has a pid field and that it has a unique
  # value on each row.
  if (!"pid" %in% colnames(primary_seed)) {
    stop("The primary seed table does not have field 'pid'.")
  }
  unique_pids <- unique(primary_seed$pid)
  if (length(unique_pids) != nrow(primary_seed)) {
    stop("The primary seed's pid field has duplicate values.")
  }
  
  # check hh tables for correctness
  for (name in names(primary_targets)) {
    tbl <- primary_targets[[name]]
    
    # Check that each target table has a geo field
    check <- grepl("geo_", colnames(tbl))
    if (!any(check)) {
      stop("primary_target table '", name, "' does not have a geo column (must start with 'geo_')")
    }
    if (sum(check) > 1) {
      stop("primary_target table '", name, "' has more than one geo column (starts with 'geo_'")
    }
    
    # Get the name of the geo field
    pos <- grep("geo_", colnames(tbl))
    geo_colname <- colnames(tbl)[pos]
    
    # Get vector of other column names
    col_names <- colnames(tbl)
    col_names <- type.convert(col_names[!col_names == geo_colname], as.is = TRUE)
    
    # Check that at least one observation of the current target is in every geo
    for (geo in unique(unlist(primary_seed[, geo_colname]))){
      test <- match(col_names, primary_seed[[name]][primary_seed[, geo_colname] == geo])
      if (any(is.na(test))) {
        prob_cat <- col_names[which(is.na(test))]
        stop(
          "Marginal ", name, ", category ", prob_cat[1], " is missing from ",
          geo_colname, " ", geo, " in the primary_seed table."
        )
      }   
    }
  }
  
  
  ## Secondary checks (if provided) ##
  
  if (!is.null(secondary_seed)) {
    # Check for NAs
    if (any(is.na(unlist(secondary_seed)))) {
      stop("secondary_seed table contains NAs")
    }
    if (any(is.na(unlist(secondary_targets)))) {
      stop("secondary_targets table contains NAs")
    }
    
    # Check that secondary seed table has a pid field
    if (!"pid" %in% colnames(secondary_seed)) {
      stop("The primary seed table does not have field 'pid'.")
    }
    
    # Check that the secondary seed table does not have any geo columns
    check <- grepl("geo_", colnames(secondary_seed))
    if (any(check)) {
      stop("Do not include geo fields in the secondary_seed table (primary_seed only).")
    }
    
    # check the per tables for correctness
    for (name in names(secondary_targets)) {
      tbl <- secondary_targets[[name]]
      
      # Check that each target table has a geo field
      check <- grepl("geo_", colnames(tbl))
      if (!any(check)) {
        stop("secondary_target table '", name, "' does not have a geo column (must start with 'geo_')")
      }
      if (sum(check) > 1) {
        stop("secondary_target table '", name, "' has more than one geo column (starts with 'geo_'")
      }
      
      # Get the name of the geo field
      pos <- grep("geo_", colnames(tbl))
      geo_colname <- colnames(tbl)[pos]
      
      # Add the geo field from the primary_seed before checking
      secondary_seed <- secondary_seed %>%
        dplyr::left_join(
          primary_seed %>% dplyr::select(pid, geo_colname),
          by = "pid"
        )
      
      # Get vector of other column names
      col_names <- colnames(tbl)
      col_names <- type.convert(col_names[!col_names == geo_colname], as.is = TRUE)
      
      # Check that at least one observation of the current target is in every geo
      for (geo in unique(unlist(secondary_seed[, geo_colname]))){  
        test <- match(col_names, secondary_seed[[name]][secondary_seed[, geo_colname] == geo])
        if (any(is.na(test))) {
          prob_cat <- col_names[which(is.na(test))]
          stop(
            "Marginal ", name, ", category ", prob_cat[1], " is missing from ",
            geo_colname, " ", geo, " in the secondary_seed table."
          )
        }   
      }
    }
  }
}


#' Compare results to targets
#'
#' @param seed \code{data.frame} Seed table with a weight column in the same
#' format required by \code{ipu()}.
#' 
#' @param targets \code{named list} of \code{data.frames} in the same format
#' required by \code{ipu()}.
#'
#' @return \code{data frame} comparing balanced results to targets

compare_results <- function(seed, targets){
  
  # Expand the target tables out into a single, long-form data frame
  comparison_tbl <- NULL
  for (name in names(targets)){
    
    # Pull out the current target table
    target <- targets[[name]]
    
    # Get the name of the geo field
    pos <- grep("geo_", colnames(target))
    geo_colname <- colnames(target)[pos]
    
    # Gather the current target table into long form
    target <- target %>%
      dplyr::mutate(geo = paste0(geo_colname, "_", !!as.name(geo_colname))) %>%
      dplyr::select(-dplyr::one_of(geo_colname)) %>%
      tidyr::gather(key = category, value = target, -geo) %>%
      dplyr::mutate(category = paste0(name, "_", category))
    
    # summarize the seed table
    result <- seed %>%
      dplyr::select(geo = !!as.name(geo_colname), category = !!as.name(name), weight) %>%
      dplyr::mutate(
        geo = paste0(geo_colname, "_", geo),
        category = paste0(name, "_", category)
      ) %>%
      dplyr::group_by(geo, category) %>%
      dplyr::summarize(result = sum(weight))
    
    # Join them together
    joined_tbl <- target %>%
      dplyr::left_join(result, by = c("geo" = "geo", "category" = "category"))
    
    # Append it to the master target df
    comparison_tbl <- dplyr::bind_rows(comparison_tbl, joined_tbl)
  }
  
  # Calculate difference and percent difference
  comparison_tbl <- comparison_tbl %>%
    dplyr::mutate(
      diff = result - target,
      pct_diff = round(diff / target * 100, 2),
      diff = round(diff, 2)
    ) %>%
    dplyr::arrange(geo, category)
  
  return(comparison_tbl)
}


#' Scale targets to ensure consistency
#' 
#' Often, different marginals may disagree on the total number of units. In the
#' context of household survey expansion, for example, one marginal might say
#' there are 100k households while another says there are 101k. This function
#' solves the problem by scaling all target tables to match the first target
#' table provided.
#' 
#' @param targets \code{named list} of \code{data.frames} in the same format
#' required by \link{ipu}. 
#' 
#' @param verbose \code{logical} Show a warning for each target scaled?
#'   Defaults to \code{FALSE}.
#' 
#' @return A \code{named list} with the scaled targets
#' 
scale_targets <- function(targets, verbose = FALSE){
  
  for (i in c(1:length(names(targets)))) {
    name <- names(targets)[i]
    target <- targets[[name]]
    
    # Get the name of the geo field
    pos <- grep("geo_", colnames(target))
    geo_colname <- colnames(target)[pos]
    
    # calculate total of table
    target <- target %>%
      tidyr::gather(key = category, value = count, -!!geo_colname)
    total <- sum(target$count)
    
    # Start a string that will be used for the warning message if targets
    # are scaled and verbose = TRUE
    warning_msg <- "Scaling target tables: "
    
    # if first iteration, set total to the global total. Otherwise, scale table
    if (i == 1) {
      global_total <- total
      show_warning <- FALSE
    } else {
      fac <- global_total / total
      # Write out warning
      if (fac != 1 & verbose) {
        show_warning <- TRUE
        warning_msg <- paste0(warning_msg, " ", name)
      }
      target <- target %>%
        dplyr::mutate(count = count * !!fac) %>%
        tidyr::spread(key = category, value = count)
      targets[[name]] <- target
    }
  }
  
  if (show_warning) {
    message(warning_msg)
    utils::flush.console()
  }
  
  return(targets)
}

#' Balances secondary targets to primary
#' 
#' The average weight per record needed to satisfy targets is computed for both
#' primary and secondary targets. Often, these can be very different, which leads
#' to poor performance. The algorithm must use extremely large or small weights
#' to match the competing goals. The secondary targets are scaled so that they
#' are consistent with the primary targets on this measurement.
#' 
#' If multiple geographies are present in the secondary_target table, then
#' balancing is done for each geography separately.
#' 
#' @inheritParams ipu
#' 
#' @return \code{named list} of the secondary targets

balance_secondary_targets <- function(primary_targets, primary_seed,
                                      secondary_targets, secondary_seed,
                                      secondary_importance){

  # Extract the first table from the primary target list and geo name
  pri_target <- primary_targets[[1]]
  pos <- grep("geo_", colnames(pri_target))
  pri_geo_colname <- colnames(pri_target)[pos]

  for (name in names(secondary_targets)){
    sec_target <- secondary_targets[[name]]

    # Get geography field
    pos <- grep("geo_", colnames(sec_target))
    sec_geo_colname <- colnames(sec_target)[pos]

    # If the geographies used aren't the same, convert the primary table
    if (pri_geo_colname != sec_geo_colname) {
      pri_target <- pri_target %>%
        dplyr::left_join(
          primary_seed %>% 
            dplyr::select(!!pri_geo_colname, sec_geo_colname) %>%
            dplyr::group_by(!!as.name(pri_geo_colname)) %>%
            dplyr::slice(1),
          by = pri_geo_colname
        ) %>%
        dplyr::select(-dplyr::one_of(pri_geo_colname))
    }

    # Summarize the primary and secondary targets by geography
    pri_target <- pri_target  %>%
      tidyr::gather(key = cat, value = count, -sec_geo_colname) %>%
      dplyr::group_by(!!as.name(sec_geo_colname)) %>%
      dplyr::summarize(total = sum(count))
    sec_target <- sec_target %>%
      tidyr::gather(key = cat, value = count, -sec_geo_colname) %>%
      dplyr::group_by(!!as.name(sec_geo_colname)) %>%
      dplyr::summarize(total = sum(count))

    # Get primary and secondary record counts
    pri_rec_count <- primary_seed %>%
      dplyr::group_by(!!as.name(sec_geo_colname)) %>%
      dplyr::summarize(recs = n())
    sec_rec_count <-secondary_seed %>%
      dplyr::left_join(
        primary_seed %>% dplyr::select(pid, dplyr::one_of(sec_geo_colname)),
        by = "pid"
      ) %>%
      dplyr::group_by(!!as.name(sec_geo_colname)) %>%
      dplyr::summarize(recs = n())

    # Calculate average weights and the secondary factor
    pri_rec_count$avg_weight <- pri_target$total / pri_rec_count$recs
    sec_rec_count$avg_weight <- sec_target$total / sec_rec_count$recs
    sec_rec_count$factor <- adjust_factor(
      pri_rec_count$avg_weight / sec_rec_count$avg_weight,
      # in this context, high importance means you want the final factor
      # in this table to be near 1. Must flip the importance variable.
      1 - secondary_importance
    )

    # Update the secondary targets by the factor
    secondary_targets[[name]] <- secondary_targets[[name]] %>%
      dplyr::left_join(
        sec_rec_count %>% dplyr::select(!!sec_geo_colname, factor),
        by = sec_geo_colname
      ) %>%
      dplyr::mutate_at(
        .vars = dplyr::vars(-factor, -dplyr::one_of(sec_geo_colname)),
        .funs = dplyr::funs(. * factor)
      ) %>%
      dplyr::select(-factor)
  }

  return(secondary_targets)
}

#' Applies an importance weight to an ipfr factor
#' 
#' @description At lower values of importance, the factor is moved closer to 1.
#' 
#' @param factor A correction factor that is calculated using target/current.
#' 
#' @param importance A \code{real} between 0 and 1 signifying the importance of
#'   the factor. A importance of 1 does not modify the factor. An importance of
#'   0.5 would shrink the factor closer to 1.0 by 50 percent.
#'
#' @return The adjusted factor.
#' 

adjust_factor <- function(factor, importance){
  
  # return the same factor if importance = 1
  if (importance == 1) {return(factor)}
  
  if (importance > 1 | importance < 0) {
    stop("`importance` argument must be between 0 and 1")
  }
  
  # Otherwise, return the adjusted factor
  adjusted <- 1 - ((1 - factor) * (importance + .0001))
  return(adjusted)
}





























