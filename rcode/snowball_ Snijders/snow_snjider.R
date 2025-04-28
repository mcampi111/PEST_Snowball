# ================================
# SNOWBALL ESTIMATION FUNCTION
# ================================

#' Estimate hidden population size using snowball sampling
#'
#' @param sample_file Path to the sample file (text format: vertex followed by its contacts)
#' @param indegree_file Optional path to the indegree file (vertex followed by its indegree)
#' @return A list containing the different estimates
snowball_estimate <- function(sample_file, indegree_file = NULL) {
  # --- Read the sample file (.dat) ---
  lines <- readLines(sample_file)
  lines <- lines[lines != ""]  # remove empty lines
  lines <- trimws(lines)       # remove extra spaces
  
  # Parse sample vertices and their outdegrees
  vertices <- c()
  outdegrees <- c()
  
  for (line in lines) {
    parts <- strsplit(line, "\\s+")[[1]]  # split by spaces
    vertex <- as.integer(parts[1])
    connections <- as.integer(parts[-1])
    
    vertices <- c(vertices, vertex)
    outdegrees <- c(outdegrees, length(connections))
  }
  
  # --- Read the indegree file if provided (.din) ---
  indegrees <- NULL
  if (!is.null(indegree_file)) {
    indegree_lines <- readLines(indegree_file)
    indegree_lines <- indegree_lines[indegree_lines != ""]
    indegree_lines <- trimws(indegree_lines)
    
    indegree_data <- do.call(rbind, strsplit(indegree_lines, "\\s+"))
    indegree_data <- as.data.frame(indegree_data, stringsAsFactors = FALSE)
    colnames(indegree_data) <- c("vertex", "indegree")
    indegree_data$vertex <- as.integer(indegree_data$vertex)
    indegree_data$indegree <- as.integer(indegree_data$indegree)
    
    indegrees <- setNames(indegree_data$indegree, indegree_data$vertex)
  }
  
  # --- Compute estimators ---
  n <- length(vertices)
  total_outdegree <- sum(outdegrees)
  
  # 1. Simple estimator
  est_simple <- n + total_outdegree
  
  # 2. Horvitz-Thompson estimator
  est_ht <- NA
  if (!is.null(indegrees)) {
    positive_indegrees <- indegrees[as.character(vertices)]
    positive_indegrees <- positive_indegrees[positive_indegrees > 0]
    
    if (length(positive_indegrees) > 0) {
      est_ht <- sum(1 / positive_indegrees)
    }
  }
  
  # Return results
  list(
    n_sampled = n,
    total_outdegree = total_outdegree,
    estimate_simple = est_simple,
    estimate_horvitz_thompson = est_ht
  )
}



# ================================
# EXAMPLE: Create a Tiny Example and Run
# ================================

# --- Create a small fake .DAT network ---
sample_lines <- c(
  "1 2 3",  # Vertex 1 points to 2 and 3
  "2 3",    # Vertex 2 points to 3
  "3"       # Vertex 3 points to nobody
)

# Write to a .dat file
sample_file <- "sample_example.dat"
writeLines(sample_lines, sample_file)

# --- Create a small fake .DIN indegree file ---
indegree_lines <- c(
  "1 0",
  "2 1",
  "3 2"
)

# Write to a .din file
indegree_file <- "sample_example.din"
writeLines(indegree_lines, indegree_file)

# ================================
# RUN THE FUNCTION
# ================================

# Run without indegree
cat("Running WITHOUT indegree file:\n")
result_no_indegree <- snowball_estimate(sample_file)
print(result_no_indegree)

cat("\n---------------------------------\n")

# Run with indegree
cat("Running WITH indegree file:\n")
result_with_indegree <- snowball_estimate(sample_file, indegree_file)
print(result_with_indegree)



# ================================
# REAL DATA
# ================================



result_SN1 <- snowball_estimate("/Users/mcampi/Downloads/snowball/SN1.DAT")
print(result_SN1)

result_SN2 <- snowball_estimate("/Users/mcampi/Downloads/snowball/SN2.DAT",
                                "/Users/mcampi/Downloads/snowball/SN2.DIN")
print(result_SN2)







# ================================
# MORE ROBUST VERSION 
# ================================




#' Estimate hidden population size using snowball sampling
#'
#' @param sample_file Path to the sample file (.dat)
#' @param indegree_file Optional path to the indegree file (.din)
#' @return A list containing estimates and statistics, or NULL if errors occur
snowball_estimate_robust <- function(sample_file, indegree_file = NULL) {
  # --- Helper function to safely convert to integer ---
  safe_as_integer <- function(x) {
    result <- tryCatch(
      as.integer(x),
      warning = function(w) {
        warning(paste("Conversion to integer failed:", w$message))
        return(NA) # Or another appropriate value
      },
      error = function(e) {
        stop(paste("Conversion to integer error:", e$message))
      }
    )
    return(result)
  }
  
  # --- Check if sample file exists ---
  if (!file.exists(sample_file)) {
    stop(paste("Sample file not found:", sample_file))
    return(NULL)
  }
  
  # --- Read and parse the sample file (.dat) ---
  lines <- readLines(sample_file, warn = FALSE) # Suppress "incomplete final line" warning for now
  lines <- trimws(lines[lines != ""])
  
  vertices <- c()
  outdegrees <- c()
  edges <- list() # Store edges for potential igraph use
  
  for (i in seq_along(lines)) {
    line <- lines[i]
    parts <- strsplit(line, "\\s+")[[1]]
    
    # --- Validate line format ---
    if (length(parts) < 1 || is.na(safe_as_integer(parts[1]))) {
      warning(paste("Invalid format in sample file, line", i, ":", line))
      return(NULL)
    }
    
    vertex <- safe_as_integer(parts[1])
    connections <- safe_as_integer(parts[-1])
    
    if (any(is.na(c(vertex,connections)))){
      warning(paste("Invalid data in sample file, line", i, ":", line))
      return(NULL)
    }
    
    vertices <- c(vertices, vertex)
    outdegrees <- c(outdegrees, length(na.omit(connections))) # Count only valid connections
    edges[[vertex]] <- connections[!is.na(connections)] # Store valid connections
  }
  
  n <- length(vertices)
  total_outdegree <- sum(outdegrees)
  
  # --- Read and parse the indegree file (.din) if provided ---
  indegrees <- NULL
  if (!is.null(indegree_file)) {
    if (!file.exists(indegree_file)) {
      warning(paste("Indegree file not found:", indegree_file))
    } else {
      indegree_lines <- readLines(indegree_file, warn = FALSE)
      indegree_lines <- trimws(indegree_lines[indegree_lines != ""])
      
      indegree_data <- do.call(rbind, strsplit(indegree_lines, "\\s+"))
      indegree_data <- as.data.frame(indegree_data, stringsAsFactors = FALSE)
      colnames(indegree_data) <- c("vertex", "indegree")
      
      indegree_data$vertex <- safe_as_integer(indegree_data$vertex)
      indegree_data$indegree <- safe_as_integer(indegree_data$indegree)
      
      if (any(is.na(indegree_data))){
        warning("Invalid data in indegree file")
        return(NULL)
      }
      
      indegrees <- setNames(indegree_data$indegree, indegree_data$vertex)
    }
  }
  
  # --- Compute estimators ---
  est_simple <- n + total_outdegree
  
  est_ht <- NA
  if (!is.null(indegrees)) {
    positive_indegrees <- indegrees[as.character(vertices)]
    positive_indegrees <- positive_indegrees[positive_indegrees > 0]
    if (length(positive_indegrees) > 0) {
      est_ht <- sum(1 / positive_indegrees)
    }
  }
  
  # --- Return results ---
  return(list(
    n_sampled = n,
    total_outdegree = total_outdegree,
    estimate_simple = est_simple,
    estimate_horvitz_thompson = est_ht,
    edges = edges # Return edges for further analysis
  ))
}

# --- Example Usage ---
# (Keep your example code as is, but use the _robust function)
result_no_indegree <- snowball_estimate_robust("sample_example.dat")
print("result_no_indegree")
print(result_no_indegree)

result_with_indegree <- snowball_estimate_robust("sample_example.dat", "sample_example.din")
print("result_with_indegree")
print(result_with_indegree)

result_SN1_robust <- snowball_estimate_robust("/Users/mcampi/Downloads/snowball/SN1.DAT")
print("result_SN1_robust")
print(result_SN1_robust)

result_SN2_robust <- snowball_estimate_robust("/Users/mcampi/Downloads/snowball/SN2.DAT",
                                              "/Users/mcampi/Downloads/snowball/SN2.DIN")
print("result_SN2_robust")
print(result_SN2_robust)








