#
# PART 0. Load libraries.
#
library(openxlsx)
library(data.table)
library(dplyr)
library(stringr)
library(minpack.lm)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)



#
# PART 1. Prepare initial data.
#
process_excel <- function(excel_sheet, start_row) {
  rep_measurements <- list()
  start_row <- start_row
  end_row <- start_row + 7
  
  cells <- c("Cell1_", "Cell2_")
  concs <- c("Conc100", "Conc10")
  reps <- paste0("_", 1:3)
  
  for (i in 0:2) {
    measurement <- excel_sheet[ (start_row+i*12):(end_row+i*12), -1] %>%
      as.data.frame()
    measurement[] <- lapply(measurement, as.numeric)
    colnames(measurement) <- apply(expand.grid(reps, concs, cells), 1, function(x) paste(rev(x), collapse = ""))
    rownames(measurement) <- c("20mikroM", "2mikroM", "200nanoM", "20nanoM", "2nanoM", "200pikoM", "20pikoM", "Vehicle")
    
    rep_measurements[[as.character(i+1)]] <- measurement
  }
  
  return(rep_measurements)
}

# Load the luminescence data from excel sheet.
excel_sheet <- read.xlsx("~/Downloads/nluc_bodipy_reads.xlsx", startRow = 55)
 
# Process each of the repeated measurements.
nluc <- process_excel(excel_sheet, 6)
bodipy <- process_excel(excel_sheet, 43)

# Take the average of each measurements to produce final averaged data table.
nluc_avg <- Reduce(`+`, nluc) / length(nluc)
bodipy_avg <- Reduce("+", bodipy) / length(bodipy)

# Check which side of the plate has higher BODIPY-cyclopamine fluorescence intensity.
# In order to determine Nluc-SMO from ΔCRD Nluc-SMO. ΔCRD Nluc-SMO has a lower
# dissociation constant for BODIPY-cyclopamine.
mean(as.numeric(unlist(bodipy_avg[ , 1:3]))) # Nluc-SMO
mean(as.numeric(unlist(bodipy_avg[ , 6:9]))) # ΔCRD Nluc-SM



#
# PART 2. Perform competitive binding analysis.
#
# Calculate the BRET ratio.
bret_ratio <- bodipy_avg / nluc_avg

# Replace outlier value with an average.
bret_ratio[1, 4] <- mean(unlist(bret_ratio[1, 5:6]))

# Add concentration information in log scale
concs <- c(20 * 1e-6, 2 * 1e-6, 200 * 1e-9, 20 * 1e-9, 2 * 1e-9, 200 * 1e-12, 20 * 1e-12)
log_concs <- c(log10(concs), -15)

# Average the BRET ratios and calculate standard deviation.
bret_data <- data.frame()

for (i in seq(1, ncol(bret_ratio), by = 3)) {
  data <- bret_ratio[ , i:(i+2) ]
  
  rowMeans(data, na.rm = FALSE)
  apply(data, 1, sd, na.rm = FALSE)
  
  bret_data <- rbind(
    bret_data,
    data.frame(
      rowMeans(data, na.rm = FALSE),
      apply(data, 1, sd, na.rm = FALSE),
      rownames(data),
      log_concs,
      str_split(colnames(data)[1], "_")[[1]][1],
      gsub("Conc", "", str_split(colnames(data)[1], "_")[[1]][2])
    )
  )
}

rownames(bret_data) <- NULL
colnames(bret_data) <- c("bret_ratio", "sd", "ligand_conc", "ligand_log_conc_M", "cell_type", "bodipy_conc")

# Add cell type names and dissociation constant based on this publication - 10.1124/mol.119.118158.
bret_data <- bret_data %>%
  mutate(kd = ifelse(cell_type == "Cell1", 105, 23)) %>%
  mutate(cell_name = ifelse(cell_type == "Cell1", "Nluc-SMO", "∆CRD Nluc-SMO"))

# Create a formula.
formula = bret_ratio ~ bottom + (top - bottom) / (1 + 10^(-logIC50 + ligand_log_conc_M))

# Generate non-linear regression fits for every cell and bodipy-cyclopamine case (4).
nls_fit_results <- list()

for (cl_tp in unique(bret_data$cell_type)) {
  for (bdpy_c in unique(bret_data$bodipy_conc)) {
    data <- bret_data %>%
      filter(cell_type == cl_tp & bodipy_conc == bdpy_c)
    
    # Fit non-linear regression.
    nls_fit <- nlsLM(formula,
                     data = data,
                     start = list(logIC50 = median(data$ligand_log_conc_M),
                                  top = max(data$bret_ratio) + 0.05,
                                  bottom = min(data$bret_ratio) - 0.05))
    
    # Estimate coefficients.
    logIC50 <- unname(coef(nls_fit)[1])
    top <- unname(coef(nls_fit)[2])
    bottom <- unname(coef(nls_fit)[3])
    
    # Calculate Ki.
    Ki <- 10^(logIC50 - log10(1 + as.numeric(bdpy_c) / data$kd[1]))
    
    # Save results.
    nls_fit_results[[paste0(cl_tp, "_", bdpy_c)]] <- list(
      nls_fit = nls_fit,
      logIC50 = logIC50,
      Ki = Ki,
      top = top,
      bottom = bottom,
      data = data
    )
  }
}



#
# PART 3. Generate dose response curve plots.
#
dose_response_curves <- lapply(names(nls_fit_results), function(key) {
  data <- nls_fit_results[[key]]$data
  
  # Create a sequence of ligand concentrations for prediction.
  ligand_log_seq <- seq(
    min(bret_data$ligand_log_conc_M),
    max(bret_data$ligand_log_conc_M),
    length.out = 100
  )
  
  pred_data <- data.frame(ligand_log_conc_M = ligand_log_seq)
  pred_data$bret_ratio <- predict(nls_fit_results[[key]]$nls_fit, 
                                  newdata = pred_data)
  pred_data$cell_name <- data$cell_name[1]
  pred_data$bodipy_conc <- paste0(data$bodipy_conc[1], " nM")
  
  p <- ggplot() +
    geom_point(data = data,
               aes(x = ligand_log_conc_M, 
                   y = bret_ratio),
               alpha = 0.8,
               color = "blue") +
    geom_errorbar(data = data,
                  aes(x = ligand_log_conc_M,
                      ymin = bret_ratio - sd,
                      ymax = bret_ratio + sd),
                  width = 0.1,
                  alpha = 0.5,
                  color = "blue") +
    geom_line(data = pred_data,
              aes(x = ligand_log_conc_M, 
                  y = bret_ratio),
              size = 0.8,
              color = "lightblue") +
    labs(x = expression(log[10]*"[ligand] M"),
         y = "raw BRET ratio",
         title = paste0(data$cell_name[1], " cells with ", data$bodipy_conc, " nM BODIPY")) +
    annotate("text",
             x = -Inf, y = Inf,
             label = paste0("Ki = ", signif(nls_fit_results[[key]]$Ki, 3), " M"),
             hjust = -0.1, vjust = 1.2,
             size = 4, fontface = "italic") +
    theme_minimal()
  
  return(p)
})

# Combine the plots together.
plot <- do.call(ggarrange, c(dose_response_curves, 
                             ncol = 2,
                             nrow = 2))
ggsave("~/Downloads/dose_response_curves.pdf", plot, device = "pdf")



#
# PART 4. Generate heatmap plots for the 96-well plate.
#
bret_matrix <- as.matrix(bret_ratio)
dimnames(bret_matrix) <- list(LETTERS[1:8], as.character(1:12))

cell_name <- rep(c("Nluc-SMO", expression(Delta * "CRD Nluc-SMO")), each = 6)
bodipy_conc <- rep(c("100 nM", "10 nM"), each = 3, times = 2)

# Make top annotation.
top_anno <- HeatmapAnnotation(
  CellType = cell_name,
  BodipyDose = dose,
  col = list(
    CellType = c("Nluc-SMO" = "#66c2a5", "∆CRD Nluc-SMO" = "#fc8d62"),
    BodipyDose = c("100 nM" = "#8da0cb", "10 nM" = "#e78ac3")
  ),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(round(c(min(bret_matrix), 
                        mean(range(bret_matrix)), 
                        max(bret_matrix)), 2), 
                      c("blue", "white", "red"))

col_fun <- circlize::colorRamp2(
  c(min(bret_matrix), max(bret_matrix)),
  c("white", "red")
)

# Generate heatmap.
heatmap <- Heatmap(
  bret_matrix,
  name = "BRET",
  col = col_fun,
  top_annotation = top_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_rot = 0,
  column_names_rot = 0,
  column_title = "96-Well BRET Assay",
  heatmap_legend_param = list(
    title = "BRET Ratio",
    at = round(c(min(bret_matrix), max(bret_matrix)), 2)
  )
)

# Save the heatmap.
pdf("96_well_heatmap.pdf", width = 6, height = 3)
draw(heatmap)
dev.off()












