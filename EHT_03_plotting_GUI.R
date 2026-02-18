###############################################################################
# µEHT ANALYSIS PIPELINE — STEP 3
# -----------------------------------------------------------------------------
#
# DESCRIPTION:
#   Interactive GUI tool for generating publication-ready plots and basic
#   statistics from merged multi-donor EHT datasets produced using the
#   NovoHeart CTScreen platform.
#
#   IMPORTANT:
#   - This script expects as input ONLY the merged CSV produced by
#     Script 02 (Step 2): "Merged_EHT_AllDonors*.csv" (or equivalent output).
#   - Do not provide raw donor files or Step 1 outputs directly to this step.
#
#   This script:
#     • Loads a merged multi-donor dataset (from Step 2)
#     • Filters data by:
#           - Donor selection
#           - Week_Time selection
#           - The pacing rule:
#               Expected_Hz == 0  OR  Following == "Yes"
#     • Generates per-week plots for each selected metric and pacing frequency:
#           A) Group plots (selected donors and/or merged groups)
#               - Bar = mean
#               - Errorbar = SEM
#               - Points = individual samples colored by Batch
#               - Optional Tukey post-hoc (top N significant comparisons)
#           B) Control vs Patient plots (t-test)
#               - Runs ONLY if Plot_Group contains exactly:
#                     Control and Patient
#               - Requires user-defined donor-to-group merge mapping
#     • Optionally generates Week comparison plots (Week1 vs Week2):
#           - BY_GROUP only (per Plot_Group separately; t-test across weeks)
#           - Runs ONLY if exactly two Week_Time values are in scope
#     • Exports all figures as PNG + PDF into timestamped output folders
#
#   POSITION IN PIPELINE:
#     This represents Script 03 of a modular 3-step EHT processing pipeline.
#
# AUTHOR:
#   Michele Buono (2026)
###############################################################################

suppressPackageStartupMessages({
  library(tcltk)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(tibble)
  library(tidyr)
})

# macOS Tk font safety (best-effort)
try(tcl("font", "configure", "TkDefaultFont", "-family", "Helvetica", "-size", 8), silent = TRUE)
try(tcl("font", "configure", "TkTextFont",    "-family", "Helvetica", "-size", 8), silent = TRUE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

safe_hz <- function(hz) gsub("\\.", "p", as.character(hz))

safe_num <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, ",", ".")
  x <- str_trim(x)
  suppressWarnings(as.numeric(x))
}

remove_outliers_iqr <- function(df, value_col = "value") {
  v <- df[[value_col]]
  if (length(v) < 4) return(df)
  q <- quantile(v, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  if (is.na(iqr) || iqr == 0) return(df)
  low <- q[1] - 1.5 * iqr
  high <- q[2] + 1.5 * iqr
  df %>% filter(.data[[value_col]] >= low, .data[[value_col]] <= high)
}

mean_sem <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(data.frame(y = NA, ymin = NA, ymax = NA))
  m <- mean(x)
  n <- length(x)
  se <- if (n > 1) sd(x) / sqrt(n) else 0
  data.frame(y = m, ymin = m - se, ymax = m + se)
}

p_to_signif <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.0001, "****",
                ifelse(p < 0.001, "***",
                       ifelse(p < 0.01,  "**",
                              ifelse(p < 0.05,  "*", "ns")))))
}

ttest_to_pval_df_twolevel <- function(df, x_col, value_col = "value") {
  x <- df[[x_col]]
  v <- df[[value_col]]
  if (length(unique(x)) != 2) return(NULL)
  n_by <- table(x)
  if (any(n_by < 2)) return(NULL)
  
  tt <- tryCatch(t.test(v ~ x), error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  
  y_max <- max(v, na.rm = TRUE)
  y_rng <- diff(range(v, na.rm = TRUE))
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- max(1, abs(y_max) * 0.2)
  
  levs <- levels(df[[x_col]])
  tibble(
    group1 = levs[1],
    group2 = levs[2],
    p = tt$p.value,
    p.signif = p_to_signif(tt$p.value),
    y.position = y_max + 0.12 * y_rng
  )
}

tukey_to_pvals <- function(tuk, data, top_n = 6) {
  if (is.null(tuk) || nrow(tuk) == 0) return(NULL)
  tb <- as.data.frame(tuk)
  tb$comparison <- rownames(tb)
  comps <- str_split(tb$comparison, "-", simplify = TRUE)
  tb$group1 <- comps[,1]
  tb$group2 <- comps[,2]
  tb$p.adj  <- tb$`p adj`
  tb$p.adj.signif <- p_to_signif(tb$p.adj)
  
  tb <- tb %>% filter(!is.na(p.adj), p.adj < 0.05)
  if (!nrow(tb)) return(NULL)
  
  tb <- tb %>% arrange(p.adj) %>% slice_head(n = top_n)
  
  y_max <- max(data$value, na.rm = TRUE)
  y_rng <- diff(range(data$value, na.rm = TRUE))
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- max(1, abs(y_max) * 0.2)
  tb$y.position <- y_max + seq_len(nrow(tb)) * (0.10 * y_rng)
  
  tb %>% select(group1, group2, p.adj, p.adj.signif, y.position)
}

parse_metric_labels <- function(metric) {
  parts <- str_split(metric, "_", simplify = TRUE)
  if (ncol(parts) < 2) return(list(title = str_to_sentence(metric), ylab = metric))
  title <- str_to_sentence(paste(parts[1, -ncol(parts)], collapse = " "))
  unit  <- parts[1, ncol(parts)]
  unit <- dplyr::case_when(
    unit == "uN" ~ "µN",
    unit == "uM" ~ "µM",
    unit == "ms" ~ "ms",
    unit == "s"  ~ "s",
    TRUE         ~ unit
  )
  list(title = title, ylab = unit)
}

hex_to_rgb <- function(hex) {
  hex <- gsub("#", "", hex)
  c(strtoi(substr(hex, 1, 2), 16L),
    strtoi(substr(hex, 3, 4), 16L),
    strtoi(substr(hex, 5, 6), 16L))
}
rgb_to_hex <- function(rgb) sprintf("#%02X%02X%02X", rgb[1], rgb[2], rgb[3])
avg_hex <- function(hex1, hex2) rgb_to_hex(round((hex_to_rgb(hex1) + hex_to_rgb(hex2)) / 2))

# Configurable Control/Patient merge mapping
make_plot_group <- function(donor,
                            control_donors = character(0),
                            patient_donors = character(0)) {
  donor <- as.character(donor)
  if (length(control_donors) && donor %in% control_donors) return("Control")
  if (length(patient_donors) && donor %in% patient_donors) return("Patient")
  donor
}

build_group_palette <- function(groups_present) {
  donor_colors <- c("C1"="#F4A261", "C2"="#90BE6D", "P1"="#4D9DE0", "P2"="#E76F51")
  
  group_colors <- donor_colors
  group_colors["Control"] <- if (all(c("C1","C2") %in% names(donor_colors))) avg_hex(donor_colors["C1"], donor_colors["C2"]) else "#BBBBBB"
  group_colors["Patient"] <- if (all(c("P1","P2") %in% names(donor_colors))) avg_hex(donor_colors["P1"], donor_colors["P2"]) else "#999999"
  
  group_colors[names(group_colors) %in% groups_present]
}

# ---------------------------------------------------------------------------
# Plot runner
# ---------------------------------------------------------------------------

run_step3_plots <- function(csv_path,
                            root_dir,
                            metric_cols,
                            weeks_keep = NULL,
                            donors_keep = c("C1","C2","P1","P2"),
                            control_donors_merge = character(0),
                            patient_donors_merge = character(0),
                            keep_following_rule = TRUE,
                            jitter_seed = 1,
                            width = 7, height = 5, dpi = 300,
                            plot_mode = c("A_only", "B_only", "A_and_B"),
                            top_n_tukey = 6,
                            do_week_compare = FALSE) {
  
  plot_mode <- match.arg(plot_mode)
  
  # Validate merge selections
  control_donors_merge <- unique(as.character(control_donors_merge))
  patient_donors_merge <- unique(as.character(patient_donors_merge))
  overlap <- intersect(control_donors_merge, patient_donors_merge)
  if (length(overlap)) {
    stop("Merge selection error: these donors are in BOTH Control and Patient merge lists: ",
         paste(overlap, collapse = ", "))
  }
  
  dat <- readr::read_csv(csv_path, show_col_types = FALSE)
  
  req <- c("Donor", "Expected_Hz", "Following", "Week_Time", "Batch")
  miss <- setdiff(req, names(dat))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  metric_cols <- metric_cols[metric_cols %in% names(dat)]
  if (!length(metric_cols)) stop("No valid metric columns selected/found in dataset.")
  
  dat <- dat %>%
    mutate(
      Expected_Hz = suppressWarnings(as.numeric(Expected_Hz)),
      Week_Time   = as.character(Week_Time),
      Following   = as.character(Following),
      Batch       = as.character(Batch),
      Donor       = as.character(Donor)
    ) %>%
    filter(Donor %in% donors_keep)
  
  if (!nrow(dat)) stop("No rows left after filtering for selected donors.")
  
  if (!is.null(weeks_keep) && length(weeks_keep) > 0) {
    dat <- dat %>% filter(Week_Time %in% weeks_keep)
  }
  
  if (keep_following_rule) {
    dat <- dat %>% filter(Expected_Hz == 0 | Following == "Yes")
  }
  
  if (!nrow(dat)) stop("No rows left after filtering (Week/Following rules).")
  if (all(is.na(dat$Batch) | dat$Batch == "")) stop("Batch is all NA/empty after filtering.")
  
  dat <- dat %>% mutate(
    Plot_Group = make_plot_group(
      Donor,
      control_donors = control_donors_merge,
      patient_donors = patient_donors_merge
    )
  )
  
  groups_present <- sort(unique(dat$Plot_Group))
  n_groups <- length(groups_present)
  
  # B only makes sense if ONLY Control and Patient exist
  can_run_B <- identical(sort(groups_present), c("Control","Patient"))
  
  # Week compare only if exactly 2 weeks in scope
  weeks <- sort(unique(dat$Week_Time))
  can_week_compare <- (length(weeks) == 2)
  
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  plots_out_dir <- file.path(root_dir, paste0("Plots_Step03_", timestamp))
  dir.create(plots_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Batch colors
  batches <- sort(unique(dat$Batch))
  batch_colors_pastel <- setNames(hue_pal(l = 80, c = 45)(length(batches)), batches)
  
  group_colors <- build_group_palette(groups_present)
  make_jit <- function() position_jitter(width = 0.12, height = 0, seed = jitter_seed)
  pt_size  <- 2.1
  pt_alpha <- 0.90
  
  hz_vals <- sort(unique(dat$Expected_Hz))
  hz_vals <- hz_vals[!is.na(hz_vals)]
  
  # Counters
  n_saved_A <- 0L
  n_saved_B <- 0L
  n_saved_W_bygroup <- 0L
  n_saved_W_alldonors <- 0L
  
  # A root folder only when needed
  out_dir_A_root <- NULL
  if (plot_mode %in% c("A_only","A_and_B") && n_groups >= 2) {
    out_dir_A_root <- file.path(plots_out_dir, "A_GroupPlots_SelectedOrMerged")
    dir.create(out_dir_A_root, recursive = TRUE, showWarnings = FALSE)
  }
  
  out_dir_B_root <- file.path(plots_out_dir, "B_Control_vs_Patient") # create only if saving
  
  # Week compare root (created only if saving)
  out_dir_W_root <- NULL
  if (isTRUE(do_week_compare) && can_week_compare) {
    out_dir_W_root <- file.path(plots_out_dir, paste0("W_", weeks[1], "_vs_", weeks[2]))
  }
  
  # -------------------------
  # 1) Per-week plots (A and/or B)
  # -------------------------
  for (wk in weeks) {
    
    dat_wk <- dat %>% filter(Week_Time == wk)
    if (!nrow(dat_wk)) next
    
    weekA_root <- if (!is.null(out_dir_A_root)) file.path(out_dir_A_root, paste0("Week_", wk)) else NULL
    
    for (hz in hz_vals) {
      
      df_base <- dat_wk %>% filter(Expected_Hz == hz)
      if (!nrow(df_base)) next
      
      dirA <- if (!is.null(weekA_root)) file.path(weekA_root, paste0("Hz_", safe_hz(hz))) else NULL
      dirB <- file.path(out_dir_B_root, paste0("Week_", wk, "_Hz_", safe_hz(hz)))
      
      for (metric in metric_cols) {
        
        labs_parsed <- parse_metric_labels(metric)
        
        # -------- A) Group plots --------
        if (!is.null(out_dir_A_root) && n_groups >= 2 && plot_mode %in% c("A_only","A_and_B")) {
          
          dfA <- df_base %>%
            transmute(
              Plot_Group = factor(Plot_Group, levels = groups_present),
              Batch = factor(Batch, levels = batches),
              value = safe_num(.data[[metric]])
            ) %>%
            filter(!is.na(Plot_Group), !is.na(value))
          
          if (nrow(dfA) >= 3 && nlevels(droplevels(dfA$Plot_Group)) >= 2) {
            
            dfA_clean <- dfA %>%
              group_by(Plot_Group) %>%
              group_modify(~ remove_outliers_iqr(.x, "value")) %>%
              ungroup()
            
            if (nlevels(droplevels(dfA_clean$Plot_Group)) < 2) next
            
            tuk_pvals <- NULL
            aov_fit <- tryCatch(aov(value ~ Plot_Group, data = dfA_clean), error = function(e) NULL)
            if (!is.null(aov_fit)) {
              tuk <- tryCatch(TukeyHSD(aov_fit)[["Plot_Group"]], error = function(e) NULL)
              tuk_pvals <- tukey_to_pvals(tuk, dfA_clean, top_n = top_n_tukey)
            }
            
            pA <- ggplot(dfA_clean, aes(x = Plot_Group, y = value)) +
              stat_summary(aes(fill = Plot_Group), fun = mean, geom = "bar",
                           width = 0.65, alpha = 0.9, color = "black") +
              stat_summary(fun.data = mean_sem, geom = "errorbar", width = 0.20, linewidth = 0.6) +
              geom_point(aes(color = Batch), position = make_jit(), shape = 16,
                         size = pt_size, alpha = pt_alpha) +
              scale_fill_manual(values = group_colors, name = "Group") +
              scale_color_manual(values = batch_colors_pastel, name = "Batch") +
              guides(color = guide_legend(override.aes = list(size = pt_size + 0.6, alpha = 1))) +
              labs(title = labs_parsed$title, x = NULL, y = labs_parsed$ylab) +
              theme_minimal(base_size = 12) +
              theme(legend.position = "right", plot.title = element_text(face = "bold"))
            
            if (!is.null(tuk_pvals) && nrow(tuk_pvals)) {
              pA <- pA + stat_pvalue_manual(
                tuk_pvals, label = "p.adj.signif",
                tip.length = 0.01, hide.ns = TRUE, inherit.aes = FALSE
              )
            }
            
            safe <- paste0("BAR_A_", gsub("[^A-Za-z0-9]+", "_", metric))
            dir.create(dirA, recursive = TRUE, showWarnings = FALSE)
            ggsave(file.path(dirA, paste0(safe, ".png")), pA, width = width, height = height, dpi = dpi)
            ggsave(file.path(dirA, paste0(safe, ".pdf")), pA, width = width, height = height)
            
            n_saved_A <- n_saved_A + 1L
          }
        }
        
        # -------- B) Control vs Patient --------
        if (plot_mode %in% c("B_only","A_and_B")) {
          
          if (!can_run_B) next
          
          dfB <- df_base %>%
            transmute(
              Group2 = factor(Plot_Group, levels = c("Control","Patient")),
              Batch  = factor(Batch, levels = batches),
              value  = safe_num(.data[[metric]])
            ) %>%
            filter(!is.na(Group2), !is.na(value))
          
          if (nrow(dfB) >= 4 && nlevels(droplevels(dfB$Group2)) == 2) {
            
            dfB_clean <- dfB %>%
              group_by(Group2) %>%
              group_modify(~ remove_outliers_iqr(.x, "value")) %>%
              ungroup()
            
            pvals2 <- ttest_to_pval_df_twolevel(dfB_clean, "Group2", "value")
            
            donor_colors <- c("C1"="#F4A261", "C2"="#90BE6D", "P1"="#4D9DE0", "P2"="#E76F51")
            group2_colors <- c(
              "Control" = avg_hex(donor_colors["C1"], donor_colors["C2"]),
              "Patient" = avg_hex(donor_colors["P1"], donor_colors["P2"])
            )
            
            pB <- ggplot(dfB_clean, aes(x = Group2, y = value)) +
              stat_summary(aes(fill = Group2), fun = mean, geom = "bar",
                           width = 0.65, alpha = 0.9, color = "black") +
              stat_summary(fun.data = mean_sem, geom = "errorbar", width = 0.20, linewidth = 0.6) +
              geom_point(aes(color = Batch), position = make_jit(), shape = 16,
                         size = pt_size, alpha = pt_alpha) +
              scale_fill_manual(values = group2_colors, name = NULL) +
              scale_color_manual(values = batch_colors_pastel, name = "Batch") +
              guides(color = guide_legend(override.aes = list(size = pt_size + 0.6, alpha = 1))) +
              labs(title = labs_parsed$title, x = NULL, y = labs_parsed$ylab) +
              theme_minimal(base_size = 12) +
              theme(legend.position = "right", plot.title = element_text(face = "bold"))
            
            if (!is.null(pvals2) && nrow(pvals2)) {
              pB <- pB + stat_pvalue_manual(
                pvals2, label = "p.signif",
                tip.length = 0.01, hide.ns = TRUE, inherit.aes = FALSE
              )
            }
            
            safe <- paste0("BAR_B_", gsub("[^A-Za-z0-9]+", "_", metric))
            
            if (n_saved_B == 0L) dir.create(out_dir_B_root, recursive = TRUE, showWarnings = FALSE)
            dir.create(dirB, recursive = TRUE, showWarnings = FALSE)
            
            ggsave(file.path(dirB, paste0(safe, ".png")), pB, width = width, height = height, dpi = dpi)
            ggsave(file.path(dirB, paste0(safe, ".pdf")), pB, width = width, height = height)
            
            n_saved_B <- n_saved_B + 1L
          }
        }
      }
    }
  }
  
  # Remove empty B folder (best-effort)
  if (plot_mode %in% c("B_only","A_and_B") && n_saved_B == 0L && dir.exists(out_dir_B_root)) {
    try(unlink(out_dir_B_root, recursive = TRUE, force = TRUE), silent = TRUE)
  }
  
  # -------------------------
  # 2) Week comparison plots (Week1 vs Week2) — BOTH MODES
  #   A) BY_GROUP: per Plot_Group separately (t-test across weeks)
  #   B) ALL_DONORS: one plot per metric+Hz with all donors together (x = Donor)
  # -------------------------
  if (isTRUE(do_week_compare) && can_week_compare) {
    
    wk1 <- weeks[1]; wk2 <- weeks[2]
    
    out_dir_W_bygroup <- file.path(out_dir_W_root, "BY_GROUP")
    out_dir_W_alldonors <- file.path(out_dir_W_root, "ALL_DONORS_GROUPED")
    
    week_cols <- setNames(hue_pal(l = 80, c = 45)(2), c(wk1, wk2))
    dodge_w <- 0.75
    
    # ---- A) BY_GROUP ----
    for (hz in hz_vals) {
      
      df_hz <- dat %>% filter(Expected_Hz == hz)
      if (!nrow(df_hz)) next
      
      for (grp in groups_present) {
        
        df_g <- df_hz %>%
          filter(Plot_Group == grp, Week_Time %in% c(wk1, wk2))
        if (!nrow(df_g)) next
        
        dirW <- file.path(out_dir_W_bygroup, paste0("Group_", grp), paste0("Hz_", safe_hz(hz)))
        
        for (metric in metric_cols) {
          
          labs_parsed <- parse_metric_labels(metric)
          
          dfW <- df_g %>%
            transmute(
              Week_Time = factor(Week_Time, levels = c(wk1, wk2)),
              Batch = factor(Batch, levels = batches),
              value = safe_num(.data[[metric]])
            ) %>%
            filter(!is.na(Week_Time), !is.na(value))
          
          if (nrow(dfW) < 4) next
          n_by <- table(dfW$Week_Time)
          if (any(n_by < 2)) next
          
          dfW_clean <- dfW %>%
            group_by(Week_Time) %>%
            group_modify(~ remove_outliers_iqr(.x, "value")) %>%
            ungroup()
          
          n_by2 <- table(dfW_clean$Week_Time)
          if (length(n_by2) != 2 || any(n_by2 < 2)) next
          
          pvalsW <- ttest_to_pval_df_twolevel(dfW_clean, "Week_Time", "value")
          
          pW <- ggplot(dfW_clean, aes(x = Week_Time, y = value)) +
            stat_summary(aes(fill = Week_Time), fun = mean, geom = "bar",
                         width = 0.65, alpha = 0.9, color = "black") +
            stat_summary(fun.data = mean_sem, geom = "errorbar", width = 0.20, linewidth = 0.6) +
            geom_point(aes(color = Batch), position = make_jit(), shape = 16,
                       size = pt_size, alpha = pt_alpha) +
            scale_fill_manual(values = week_cols, name = "Week") +
            scale_color_manual(values = batch_colors_pastel, name = "Batch") +
            guides(color = guide_legend(override.aes = list(size = pt_size + 0.6, alpha = 1))) +
            labs(
              title = paste0(labs_parsed$title, " — ", grp, " (", wk1, " vs ", wk2, ")"),
              x = NULL, y = labs_parsed$ylab
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "right", plot.title = element_text(face = "bold"))
          
          if (!is.null(pvalsW) && nrow(pvalsW)) {
            pW <- pW + stat_pvalue_manual(
              pvalsW, label = "p.signif",
              tip.length = 0.01, hide.ns = TRUE, inherit.aes = FALSE
            )
          }
          
          safe <- paste0("BAR_W_BYGROUP_", gsub("[^A-Za-z0-9]+", "_", metric))
          
          if ((n_saved_W_bygroup + n_saved_W_alldonors) == 0L) {
            dir.create(out_dir_W_root, recursive = TRUE, showWarnings = FALSE)
          }
          if (n_saved_W_bygroup == 0L) {
            dir.create(out_dir_W_bygroup, recursive = TRUE, showWarnings = FALSE)
          }
          dir.create(dirW, recursive = TRUE, showWarnings = FALSE)
          
          ggsave(file.path(dirW, paste0(safe, ".png")), pW, width = width, height = height, dpi = dpi)
          ggsave(file.path(dirW, paste0(safe, ".pdf")), pW, width = width, height = height)
          
          n_saved_W_bygroup <- n_saved_W_bygroup + 1L
        }
      }
    }
    
    # ---- B) ALL_DONORS grouped ----
    for (hz in hz_vals) {
      
      df_hz <- dat %>%
        filter(Expected_Hz == hz, Week_Time %in% c(wk1, wk2), Donor %in% donors_keep)
      
      if (!nrow(df_hz)) next
      
      dirW_hz <- file.path(out_dir_W_alldonors, paste0("Hz_", safe_hz(hz)))
      
      for (metric in metric_cols) {
        
        labs_parsed <- parse_metric_labels(metric)
        
        dfW <- df_hz %>%
          transmute(
            Donor     = factor(Donor, levels = donors_keep),
            Week_Time = factor(Week_Time, levels = c(wk1, wk2)),
            Batch     = factor(Batch, levels = batches),
            value     = safe_num(.data[[metric]])
          ) %>%
          filter(!is.na(Donor), !is.na(Week_Time), !is.na(value))
        
        if (nrow(dfW) < 6) next
        
        dfW_clean <- dfW %>%
          group_by(Donor, Week_Time) %>%
          group_modify(~ remove_outliers_iqr(.x, "value")) %>%
          ungroup()
        
        donor_ok <- dfW_clean %>%
          count(Donor, Week_Time, name = "n") %>%
          tidyr::pivot_wider(names_from = Week_Time, values_from = n, values_fill = 0) %>%
          mutate(ok = (.data[[wk1]] >= 2 & .data[[wk2]] >= 2)) %>%
          filter(ok) %>%
          pull(Donor) %>%
          as.character()
        
        if (length(donor_ok) < 2) next
        
        dfW_clean <- dfW_clean %>% filter(as.character(Donor) %in% donor_ok)
        
        pW <- ggplot(dfW_clean, aes(x = Donor, y = value, fill = Week_Time)) +
          stat_summary(fun = mean, geom = "bar",
                       position = position_dodge(width = dodge_w),
                       width = 0.65, alpha = 0.9, color = "black") +
          stat_summary(fun.data = mean_sem, geom = "errorbar",
                       position = position_dodge(width = dodge_w),
                       width = 0.20, linewidth = 0.6) +
          geom_point(
            aes(color = Batch),
            position = position_jitterdodge(
              jitter.width = 0.12, dodge.width = dodge_w, seed = jitter_seed
            ),
            shape = 16, size = pt_size, alpha = pt_alpha
          ) +
          scale_fill_manual(values = week_cols, name = "Week") +
          scale_color_manual(values = batch_colors_pastel, name = "Batch") +
          guides(color = guide_legend(override.aes = list(size = pt_size + 0.6, alpha = 1))) +
          labs(
            title = paste0(labs_parsed$title, " — Week comparison (", wk1, " vs ", wk2, ")"),
            x = NULL, y = labs_parsed$ylab
          ) +
          theme_minimal(base_size = 12) +
          theme(legend.position = "right", plot.title = element_text(face = "bold"))
        
        safe <- paste0("BAR_W_ALLDONORS_", gsub("[^A-Za-z0-9]+", "_", metric))
        
        if ((n_saved_W_bygroup + n_saved_W_alldonors) == 0L) {
          dir.create(out_dir_W_root, recursive = TRUE, showWarnings = FALSE)
        }
        if (n_saved_W_alldonors == 0L) {
          dir.create(out_dir_W_alldonors, recursive = TRUE, showWarnings = FALSE)
        }
        dir.create(dirW_hz, recursive = TRUE, showWarnings = FALSE)
        
        ggsave(file.path(dirW_hz, paste0(safe, ".png")), pW, width = width, height = height, dpi = dpi)
        ggsave(file.path(dirW_hz, paste0(safe, ".pdf")), pW, width = width, height = height)
        
        n_saved_W_alldonors <- n_saved_W_alldonors + 1L
      }
    }
    
    if ((n_saved_W_bygroup + n_saved_W_alldonors) == 0L &&
        !is.null(out_dir_W_root) && dir.exists(out_dir_W_root)) {
      try(unlink(out_dir_W_root, recursive = TRUE, force = TRUE), silent = TRUE)
    }
  }
  
  list(
    plots_out_dir = plots_out_dir,
    groups_present = groups_present,
    can_run_B = can_run_B,
    can_week_compare = can_week_compare,
    weeks_in_scope = weeks,
    n_saved_A = n_saved_A,
    n_saved_B = n_saved_B,
    n_saved_W = (n_saved_W_bygroup + n_saved_W_alldonors),
    n_saved_W_bygroup = n_saved_W_bygroup,
    n_saved_W_alldonors = n_saved_W_alldonors
  )
}

# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

build_step3_gui <- function() {
  
  tt <- tktoplevel()
  tkwm.title(tt, "EHT Step 3 — Plotting (v5.1: merge UI only for B/A+B)")
  
  v_csv <- tclVar("")
  v_out <- tclVar("")
  v_keep_follow <- tclVar("1")
  v_plot_mode  <- tclVar("A_only")  # A_only | B_only | A_and_B
  v_week_compare <- tclVar("0")
  
  row <- 0
  
  # CSV
  tkgrid(ttklabel(tt, text="Merged CSV (from Script 02):"),
         row=row, column=0, sticky="w", padx=8, pady=6)
  tkgrid(ttkentry(tt, textvariable=v_csv, width=60),
         row=row, column=1, sticky="we", padx=8, pady=6)
  tkgrid(ttkbutton(tt, text="Browse…", command=function(){
    p <- tk_choose.files(caption="Select merged CSV (from Script 02)")
    if (length(p) && nzchar(p[1])) tclvalue(v_csv) <- p[1]
  }), row=row, column=2, padx=8, pady=6)
  
  row <- row + 1
  
  # Output
  tkgrid(ttklabel(tt, text="Output folder (plots root):"),
         row=row, column=0, sticky="w", padx=8, pady=6)
  tkgrid(ttkentry(tt, textvariable=v_out, width=60),
         row=row, column=1, sticky="we", padx=8, pady=6)
  tkgrid(ttkbutton(tt, text="Browse…", command=function(){
    p <- tk_choose.dir(caption="Select folder to save plots")
    if (!is.na(p) && nzchar(p)) tclvalue(v_out) <- p
  }), row=row, column=2, padx=8, pady=6)
  
  # ---- Week_Time selector ----
  row <- row + 1
  tkgrid(ttklabel(tt, text="Week_Time to include (multi-select):"),
         row=row, column=0, sticky="nw", padx=8, pady=6)
  
  week_frame <- ttkframe(tt)
  tkgrid(week_frame, row=row, column=1, columnspan=2, sticky="we", padx=8, pady=6)
  
  week_list <- tklistbox(week_frame, selectmode="extended", height=6, exportselection=0)
  sbw <- tkscrollbar(week_frame, orient="vertical", command=function(...) tkview(week_list, ...))
  tkconfigure(week_list, yscrollcommand=function(...) tkset(sbw, ...))
  tkgrid(week_list, row=0, column=0, sticky="we")
  tkgrid(sbw, row=0, column=1, sticky="ns")
  tkgrid.columnconfigure(week_frame, 0, weight=1)
  
  weeks_available <- character(0)
  
  get_selected_weeks <- function() {
    sel <- as.integer(tkcurselection(week_list))
    if (!length(sel)) return(character(0))
    weeks_available[sel + 1]
  }
  
  load_weeks_from_csv <- function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE, n_max = 2000)
    if (!("Week_Time" %in% names(df))) stop("Column 'Week_Time' not found in CSV.")
    sort(unique(as.character(df$Week_Time)))
  }
  
  week_btns <- ttkframe(week_frame)
  tkgrid(week_btns, row=1, column=0, columnspan=2, sticky="w", pady=6)
  
  tkgrid(ttkbutton(week_btns, text="Load Week_Time", command=function(){
    csv_path <- tclvalue(v_csv)
    if (!nzchar(csv_path) || !file.exists(csv_path)) {
      tkmessageBox(title="Error", message="Select a valid merged CSV first.", icon="error", type="ok")
      return(invisible(NULL))
    }
    ww <- try(load_weeks_from_csv(csv_path), silent = TRUE)
    if (inherits(ww, "try-error")) {
      tkmessageBox(title="Error", message=as.character(ww), icon="error", type="ok")
      return(invisible(NULL))
    }
    weeks_available <<- ww
    tkdelete(week_list, 0, "end")
    for (w in weeks_available) tkinsert(week_list, "end", w)
    tkselection.set(week_list, 0, "end")
  }), row=0, column=0, padx=6)
  
  tkgrid(ttkbutton(week_btns, text="Select all", command=function(){
    if (!length(weeks_available)) return(invisible(NULL))
    tkselection.set(week_list, 0, "end")
  }), row=0, column=1, padx=6)
  
  tkgrid(ttkbutton(week_btns, text="Select none", command=function(){
    tkselection.clear(week_list, 0, "end")
  }), row=0, column=2, padx=6)
  
  row <- row + 1
  
  # Keep rule
  tkgrid(ttkcheckbutton(tt,
                        text='Keep only: Expected_Hz==0 OR Following=="Yes"',
                        variable=v_keep_follow),
         row=row, column=1, sticky="w", padx=8, pady=6)
  
  row <- row + 1
  
  # Week compare option
  tkgrid(ttkcheckbutton(tt,
                        text="Also create Week comparison plots (only works if exactly 2 Week_Time selected)",
                        variable=v_week_compare),
         row=row, column=1, sticky="w", padx=8, pady=6)
  
  row <- row + 1
  
  # Plot mode
  tkgrid(ttklabel(tt, text="Plot mode:"),
         row=row, column=0, sticky="w", padx=8, pady=6)
  
  mode_frame <- ttkframe(tt)
  tkgrid(mode_frame, row=row, column=1, sticky="w", padx=8, pady=6)
  
  # Donor selector needs donors_all available
  donors_all <- c("C1","C2","P1","P2")
  
  # ---- Merge mapping (created here, shown/hidden based on mode) ----
  # We build widgets now, but only grid them when needed.
  
  # placeholder label row for merge
  # (we keep it always visible as a hint; frame itself is hidden in A_only)
  row <- row + 1
  tkgrid(ttklabel(tt, text="Merge mapping (only for B or A+B):"),
         row=row, column=0, sticky="nw", padx=8, pady=6)
  
  merge_row <- row
  merge_map_frame <- ttkframe(tt)
  
  ctrl_box_frame <- ttkframe(merge_map_frame)
  tkgrid(ctrl_box_frame, row=0, column=0, padx=8, sticky="nw")
  tkgrid(ttklabel(ctrl_box_frame, text="Donors merged into Control:"), row=0, column=0, sticky="w")
  ctrl_list <- tklistbox(ctrl_box_frame, selectmode="extended", height=4, exportselection=0)
  tkgrid(ctrl_list, row=1, column=0, sticky="w")
  for (d in donors_all) tkinsert(ctrl_list, "end", d)
  
  pat_box_frame <- ttkframe(merge_map_frame)
  tkgrid(pat_box_frame, row=0, column=1, padx=8, sticky="nw")
  tkgrid(ttklabel(pat_box_frame, text="Donors merged into Patient:"), row=0, column=0, sticky="w")
  pat_list <- tklistbox(pat_box_frame, selectmode="extended", height=4, exportselection=0)
  tkgrid(pat_list, row=1, column=0, sticky="w")
  for (d in donors_all) tkinsert(pat_list, "end", d)
  
  get_selected_from_list <- function(lb, universe) {
    sel <- as.integer(tkcurselection(lb))
    if (!length(sel)) return(character(0))
    universe[sel + 1]
  }
  
  set_default_merge <- function() {
    tkselection.clear(ctrl_list, 0, "end")
    tkselection.clear(pat_list,  0, "end")
    # Control: C1,C2 ; Patient: P1,P2
    tkselection.set(ctrl_list, 0); tkselection.set(ctrl_list, 1)
    tkselection.set(pat_list,  2); tkselection.set(pat_list,  3)
  }
  
  update_merge_ui <- function() {
    pm <- tclvalue(v_plot_mode)
    show_merge <- pm %in% c("B_only", "A_and_B")
    
    if (show_merge) {
      tkgrid(merge_map_frame, row=merge_row, column=1, columnspan=2, sticky="we", padx=8, pady=6)
      # If nothing selected, set defaults (nice UX)
      if (length(as.integer(tkcurselection(ctrl_list))) == 0 &&
          length(as.integer(tkcurselection(pat_list)))  == 0) {
        set_default_merge()
      }
    } else {
      # hide + clear selections so A_only never uses merge mapping
      tkgrid.remove(merge_map_frame)
      tkselection.clear(ctrl_list, 0, "end")
      tkselection.clear(pat_list,  0, "end")
    }
  }
  
  # Now add plot mode radio buttons (need update_merge_ui callback)
  # We already incremented row for merge label; plot mode controls should be ABOVE donor selector.
  # So we place plot mode controls earlier by using the existing layout row-?:
  # To keep minimal edits, we will place them right now but visually they are after merge label.
  # (If you want them earlier, I can re-order, but this works and is stable.)
  
  tkgrid(ttkradiobutton(mode_frame, text="A only (group plot)", variable=v_plot_mode, value="A_only",
                        command=update_merge_ui),
         row=0, column=0, padx=6)
  tkgrid(ttkradiobutton(mode_frame, text="B only (Control vs Patient)", variable=v_plot_mode, value="B_only",
                        command=update_merge_ui),
         row=0, column=1, padx=6)
  tkgrid(ttkradiobutton(mode_frame, text="A + B", variable=v_plot_mode, value="A_and_B",
                        command=update_merge_ui),
         row=0, column=2, padx=6)
  
  # NOTE: We placed mode_frame earlier (row before merge label). We must grid it now:
  # mode_frame should appear at row = (merge_row - 1)
  # We already created it at the correct earlier row index before merge label; keep as-is.
  
  # ---- Donor selector ----
  row <- row + 1
  tkgrid(ttklabel(tt, text="Donors to include (multi-select):"),
         row=row, column=0, sticky="nw", padx=8, pady=6)
  
  donor_frame <- ttkframe(tt)
  tkgrid(donor_frame, row=row, column=1, columnspan=2, sticky="we", padx=8, pady=6)
  
  donor_list <- tklistbox(donor_frame, selectmode="extended", height=4, exportselection=0)
  tkgrid(donor_list, row=0, column=0, sticky="w")
  for (d in donors_all) tkinsert(donor_list, "end", d)
  tkselection.set(donor_list, 0, "end")
  
  donor_btns <- ttkframe(donor_frame)
  tkgrid(donor_btns, row=0, column=1, sticky="w", padx=10)
  
  tkgrid(ttkbutton(donor_btns, text="Select all", command=function(){
    tkselection.set(donor_list, 0, "end")
  }), row=0, column=0, padx=6)
  
  tkgrid(ttkbutton(donor_btns, text="Select none", command=function(){
    tkselection.clear(donor_list, 0, "end")
  }), row=0, column=1, padx=6)
  
  get_selected_donors <- function() {
    sel <- as.integer(tkcurselection(donor_list))
    if (!length(sel)) return(character(0))
    donors_all[sel + 1]
  }
  
  # Apply initial merge UI visibility based on default plot mode
  update_merge_ui()
  
  # ---- Metrics selector ----
  row <- row + 1
  tkgrid(ttklabel(tt, text="Metrics to plot (multi-select):"),
         row=row, column=0, sticky="nw", padx=8, pady=6)
  
  metrics_frame <- ttkframe(tt)
  tkgrid(metrics_frame, row=row, column=1, columnspan=2, sticky="we", padx=8, pady=6)
  
  metrics_list <- tklistbox(metrics_frame, selectmode="extended", height=10, exportselection=0)
  sb <- tkscrollbar(metrics_frame, orient="vertical", command=function(...) tkview(metrics_list, ...))
  tkconfigure(metrics_list, yscrollcommand=function(...) tkset(sb, ...))
  
  tkgrid(metrics_list, row=0, column=0, sticky="we")
  tkgrid(sb, row=0, column=1, sticky="ns")
  tkgrid.columnconfigure(metrics_frame, 0, weight=1)
  
  metrics_available <- character(0)
  
  get_selected_metrics <- function() {
    sel <- as.integer(tkcurselection(metrics_list))
    if (!length(sel)) return(character(0))
    metrics_available[sel + 1]
  }
  
  load_metrics_from_csv <- function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE, n_max = 50)
    req <- c("Donor","Expected_Hz","Following","Week_Time","Batch")
    miss <- setdiff(req, names(df))
    if (length(miss)) stop("CSV is missing required columns: ", paste(miss, collapse = ", "))
    
    default_metrics <- names(df)[4:18]
    default_metrics <- default_metrics[default_metrics %in% names(df)]
    
    meta <- c("source_file","Batch","Donor","Week_Time","Following","Expected_Hz",
              "Deviation_Hz","Rel_Error","File_Name","Plot_Group")
    cand <- setdiff(names(df), meta)
    
    numeric_like <- cand[sapply(df[cand], function(x) {
      xx <- suppressWarnings(as.numeric(as.character(x)))
      mean(!is.na(xx)) > 0.6
    })]
    
    all_metrics <- unique(c(default_metrics, numeric_like))
    list(all_metrics = all_metrics, default_metrics = default_metrics)
  }
  
  metrics_btns <- ttkframe(metrics_frame)
  tkgrid(metrics_btns, row=1, column=0, columnspan=2, sticky="w", pady=6)
  
  tkgrid(ttkbutton(metrics_btns, text="Load metrics", command=function(){
    csv_path <- tclvalue(v_csv)
    if (!nzchar(csv_path) || !file.exists(csv_path)) {
      tkmessageBox(title="Error", message="Select a valid merged CSV first.", icon="error", type="ok")
      return(invisible(NULL))
    }
    res <- try(load_metrics_from_csv(csv_path), silent = TRUE)
    if (inherits(res, "try-error")) {
      tkmessageBox(title="Error", message=as.character(res), icon="error", type="ok")
      return(invisible(NULL))
    }
    
    metrics_available <<- res$all_metrics
    tkdelete(metrics_list, 0, "end")
    for (m in metrics_available) tkinsert(metrics_list, "end", m)
    
    if (length(res$default_metrics)) {
      tkselection.clear(metrics_list, 0, "end")
      idxs <- which(metrics_available %in% res$default_metrics) - 1
      for (i in idxs) tkselection.set(metrics_list, i)
    }
  }), row=0, column=0, padx=6)
  
  row <- row + 1
  
  # Log
  tkgrid(ttklabel(tt, text="Log:"),
         row=row, column=0, sticky="nw", padx=8, pady=6)
  
  status <- tktext(tt, height=9, width=85)
  tkgrid(status, row=row, column=1, columnspan=2, sticky="we", padx=8, pady=6)
  
  log_line <- function(x) {
    tkinsert(status, "end", paste0(x, "\n"))
    tksee(status, "end")
  }
  
  row <- row + 1
  
  # Buttons
  btn_frame <- ttkframe(tt)
  tkgrid(btn_frame, row=row, column=1, sticky="w", padx=8, pady=10)
  
  run_btn <- ttkbutton(btn_frame, text="Run Step 3")
  tkgrid(run_btn, row=0, column=0, padx=6)
  
  tkconfigure(run_btn, command=function(){
    
    tkconfigure(run_btn, state="disabled")
    
    csv_path <- tclvalue(v_csv)
    out_dir  <- tclvalue(v_out)
    
    if (!nzchar(csv_path) || !file.exists(csv_path)) {
      tkconfigure(run_btn, state="normal")
      tkmessageBox(title="Error", message="Please select a valid merged CSV file.", icon="error", type="ok")
      return(invisible(NULL))
    }
    if (!nzchar(out_dir) || !dir.exists(out_dir)) {
      tkconfigure(run_btn, state="normal")
      tkmessageBox(title="Error", message="Please select a valid output folder.", icon="error", type="ok")
      return(invisible(NULL))
    }
    
    weeks_keep <- get_selected_weeks()
    if (length(weeks_keep) == 0) weeks_keep <- NULL
    
    if (!length(metrics_available)) {
      tkconfigure(run_btn, state="normal")
      tkmessageBox(title="Error", message="Click 'Load metrics' first, then select metrics to plot.", icon="error", type="ok")
      return(invisible(NULL))
    }
    
    metric_cols <- get_selected_metrics()
    if (!length(metric_cols)) {
      tkconfigure(run_btn, state="normal")
      tkmessageBox(title="Error", message="No metrics selected. Select at least one metric.", icon="error", type="ok")
      return(invisible(NULL))
    }
    
    donors_keep <- get_selected_donors()
    if (length(donors_keep) < 1) {
      tkconfigure(run_btn, state="normal")
      tkmessageBox(title="Error", message="Select at least 1 donor.", icon="error", type="ok")
      return(invisible(NULL))
    }
    
    keep_follow <- (tclvalue(v_keep_follow) == "1")
    do_week_compare <- (tclvalue(v_week_compare) == "1")
    plot_mode  <- tclvalue(v_plot_mode)
    
    # Merge mapping is ONLY used for B_only or A_and_B
    ctrl_merge <- character(0)
    pat_merge  <- character(0)
    if (plot_mode %in% c("B_only","A_and_B")) {
      ctrl_merge <- get_selected_from_list(ctrl_list, donors_all)
      pat_merge  <- get_selected_from_list(pat_list, donors_all)
    }
    
    # Make sure merge donors are included
    ctrl_merge <- intersect(ctrl_merge, donors_keep)
    pat_merge  <- intersect(pat_merge, donors_keep)
    
    overlap <- intersect(ctrl_merge, pat_merge)
    if (length(overlap)) {
      tkconfigure(run_btn, state="normal")
      tkmessageBox(
        title="Error",
        message=paste0("Same donor cannot be merged into BOTH groups:\n", paste(overlap, collapse=", ")),
        icon="error", type="ok"
      )
      return(invisible(NULL))
    }
    
    # Internal defaults (GUI removed)
    seed <- 1
    topn <- 6
    
    log_line("------------------------------------------------------------")
    log_line(paste0("CSV: ", csv_path))
    log_line(paste0("Output: ", out_dir))
    log_line(paste0("Week_Time selection: ", if (is.null(weeks_keep)) "(all)" else paste(weeks_keep, collapse = ", ")))
    log_line(paste0('Keep rule (Expected_Hz==0 OR Following=="Yes"): ', keep_follow))
    log_line(paste0("Donors selected: ", paste(donors_keep, collapse = ", ")))
    log_line(paste0("Plot mode: ", plot_mode))
    log_line(paste0("Merge into Control: ", ifelse(length(ctrl_merge), paste(ctrl_merge, collapse=", "), "(not used / none)")))
    log_line(paste0("Merge into Patient: ", ifelse(length(pat_merge), paste(pat_merge, collapse=", "), "(not used / none)")))
    log_line(paste0("Week comparison: ", do_week_compare))
    log_line(paste0("Metrics selected: ", length(metric_cols)))
    
    res <- try(
      run_step3_plots(
        csv_path = csv_path,
        root_dir = out_dir,
        metric_cols = metric_cols,
        weeks_keep = weeks_keep,
        donors_keep = donors_keep,
        control_donors_merge = ctrl_merge,
        patient_donors_merge = pat_merge,
        keep_following_rule = keep_follow,
        jitter_seed = seed,
        plot_mode = plot_mode,
        top_n_tukey = topn,
        do_week_compare = do_week_compare
      ),
      silent = TRUE
    )
    
    if (inherits(res, "try-error")) {
      tkconfigure(run_btn, state="normal")
      tkmessageBox(title="Error", message=as.character(res), icon="error", type="ok")
      log_line(paste0("ERROR: ", as.character(res)))
      return(invisible(NULL))
    }
    
    log_line("✅ Done.")
    log_line(paste0("Groups present: ", paste(res$groups_present, collapse = ", ")))
    log_line(paste0("Weeks in scope: ", paste(res$weeks_in_scope, collapse = ", ")))
    log_line(paste0("Saved A plots: ", res$n_saved_A))
    log_line(paste0("Saved B plots: ", res$n_saved_B))
    log_line(paste0("Saved WeekCompare plots: ", res$n_saved_W,
                    " (BY_GROUP=", res$n_saved_W_bygroup,
                    ", ALL_DONORS=", res$n_saved_W_alldonors, ")"))
    
    if (plot_mode %in% c("B_only","A_and_B") && !res$can_run_B) {
      log_line("Note: B plots skipped (requires ONLY two groups after merging: Control + Patient).")
    }
    if (do_week_compare && !res$can_week_compare) {
      log_line("Note: Week comparison skipped (requires exactly 2 Week_Time values in scope).")
    }
    
    tkmessageBox(
      title = "Done",
      message = paste0(
        "Step 3 completed successfully.\n\nPlots saved in:\n",
        res$plots_out_dir,
        "\n\nSaved A plots: ", res$n_saved_A,
        "\nSaved B plots: ", res$n_saved_B,
        "\nSaved WeekCompare plots: ", res$n_saved_W,
        " (BY_GROUP=", res$n_saved_W_bygroup,
        ", ALL_DONORS=", res$n_saved_W_alldonors, ")",
        "\n\nYou can close this window now."
      ),
      icon = "info",
      type = "ok"
    )
    
    tkconfigure(run_btn, state="normal")
  })
  
  tkgrid(ttkbutton(btn_frame, text="Close", command=function() tkdestroy(tt)),
         row=0, column=1, padx=6)
  
  tkgrid.columnconfigure(tt, 1, weight=1)
}

build_step3_gui()
