# import packages
packages <- c(
  "tidyverse",
  "cluster",
  "factoextra",
  "arules",
  "arulesViz",
  "e1071"
)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

set.seed(123)

# import dataset
candidate_files <- c(
  "C:/Users/HP/Downloads/heart_disease_uci.csv"
)

file_path <- candidate_files[file.exists(candidate_files)][1]

if(is.na(file_path)){
  message("File tidak ditemukan")
  file_path <- file.choose()
}

# read datasett
df_raw <- read.csv(
  file_path,
  na.strings = c("", " ", "NA", "?", "NaN"),
  stringsAsFactors = FALSE
)

cat("\nBagian import dataset\n")
cat("Dimensi data awal:", nrow(df_raw), "baris dan", ncol(df_raw), "kolom\n")

names(df_raw) <- names(df_raw) %>%
  tolower() %>%
  gsub("[^a-z0-9]+", "_", .) %>%
  gsub("_$", "", .)

if("thalach" %in% names(df_raw) && !"thalch" %in% names(df_raw)){
  df_raw <- df_raw %>%
    dplyr::rename(thalch = thalach)
}

cat("\nNama kolom yang terbaca\n")
print(names(df_raw))

cat("\nMissing value awal\n")
print(colSums(is.na(df_raw)))

# definisi variabel
# nums var
num_vars <- c(
  "age",
  "trestbps",
  "chol",
  "thalch",
  "oldpeak"
)

# cats var
cat_vars <- c(
  "sex",
  "cp",
  "fbs",
  "restecg",
  "exang",
  "slope",
  "ca",
  "thal"
)

target_var <- "num" # target var

needed_vars <- c(num_vars, cat_vars, target_var)
missing_vars <- setdiff(needed_vars, names(df_raw))

if(length(missing_vars) > 0){
  stop(
    paste(
      "Kolom berikut tidak ditemukan:",
      paste(missing_vars, collapse = ", ")
    )
  )
}

# helper function
mode_value <- function(x){
  x <- na.omit(x)
  if(length(x) == 0){
    return("Unknown")
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# function scaling
rescale01 <- function(x){
  if(all(is.na(x))){
    return(rep(0.5, length(x)))
  }
  
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  
  if(isTRUE(all.equal(min_x, max_x))){
    return(rep(0.5, length(x)))
  }
  
  (x - min_x) / (max_x - min_x)
}

# function silh
safe_silhouette_mean <- function(labels, dist_obj){
  labels <- as.integer(as.factor(labels))
  
  if(length(unique(labels)) < 2){
    return(NA_real_)
  }
  
  sil <- cluster::silhouette(labels, dist_obj)
  mean(sil[, 3])
}

# function dunn 
dunn_index_from_dist <- function(labels, dist_obj){
  labels <- as.factor(labels)
  dm <- as.matrix(dist_obj)
  cluster_ids <- levels(labels)
  
  if(length(cluster_ids) < 2){
    return(NA_real_)
  }
  
  diameters <- c()
  
  for(cl in cluster_ids){
    idx <- which(labels == cl)
    
    if(length(idx) <= 1){
      diameters <- c(diameters, 0)
    } else {
      diameters <- c(diameters, max(dm[idx, idx, drop = FALSE]))
    }
  }
  
  max_diameter <- max(diameters, na.rm = TRUE)
  
  inter_distances <- c()
  
  for(i in seq_along(cluster_ids)){
    for(j in seq_along(cluster_ids)){
      if(j > i){
        idx_i <- which(labels == cluster_ids[i])
        idx_j <- which(labels == cluster_ids[j])
        inter_distances <- c(
          inter_distances,
          min(dm[idx_i, idx_j, drop = FALSE], na.rm = TRUE)
        )
      }
    }
  }
  
  if(max_diameter == 0){
    return(NA_real_)
  }
  
  min(inter_distances, na.rm = TRUE) / max_diameter
}

# function calc rasio rata-rata jarak intra-cluster
within_between_ratio_from_dist <- function(labels, dist_obj){
  labels <- as.factor(labels)
  dm <- as.matrix(dist_obj)
  cluster_ids <- levels(labels)
  
  within_values <- c()
  between_values <- c()
  
  for(cl in cluster_ids){
    idx <- which(labels == cl)
    
    if(length(idx) > 1){
      temp <- dm[idx, idx, drop = FALSE]
      within_values <- c(within_values, temp[upper.tri(temp)])
    }
  }
  
  for(i in seq_along(cluster_ids)){
    for(j in seq_along(cluster_ids)){
      if(j > i){
        idx_i <- which(labels == cluster_ids[i])
        idx_j <- which(labels == cluster_ids[j])
        between_values <- c(
          between_values,
          as.vector(dm[idx_i, idx_j, drop = FALSE])
        )
      }
    }
  }
  
  if(length(within_values) == 0 || length(between_values) == 0){
    return(NA_real_)
  }
  
  mean(within_values, na.rm = TRUE) / mean(between_values, na.rm = TRUE)
}

# function score" lain (purity, ari dll)
purity_score <- function(cluster_label, true_label){
  tab <- table(cluster_label, true_label)
  sum(apply(tab, 1, max)) / sum(tab)
}

ari_score <- function(cluster_label, true_label){
  tab <- table(cluster_label, true_label)
  
  comb2 <- function(x){
    ifelse(x < 2, 0, x * (x - 1) / 2)
  }
  
  sum_comb <- sum(comb2(tab))
  row_comb <- sum(comb2(rowSums(tab)))
  col_comb <- sum(comb2(colSums(tab)))
  total_comb <- comb2(sum(tab))
  
  expected_index <- row_comb * col_comb / total_comb
  max_index <- 0.5 * (row_comb + col_comb)
  
  if(isTRUE(all.equal(max_index, expected_index))){
    return(0)
  }
  
  (sum_comb - expected_index) / (max_index - expected_index)
}

nmi_score <- function(cluster_label, true_label){
  tab <- table(cluster_label, true_label)
  n <- sum(tab)
  
  p_ij <- tab / n
  p_i <- rowSums(p_ij)
  p_j <- colSums(p_ij)
  
  expected <- outer(p_i, p_j)
  valid <- p_ij > 0 & expected > 0
  
  mi <- sum(p_ij[valid] * log(p_ij[valid] / expected[valid]))
  h_i <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
  h_j <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))
  
  if(h_i == 0 || h_j == 0){
    return(0)
  }
  
  mi / sqrt(h_i * h_j)
}

db_index_numeric <- function(data_scaled, labels){
  labels <- as.factor(labels)
  cluster_ids <- levels(labels)
  
  if(length(cluster_ids) < 2){
    return(NA_real_)
  }
  
  centers <- matrix(
    NA_real_,
    nrow = length(cluster_ids),
    ncol = ncol(data_scaled)
  )
  
  s_values <- c()
  
  for(i in seq_along(cluster_ids)){
    idx <- which(labels == cluster_ids[i])
    centers[i, ] <- colMeans(data_scaled[idx, , drop = FALSE])
    
    if(length(idx) <= 1){
      s_values <- c(s_values, 0)
    } else {
      distances <- sqrt(rowSums(
        sweep(data_scaled[idx, , drop = FALSE], 2, centers[i, ])^2
      ))
      s_values <- c(s_values, mean(distances))
    }
  }
  
  center_dist <- as.matrix(dist(centers))
  db_values <- c()
  
  for(i in seq_along(cluster_ids)){
    ratio_values <- c()
    
    for(j in seq_along(cluster_ids)){
      if(i != j){
        if(center_dist[i, j] == 0){
          ratio_values <- c(ratio_values, NA_real_)
        } else {
          ratio_values <- c(
            ratio_values,
            (s_values[i] + s_values[j]) / center_dist[i, j]
          )
        }
      }
    }
    
    db_values <- c(db_values, max(ratio_values, na.rm = TRUE))
  }
  
  mean(db_values, na.rm = TRUE)
}

partition_metrics <- function(membership_matrix){
  U <- as.matrix(membership_matrix)
  cluster_count <- ncol(U)
  
  pc <- sum(U^2) / nrow(U)
  pe <- -sum(U * log(U + 1e-10)) / nrow(U)
  
  if(cluster_count <= 1){
    mpc <- NA_real_
  } else {
    mpc <- (pc - (1 / cluster_count)) / (1 - (1 / cluster_count))
  }
  
  data.frame(
    partition_coefficient = pc,
    partition_entropy = pe,
    modified_partition_coefficient = mpc
  )
}

xie_beni_index <- function(data_scaled, centers, membership_matrix, m_value){
  U <- as.matrix(membership_matrix)
  centers <- as.matrix(centers)
  
  dist_sq <- matrix(
    NA_real_,
    nrow = nrow(data_scaled),
    ncol = nrow(centers)
  )
  
  for(j in 1:nrow(centers)){
    dist_sq[, j] <- rowSums(
      sweep(data_scaled, 2, centers[j, ])^2
    )
  }
  
  numerator <- sum((U^m_value) * dist_sq)
  
  center_dist_sq <- as.matrix(dist(centers))^2
  min_center_dist_sq <- min(center_dist_sq[upper.tri(center_dist_sq)], na.rm = TRUE)
  
  if(min_center_dist_sq == 0){
    return(NA_real_)
  }
  
  numerator / (nrow(data_scaled) * min_center_dist_sq)
}

make_bin <- function(x, label_name){
  x <- as.numeric(x)
  
  q <- unique(
    quantile(
      x,
      probs = c(0, 1/3, 2/3, 1),
      na.rm = TRUE
    )
  )
  
  if(length(q) < 4){
    med <- median(x, na.rm = TRUE)
    
    return(
      factor(
        ifelse(
          x <= med,
          paste0(label_name, "_Low"),
          paste0(label_name, "_High")
        )
      )
    )
  }
  
  cut(
    x,
    breaks = q,
    labels = c(
      paste0(label_name, "_Low"),
      paste0(label_name, "_Medium"),
      paste0(label_name, "_High")
    ),
    include.lowest = TRUE
  )
}

clean_item <- function(x, prefix){
  x <- tolower(as.character(x))
  x <- stringr::str_squish(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_$", "", x)
  paste0(prefix, "_", x)
}

rules_to_df <- function(rules_obj, transactions_obj = NULL){
  if(length(rules_obj) == 0){
    return(
      data.frame(
        rules = character(),
        support = numeric(),
        confidence = numeric(),
        coverage = numeric(),
        lift = numeric(),
        count = numeric(),
        leverage = numeric(),
        conviction = numeric()
      )
    )
  }
  
  df_rules <- as(rules_obj, "data.frame")
  
  if(!"count" %in% names(df_rules)){
    df_rules$count <- round(df_rules$support * length(transactions_obj))
  }
  
  if(!is.null(transactions_obj)){
    extra_quality <- arules::interestMeasure(
      rules_obj,
      measure = c("leverage", "conviction"),
      transactions = transactions_obj
    )
    
    df_rules$leverage <- extra_quality$leverage
    df_rules$conviction <- extra_quality$conviction
  }
  
  if(!"leverage" %in% names(df_rules)){
    df_rules$leverage <- NA_real_
  }
  
  if(!"conviction" %in% names(df_rules)){
    df_rules$conviction <- NA_real_
  }
  
  df_rules %>%
    dplyr::mutate(
      support = round(as.numeric(support), 4),
      confidence = round(as.numeric(confidence), 4),
      coverage = round(as.numeric(coverage), 4),
      lift = round(as.numeric(lift), 4),
      count = as.numeric(count),
      leverage = round(as.numeric(leverage), 4),
      conviction = round(as.numeric(conviction), 4)
    )
}

# function rules apriori
inspect_rules_safe <- function(rules_obj, title, transactions_obj = NULL, n = 10){
  cat("\n")
  cat(title, "\n")
  
  if(length(rules_obj) == 0){
    cat("Tidak ada rules yang terbentuk.\n")
  } else {
    print(head(rules_to_df(rules_obj, transactions_obj), n))
  }
}

profil_cluster <- function(data_input, cluster_col){
  data_input %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::summarise(
      jumlah_pasien = dplyr::n(),
      age_mean = round(mean(age, na.rm = TRUE), 2),
      trestbps_mean = round(mean(trestbps, na.rm = TRUE), 2),
      chol_mean = round(mean(chol, na.rm = TRUE), 2),
      thalch_mean = round(mean(thalch, na.rm = TRUE), 2),
      oldpeak_mean = round(mean(oldpeak, na.rm = TRUE), 2),
      persen_disease = round(mean(disease_bin == "Disease", na.rm = TRUE) * 100, 2),
      .groups = "drop"
    )
}

buat_data_apriori <- function(data_input, cluster_col){
  data_input %>%
    dplyr::transmute(
      Age = make_bin(age, "Age"),
      RestingBP = make_bin(trestbps, "RestingBP"),
      Cholesterol = make_bin(chol, "Cholesterol"),
      MaxHeartRate = make_bin(thalch, "MaxHeartRate"),
      Oldpeak = make_bin(oldpeak, "Oldpeak"),
      Sex = factor(clean_item(sex, "Sex")),
      ChestPain = factor(clean_item(cp, "ChestPain")),
      FastingBloodSugar = factor(clean_item(fbs, "FBS")),
      RestECG = factor(clean_item(restecg, "RestECG")),
      ExerciseAngina = factor(clean_item(exang, "ExerciseAngina")),
      Slope = factor(clean_item(slope, "Slope")),
      CA = factor(clean_item(ca, "CA")),
      Thal = factor(clean_item(thal, "Thal")),
      Disease = disease_bin,
      Cluster = .data[[cluster_col]]
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::everything(), as.factor)
    )
}

# Preprocessing dataset
cat("\nBagian preprocessing data\n")

df <- df_raw %>%
  dplyr::select(-dplyr::any_of(c("id", "dataset"))) %>%
  dplyr::mutate(
    dplyr::across(
      where(is.character),
      ~ stringr::str_squish(.x)
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(num_vars), ~ as.numeric(.x)),
    num = as.numeric(num)
  ) %>%
  dplyr::mutate(
    trestbps = ifelse(!is.na(trestbps) & trestbps <= 0, NA_real_, trestbps),
    chol = ifelse(!is.na(chol) & chol <= 0, NA_real_, chol)
  ) %>%
  dplyr::filter(!is.na(num)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(num_vars),
      ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(cat_vars),
      ~ factor(ifelse(is.na(.x) | .x == "", "Unknown", as.character(.x)))
    )
  ) %>%
  dplyr::mutate(
    disease_bin = ifelse(num == 0, "No_Disease", "Disease"),
    disease_bin = factor(disease_bin, levels = c("No_Disease", "Disease"))
  )

cat("Dimensi data setelah preprocessing:", nrow(df), "baris dan", ncol(df), "kolom\n")

cat("\nMissing value setelah preprocessing\n")
print(colSums(is.na(df)))

cat("\nDistribusi status penyakit\n")
print(table(df$disease_bin))
print(round(prop.table(table(df$disease_bin)) * 100, 2))

cat("\nJumlah nilai Unknown pada variabel kategorikal\n")
unknown_summary <- sapply(cat_vars, function(v){
  sum(as.character(df[[v]]) == "Unknown")
})
print(unknown_summary)

# Hierarchical Clustering
cat("\nBagian Hierarchical Clustering (Baseline)\n")

hc_awal_features <- c(
  "age",
  "sex",
  "cp",
  "trestbps",
  "chol",
  "fbs",
  "restecg",
  "thalch",
  "exang",
  "oldpeak",
  "slope",
  "ca",
  "thal"
)

data_hc_awal <- df %>%
  dplyr::select(dplyr::all_of(hc_awal_features))

dist_hc_awal <- cluster::daisy(data_hc_awal, metric = "gower")
hc_awal <- hclust(dist_hc_awal, method = "average")

k_awal <- 2
label_hc_awal <- cutree(hc_awal, k = k_awal)

df$cluster_hc_awal <- factor(paste0("HC_Awal_", label_hc_awal))

sil_hc_awal <- cluster::silhouette(label_hc_awal, dist_hc_awal)
sil_hc_awal_mean <- mean(sil_hc_awal[, 3])
dunn_hc_awal <- dunn_index_from_dist(label_hc_awal, dist_hc_awal)
ratio_hc_awal <- within_between_ratio_from_dist(label_hc_awal, dist_hc_awal)
coph_hc_awal <- cor(as.vector(dist_hc_awal), as.vector(cophenetic(hc_awal)))

hc_awal_external <- data.frame(
  purity = purity_score(df$cluster_hc_awal, df$disease_bin),
  ari = ari_score(df$cluster_hc_awal, df$disease_bin),
  nmi = nmi_score(df$cluster_hc_awal, df$disease_bin)
)

cat("Silhouette HC awal:", round(sil_hc_awal_mean, 4), "\n")
cat("Dunn Index HC awal:", round(dunn_hc_awal, 4), "\n")
cat("Within between ratio HC awal:", round(ratio_hc_awal, 4), "\n")
cat("Cophenetic correlation HC awal:", round(coph_hc_awal, 4), "\n")

cat("\nUkuran cluster HC awal\n")
print(table(df$cluster_hc_awal))

cat("\nEvaluasi eksternal HC awal terhadap label disease\n")
print(round(hc_awal_external, 4))

cat("\nTabulasi HC awal terhadap disease\n")
tab_hc_awal <- table(df$cluster_hc_awal, df$disease_bin)
print(tab_hc_awal)
print(round(prop.table(tab_hc_awal, margin = 1) * 100, 2))

cat("\nProfil cluster HC awal\n")
profil_hc_awal <- profil_cluster(df, "cluster_hc_awal")
print(profil_hc_awal)

# Hierarchical Clustering model (opttimized)
cat("\nBagian Hierarchical Clustering model optimized\n")

# seleksi fitur
hc_feature_sets <- list(
  list(
    name = "campuran_lengkap",
    type = "mixed",
    features = c(
      "age", "sex", "cp", "trestbps", "chol", "fbs",
      "restecg", "thalch", "exang", "oldpeak",
      "slope", "ca", "thal"
    )
  ),
  list(
    name = "campuran_tanpa_missing",
    type = "mixed",
    features = c(
      "age", "sex", "cp", "trestbps", "chol",
      "fbs", "restecg", "thalch", "exang", "oldpeak"
    )
  ),
  list(
    name = "campuran_risiko_inti",
    type = "mixed",
    features = c(
      "age", "sex", "cp", "thalch", "exang", "oldpeak"
    )
  ),
  list(
    name = "numerik_lengkap",
    type = "numeric",
    features = c(
      "age", "trestbps", "chol", "thalch", "oldpeak"
    )
  ),
  list(
    name = "numerik_tanpa_chol",
    type = "numeric",
    features = c(
      "age", "trestbps", "thalch", "oldpeak"
    )
  ),
  list(
    name = "numerik_risiko_inti",
    type = "numeric",
    features = c(
      "age", "thalch", "oldpeak"
    )
  )
)

hc_compare <- data.frame()
hc_objects <- list()

for(feature_config in hc_feature_sets){
  fs_name <- feature_config$name
  fs_type <- feature_config$type
  fitur_pakai <- feature_config$features
  
  data_temp <- df %>%
    dplyr::select(dplyr::all_of(fitur_pakai))
  
  if(fs_type == "numeric"){
    data_model <- scale(data_temp)
    dist_temp <- dist(data_model, method = "euclidean")
    methods_use <- c("average", "complete", "mcquitty", "ward.D2")
  } else {
    data_model <- data_temp
    dist_temp <- cluster::daisy(data_model, metric = "gower")
    methods_use <- c("average", "complete", "mcquitty", "single")
  }
  
  for(method_name in methods_use){
    hc_temp <- hclust(dist_temp, method = method_name)
    
    for(k in 2:8){
      label_temp <- cutree(hc_temp, k = k)
      cluster_size_temp <- table(label_temp)
      min_cluster_size <- min(cluster_size_temp)
      max_cluster_size <- max(cluster_size_temp)
      balance_ratio <- min_cluster_size / max_cluster_size
      
      sil_temp <- safe_silhouette_mean(label_temp, dist_temp)
      dunn_temp <- dunn_index_from_dist(label_temp, dist_temp)
      ratio_temp <- within_between_ratio_from_dist(label_temp, dist_temp)
      coph_temp <- cor(as.vector(dist_temp), as.vector(cophenetic(hc_temp)))
      
      key_name <- paste(fs_name, fs_type, method_name, k, sep = "_")
      
      hc_compare <- rbind(
        hc_compare,
        data.frame(
          key = key_name,
          feature_set = fs_name,
          feature_type = fs_type,
          method = method_name,
          k = k,
          avg_silhouette = sil_temp,
          dunn_index = dunn_temp,
          within_between_ratio = ratio_temp,
          cophenetic_correlation = coph_temp,
          min_cluster_size = as.integer(min_cluster_size),
          balance_ratio = as.numeric(balance_ratio),
          cluster_size = paste(
            names(cluster_size_temp),
            as.integer(cluster_size_temp),
            sep = ":",
            collapse = "; "
          )
        )
      )
      
      hc_objects[[key_name]] <- list(
        fitur = fitur_pakai,
        type = fs_type,
        dist = dist_temp,
        hc_model = hc_temp,
        label = label_temp,
        data_model = data_model
      )
    }
  }
}

hc_compare <- hc_compare %>%
  dplyr::mutate(
    dunn_score = rescale01(dunn_index),
    ratio_score = 1 - rescale01(within_between_ratio),
    coph_score = rescale01(cophenetic_correlation),
    score = avg_silhouette +
      0.15 * dunn_score +
      0.10 * ratio_score +
      0.05 * coph_score +
      0.05 * balance_ratio
  ) %>%
  dplyr::arrange(dplyr::desc(score))

min_size_rule <- max(10, ceiling(0.03 * nrow(df)))

hc_compare_filtered <- hc_compare %>%
  dplyr::filter(
    min_cluster_size >= min_size_rule,
    !is.na(avg_silhouette)
  ) %>%
  dplyr::arrange(dplyr::desc(score))

cat("\nRanking kandidat HC terbaik sebelum filter\n")
print(head(hc_compare, 15))

cat("\nRanking kandidat HC terbaik setelah filter ukuran cluster\n")
print(head(hc_compare_filtered, 15))

if(nrow(hc_compare_filtered) == 0){
  best_hc <- hc_compare[1, ]
} else {
  best_hc <- hc_compare_filtered[1, ]
}

cat("\nModel HC optimal yang dipilih\n")
print(best_hc)

hc_opt_obj <- hc_objects[[best_hc$key]]
label_hc_opt <- hc_opt_obj$label
dist_hc_opt <- hc_opt_obj$dist
hc_opt <- hc_opt_obj$hc_model

df$cluster_hc_opt <- factor(paste0("HC_Opt_", label_hc_opt))

sil_hc_opt <- cluster::silhouette(label_hc_opt, dist_hc_opt)
sil_hc_opt_mean <- mean(sil_hc_opt[, 3])
dunn_hc_opt <- dunn_index_from_dist(label_hc_opt, dist_hc_opt)
ratio_hc_opt <- within_between_ratio_from_dist(label_hc_opt, dist_hc_opt)
coph_hc_opt <- cor(as.vector(dist_hc_opt), as.vector(cophenetic(hc_opt)))

hc_opt_external <- data.frame(
  purity = purity_score(df$cluster_hc_opt, df$disease_bin),
  ari = ari_score(df$cluster_hc_opt, df$disease_bin),
  nmi = nmi_score(df$cluster_hc_opt, df$disease_bin)
)

cat("\nSilhouette HC optimal:", round(sil_hc_opt_mean, 4), "\n")
cat("Dunn Index HC optimal:", round(dunn_hc_opt, 4), "\n")
cat("Within between ratio HC optimal:", round(ratio_hc_opt, 4), "\n")
cat("Cophenetic correlation HC optimal:", round(coph_hc_opt, 4), "\n")

cat("\nUkuran cluster HC optimal\n")
print(table(df$cluster_hc_opt))

cat("\nEvaluasi eksternal HC optimal terhadap label disease\n")
print(round(hc_opt_external, 4))

cat("\nTabulasi HC optimal terhadap disease\n")
tab_hc_opt <- table(df$cluster_hc_opt, df$disease_bin)
print(tab_hc_opt)
print(round(prop.table(tab_hc_opt, margin = 1) * 100, 2))

cat("\nProfil cluster HC optimal\n")
profil_hc_opt <- profil_cluster(df, "cluster_hc_opt")
print(profil_hc_opt)

plot(
  hc_opt,
  labels = FALSE,
  main = paste0("Dendrogram HC optimal dengan k ", best_hc$k),
  xlab = "Observasi",
  ylab = "Jarak"
)

rect.hclust(
  hc_opt,
  k = best_hc$k,
  border = 2:(best_hc$k + 1)
)

factoextra::fviz_silhouette(
  sil_hc_opt,
  main = "Silhouette Plot HC Optimal"
)

# Apriori (Baseline)
cat("\nBagian Apriori baseline\n")

df_ap_awal <- buat_data_apriori(
  data_input = df,
  cluster_col = "cluster_hc_awal"
)

trans_awal <- as(df_ap_awal, "transactions")

cat("\nRingkasan transaksi Apriori awal\n")
print(summary(trans_awal))

rules_awal <- arules::apriori(
  trans_awal,
  parameter = list(
    supp = 0.10,
    conf = 0.70,
    minlen = 2,
    maxlen = 5
  ),
  control = list(verbose = FALSE)
)

rules_awal <- sort(rules_awal, by = "lift", decreasing = TRUE)

if(length(rules_awal) > 0){
  rules_awal_nonred <- rules_awal[!is.redundant(rules_awal)]
} else {
  rules_awal_nonred <- rules_awal
}

rules_awal_disease <- subset(
  rules_awal_nonred,
  rhs %in% "Disease=Disease"
)

rules_awal_no_disease <- subset(
  rules_awal_nonred,
  rhs %in% "Disease=No_Disease"
)

cat("\nJumlah rules Apriori awal:", length(rules_awal), "\n")
cat("Jumlah rules non redundant Apriori awal:", length(rules_awal_nonred), "\n")
cat("Jumlah rules menuju Disease:", length(rules_awal_disease), "\n")
cat("Jumlah rules menuju No Disease:", length(rules_awal_no_disease), "\n")

inspect_rules_safe(
  rules_awal_nonred,
  "Top rules Apriori awal berdasarkan lift",
  trans_awal,
  n = 10
)

inspect_rules_safe(
  sort(rules_awal_disease, by = "lift", decreasing = TRUE),
  "Top rules Apriori awal menuju Disease",
  trans_awal,
  n = 10
)

# Apriori model (Optimized)
cat("\nBagian Apriori model optimized\n")

df_ap_opt <- buat_data_apriori(
  data_input = df,
  cluster_col = "cluster_hc_opt"
)

trans_opt <- as(df_ap_opt, "transactions")

support_values <- c(0.10, 0.07, 0.05, 0.03, 0.02, 0.01)
confidence_values <- c(0.80, 0.70, 0.60, 0.50)

apriori_compare <- data.frame()
apriori_objects <- list()

for(supp_val in support_values){
  for(conf_val in confidence_values){
    rules_temp <- arules::apriori(
      trans_opt,
      parameter = list(
        supp = supp_val,
        conf = conf_val,
        minlen = 2,
        maxlen = 6
      ),
      control = list(verbose = FALSE)
    )
    
    if(length(rules_temp) > 0){
      rules_temp <- sort(rules_temp, by = "lift", decreasing = TRUE)
      rules_temp <- rules_temp[!is.redundant(rules_temp)]
    }
    
    rules_disease_temp <- subset(
      rules_temp,
      rhs %in% "Disease=Disease"
    )
    
    rules_no_disease_temp <- subset(
      rules_temp,
      rhs %in% "Disease=No_Disease"
    )
    
    rules_cluster_temp <- subset(
      rules_temp,
      rhs %pin% "Cluster="
    )
    
    top_lift <- ifelse(
      length(rules_temp) > 0,
      mean(head(quality(rules_temp)$lift, 10), na.rm = TRUE),
      0
    )
    
    disease_lift <- ifelse(
      length(rules_disease_temp) > 0,
      mean(head(quality(sort(rules_disease_temp, by = "lift", decreasing = TRUE))$lift, 10), na.rm = TRUE),
      0
    )
    
    rule_count <- length(rules_temp)
    disease_count <- length(rules_disease_temp)
    no_disease_count <- length(rules_no_disease_temp)
    cluster_count <- length(rules_cluster_temp)
    
    score_temp <- top_lift +
      0.50 * disease_lift +
      0.10 * log1p(rule_count) +
      0.20 * log1p(disease_count) +
      0.05 * log1p(cluster_count)
    
    if(rule_count > 5000){
      score_temp <- score_temp * 0.75
    }
    
    key_name <- paste0(
      "supp_", supp_val,
      "_conf_", conf_val
    )
    
    apriori_compare <- rbind(
      apriori_compare,
      data.frame(
        key = key_name,
        support = supp_val,
        confidence = conf_val,
        jumlah_rules = rule_count,
        rules_disease = disease_count,
        rules_no_disease = no_disease_count,
        rules_cluster = cluster_count,
        rata_lift_top10 = top_lift,
        rata_lift_disease_top10 = disease_lift,
        score = score_temp
      )
    )
    
    apriori_objects[[key_name]] <- rules_temp
  }
}

apriori_compare <- apriori_compare %>%
  dplyr::arrange(dplyr::desc(score))

cat("\nRanking kandidat Apriori optimalisasi\n")
print(apriori_compare)

apriori_compare_filtered <- apriori_compare %>%
  dplyr::filter(
    jumlah_rules >= 10,
    rules_disease >= 1
  ) %>%
  dplyr::arrange(dplyr::desc(score))

if(nrow(apriori_compare_filtered) == 0){
  best_apriori <- apriori_compare[1, ]
} else {
  best_apriori <- apriori_compare_filtered[1, ]
}

cat("\nParameter Apriori optimal yang dipilih\n")
print(best_apriori)

rules_opt <- apriori_objects[[best_apriori$key]]
rules_opt <- sort(rules_opt, by = "lift", decreasing = TRUE)

rules_opt_disease <- subset(
  rules_opt,
  rhs %in% "Disease=Disease"
)

rules_opt_no_disease <- subset(
  rules_opt,
  rhs %in% "Disease=No_Disease"
)

rules_opt_cluster <- subset(
  rules_opt,
  rhs %pin% "Cluster="
)

inspect_rules_safe(
  rules_opt,
  "Top rules Apriori optimal berdasarkan lift",
  trans_opt,
  n = 10
)

inspect_rules_safe(
  sort(rules_opt_disease, by = "lift", decreasing = TRUE),
  "Top rules Apriori optimal menuju Disease",
  trans_opt,
  n = 10
)

inspect_rules_safe(
  sort(rules_opt_no_disease, by = "lift", decreasing = TRUE),
  "Top rules Apriori optimal menuju No Disease",
  trans_opt,
  n = 10
)

inspect_rules_safe(
  sort(rules_opt_cluster, by = "lift", decreasing = TRUE),
  "Top rules Apriori optimal menuju Cluster",
  trans_opt,
  n = 10
)

prop_disease_cluster_opt <- prop.table(
  table(df$cluster_hc_opt, df$disease_bin),
  margin = 1
)

risk_cluster_opt <- rownames(prop_disease_cluster_opt)[
  which.max(prop_disease_cluster_opt[, "Disease"])
]

risk_cluster_item <- paste0("Cluster=", risk_cluster_opt)

cat("\nCluster HC optimal dengan proporsi Disease paling tinggi:", risk_cluster_opt, "\n")

rules_risk_cluster <- rules_opt[0]

if(risk_cluster_item %in% itemLabels(trans_opt)){
  rules_risk_cluster <- arules::apriori(
    trans_opt,
    parameter = list(
      supp = 0.005,
      conf = 0.30,
      minlen = 2,
      maxlen = 5
    ),
    appearance = list(
      rhs = risk_cluster_item,
      default = "lhs"
    ),
    control = list(verbose = FALSE)
  )
  
  if(length(rules_risk_cluster) > 0){
    rules_risk_cluster <- sort(rules_risk_cluster, by = "lift", decreasing = TRUE)
    rules_risk_cluster <- rules_risk_cluster[!is.redundant(rules_risk_cluster)]
  }
}

inspect_rules_safe(
  rules_risk_cluster,
  paste0("Top rules menuju cluster risiko ", risk_cluster_opt),
  trans_opt,
  n = 10
)

if(length(rules_opt) > 0){
  plot(
    rules_opt,
    method = "scatterplot",
    measure = c("support", "confidence"),
    shading = "lift",
    main = "Apriori optimal support confidence dan lift"
  )
}

if(length(rules_opt) >= 2){
  plot(
    head(rules_opt, 30),
    method = "grouped",
    main = "Apriori optimal top rules"
  )
}

# Fuzzy C Means (Baseline)
cat("\nBagian Fuzzy C Means baseline\n")

fcm_awal_features <- c(
  "age",
  "trestbps",
  "chol",
  "thalch",
  "oldpeak"
)

data_fcm_awal <- df %>%
  dplyr::select(dplyr::all_of(fcm_awal_features))

data_fcm_awal_scaled <- scale(data_fcm_awal)
dist_fcm_awal <- dist(data_fcm_awal_scaled, method = "euclidean")

c_fcm_awal <- 2
m_fcm_awal <- 2

set.seed(123)

fcm_awal <- e1071::cmeans(
  data_fcm_awal_scaled,
  centers = c_fcm_awal,
  iter.max = 300,
  m = m_fcm_awal,
  method = "cmeans",
  dist = "euclidean"
)

label_fcm_awal <- fcm_awal$cluster
membership_fcm_awal <- as.matrix(fcm_awal$membership)
max_memb_awal <- apply(membership_fcm_awal, 1, max)

sil_fcm_awal <- cluster::silhouette(label_fcm_awal, dist_fcm_awal)
sil_fcm_awal_mean <- mean(sil_fcm_awal[, 3])
dunn_fcm_awal <- dunn_index_from_dist(label_fcm_awal, dist_fcm_awal)
dbi_fcm_awal <- db_index_numeric(data_fcm_awal_scaled, label_fcm_awal)
pc_pe_awal <- partition_metrics(membership_fcm_awal)
xb_fcm_awal <- xie_beni_index(
  data_fcm_awal_scaled,
  fcm_awal$centers,
  membership_fcm_awal,
  m_fcm_awal
)

df$cluster_fcm_awal <- factor(paste0("FCM_Awal_", label_fcm_awal))
df$max_memb_awal <- max_memb_awal
df$status_memb_awal <- dplyr::case_when(
  max_memb_awal < 0.60 ~ "Ambiguous",
  max_memb_awal < 0.75 ~ "Moderate",
  TRUE ~ "Strong"
)

fcm_awal_external <- data.frame(
  purity = purity_score(df$cluster_fcm_awal, df$disease_bin),
  ari = ari_score(df$cluster_fcm_awal, df$disease_bin),
  nmi = nmi_score(df$cluster_fcm_awal, df$disease_bin)
)

cat("Silhouette FCM awal:", round(sil_fcm_awal_mean, 4), "\n")
cat("Dunn Index FCM awal:", round(dunn_fcm_awal, 4), "\n")
cat("Davies Bouldin Index FCM awal:", round(dbi_fcm_awal, 4), "\n")
cat("Partition Coefficient FCM awal:", round(pc_pe_awal$partition_coefficient, 4), "\n")
cat("Partition Entropy FCM awal:", round(pc_pe_awal$partition_entropy, 4), "\n")
cat("Modified Partition Coefficient FCM awal:", round(pc_pe_awal$modified_partition_coefficient, 4), "\n")
cat("Xie Beni Index FCM awal:", round(xb_fcm_awal, 4), "\n")

cat("\nEvaluasi eksternal FCM awal terhadap label disease\n")
print(round(fcm_awal_external, 4))

cat("\nStatus membership FCM awal\n")
print(table(df$status_memb_awal))

tab_fcm_awal <- table(df$cluster_fcm_awal, df$disease_bin)

cat("\nTabulasi FCM awal terhadap disease\n")
print(tab_fcm_awal)
print(round(prop.table(tab_fcm_awal, margin = 1) * 100, 2))

cat("\nProfil cluster FCM awal\n")
profil_fcm_awal <- profil_cluster(df, "cluster_fcm_awal")
print(profil_fcm_awal)

# Fuzzy C Means model (Optimzed)
cat("\nBagian Fuzzy C Means model optimized\n")

fcm_feature_sets <- list(
  numerik_lengkap = c(
    "age", "trestbps", "chol", "thalch", "oldpeak"
  ),
  numerik_tanpa_chol = c(
    "age", "trestbps", "thalch", "oldpeak"
  ),
  numerik_risiko_inti = c(
    "age", "thalch", "oldpeak"
  )
)

c_values <- 2:6
m_values <- c(1.5, 2.0, 2.5)

eval_fcm <- data.frame()
fcm_objects <- list()

for(fs_name in names(fcm_feature_sets)){
  fitur_pakai <- fcm_feature_sets[[fs_name]]
  
  data_fcm_temp <- df %>%
    dplyr::select(dplyr::all_of(fitur_pakai))
  
  data_fcm_temp_scaled <- scale(data_fcm_temp)
  dist_fcm_temp <- dist(data_fcm_temp_scaled, method = "euclidean")
  
  for(c_val in c_values){
    for(m_val in m_values){
      set.seed(123)
      
      fcm_temp <- e1071::cmeans(
        data_fcm_temp_scaled,
        centers = c_val,
        iter.max = 300,
        m = m_val,
        method = "cmeans",
        dist = "euclidean"
      )
      
      label_temp <- fcm_temp$cluster
      membership_temp <- as.matrix(fcm_temp$membership)
      max_memb_temp <- apply(membership_temp, 1, max)
      
      cluster_size_temp <- table(label_temp)
      min_cluster_size <- min(cluster_size_temp)
      max_cluster_size <- max(cluster_size_temp)
      balance_ratio <- min_cluster_size / max_cluster_size
      
      sil_temp <- safe_silhouette_mean(label_temp, dist_fcm_temp)
      dunn_temp <- dunn_index_from_dist(label_temp, dist_fcm_temp)
      dbi_temp <- db_index_numeric(data_fcm_temp_scaled, label_temp)
      pc_pe_temp <- partition_metrics(membership_temp)
      xb_temp <- xie_beni_index(
        data_fcm_temp_scaled,
        fcm_temp$centers,
        membership_temp,
        m_val
      )
      
      ambiguous_temp <- sum(max_memb_temp < 0.60)
      moderate_temp <- sum(max_memb_temp >= 0.60 & max_memb_temp < 0.75)
      strong_temp <- sum(max_memb_temp >= 0.75)
      
      key_name <- paste(fs_name, c_val, m_val, sep = "_")
      
      eval_fcm <- rbind(
        eval_fcm,
        data.frame(
          key = key_name,
          feature_set = fs_name,
          c = c_val,
          m = m_val,
          silhouette = sil_temp,
          dunn_index = dunn_temp,
          davies_bouldin_index = dbi_temp,
          partition_coefficient = pc_pe_temp$partition_coefficient,
          partition_entropy = pc_pe_temp$partition_entropy,
          modified_partition_coefficient = pc_pe_temp$modified_partition_coefficient,
          xie_beni_index = xb_temp,
          ambiguous = ambiguous_temp,
          moderate = moderate_temp,
          strong = strong_temp,
          balance_ratio = as.numeric(balance_ratio),
          min_cluster_size = as.integer(min_cluster_size),
          cluster_size = paste(
            names(cluster_size_temp),
            as.integer(cluster_size_temp),
            sep = ":",
            collapse = "; "
          )
        )
      )
      
      fcm_objects[[key_name]] <- list(
        fitur = fitur_pakai,
        data_scaled = data_fcm_temp_scaled,
        dist = dist_fcm_temp,
        model = fcm_temp,
        label = label_temp,
        membership = membership_temp,
        max_membership = max_memb_temp
      )
    }
  }
}

eval_fcm <- eval_fcm %>%
  dplyr::mutate(
    dunn_score = rescale01(dunn_index),
    dbi_score = 1 - rescale01(davies_bouldin_index),
    xb_score = 1 - rescale01(xie_beni_index),
    pe_score = 1 - rescale01(partition_entropy),
    ambiguous_rate = ambiguous / nrow(df),
    score = silhouette +
      0.20 * modified_partition_coefficient +
      0.15 * dunn_score +
      0.15 * dbi_score +
      0.15 * xb_score +
      0.10 * pe_score -
      0.20 * ambiguous_rate +
      0.05 * balance_ratio
  ) %>%
  dplyr::arrange(dplyr::desc(score))

cat("\nEvaluasi kombinasi FCM\n")
print(eval_fcm)

eval_fcm_filtered <- eval_fcm %>%
  dplyr::filter(
    min_cluster_size >= min_size_rule,
    !is.na(silhouette)
  ) %>%
  dplyr::arrange(dplyr::desc(score))

if(nrow(eval_fcm_filtered) == 0){
  best_fcm <- eval_fcm[1, ]
} else {
  best_fcm <- eval_fcm_filtered[1, ]
}

cat("\nModel FCM optimal yang dipilih\n")
print(best_fcm)

fcm_opt_obj <- fcm_objects[[best_fcm$key]]

fcm_opt <- fcm_opt_obj$model
label_fcm_opt <- fcm_opt_obj$label
membership_fcm_opt <- fcm_opt_obj$membership
max_memb_opt <- fcm_opt_obj$max_membership
data_fcm_opt_scaled <- fcm_opt_obj$data_scaled
dist_fcm_opt <- fcm_opt_obj$dist

sil_fcm_opt <- cluster::silhouette(label_fcm_opt, dist_fcm_opt)
sil_fcm_opt_mean <- mean(sil_fcm_opt[, 3])
dunn_fcm_opt <- dunn_index_from_dist(label_fcm_opt, dist_fcm_opt)
dbi_fcm_opt <- db_index_numeric(data_fcm_opt_scaled, label_fcm_opt)
pc_pe_opt <- partition_metrics(membership_fcm_opt)
xb_fcm_opt <- xie_beni_index(
  data_fcm_opt_scaled,
  fcm_opt$centers,
  membership_fcm_opt,
  best_fcm$m
)

df$cluster_fcm_opt <- factor(paste0("FCM_Opt_", label_fcm_opt))
df$max_memb_opt <- max_memb_opt
df$status_memb_opt <- dplyr::case_when(
  max_memb_opt < 0.60 ~ "Ambiguous",
  max_memb_opt < 0.75 ~ "Moderate",
  TRUE ~ "Strong"
)

fcm_opt_external <- data.frame(
  purity = purity_score(df$cluster_fcm_opt, df$disease_bin),
  ari = ari_score(df$cluster_fcm_opt, df$disease_bin),
  nmi = nmi_score(df$cluster_fcm_opt, df$disease_bin)
)

cat("\nSilhouette FCM optimal:", round(sil_fcm_opt_mean, 4), "\n")
cat("Dunn Index FCM optimal:", round(dunn_fcm_opt, 4), "\n")
cat("Davies Bouldin Index FCM optimal:", round(dbi_fcm_opt, 4), "\n")
cat("Partition Coefficient FCM optimal:", round(pc_pe_opt$partition_coefficient, 4), "\n")
cat("Partition Entropy FCM optimal:", round(pc_pe_opt$partition_entropy, 4), "\n")
cat("Modified Partition Coefficient FCM optimal:", round(pc_pe_opt$modified_partition_coefficient, 4), "\n")
cat("Xie Beni Index FCM optimal:", round(xb_fcm_opt, 4), "\n")

cat("\nEvaluasi eksternal FCM optimal terhadap label disease\n")
print(round(fcm_opt_external, 4))

cat("\nStatus membership FCM optimal\n")
print(table(df$status_memb_opt))

tab_fcm_opt <- table(df$cluster_fcm_opt, df$disease_bin)

cat("\nTabulasi FCM optimal terhadap disease\n")
print(tab_fcm_opt)
print(round(prop.table(tab_fcm_opt, margin = 1) * 100, 2))

cat("\nProfil cluster FCM optimal\n")
profil_fcm_opt <- profil_cluster(df, "cluster_fcm_opt")
print(profil_fcm_opt)

# Bagian 11
# Ringkasan perbandingan model

cat("\nBagian ringkasan perbandingan model\n")

perbandingan_hc <- data.frame(
  Metode = c("HC Awal", "HC Optimal"),
  Jumlah_Cluster = c(k_awal, best_hc$k),
  Silhouette = c(
    round(sil_hc_awal_mean, 4),
    round(sil_hc_opt_mean, 4)
  ),
  Dunn_Index = c(
    round(dunn_hc_awal, 4),
    round(dunn_hc_opt, 4)
  ),
  Within_Between_Ratio = c(
    round(ratio_hc_awal, 4),
    round(ratio_hc_opt, 4)
  ),
  Cophenetic_Correlation = c(
    round(coph_hc_awal, 4),
    round(coph_hc_opt, 4)
  ),
  Purity = c(
    round(hc_awal_external$purity, 4),
    round(hc_opt_external$purity, 4)
  ),
  ARI = c(
    round(hc_awal_external$ari, 4),
    round(hc_opt_external$ari, 4)
  ),
  NMI = c(
    round(hc_awal_external$nmi, 4),
    round(hc_opt_external$nmi, 4)
  ),
  Keterangan = c(
    "Fitur campuran lengkap dengan Gower dan average linkage",
    paste0(
      best_hc$feature_set,
      " dengan ",
      best_hc$method,
      " dan k ",
      best_hc$k
    )
  )
)

perbandingan_apriori <- data.frame(
  Metode = c("Apriori Awal", "Apriori Optimal"),
  Support = c(0.10, best_apriori$support),
  Confidence = c(0.70, best_apriori$confidence),
  Jumlah_Rules = c(
    length(rules_awal),
    length(rules_opt)
  ),
  Rules_Non_Redundant = c(
    length(rules_awal_nonred),
    length(rules_opt)
  ),
  Rules_Disease = c(
    length(rules_awal_disease),
    length(rules_opt_disease)
  ),
  Rules_No_Disease = c(
    length(rules_awal_no_disease),
    length(rules_opt_no_disease)
  ),
  Rules_Cluster = c(
    NA,
    length(rules_opt_cluster)
  ),
  Rules_Risk_Cluster = c(
    NA,
    length(rules_risk_cluster)
  )
)

perbandingan_fcm <- data.frame(
  Metode = c("FCM Awal", "FCM Optimal"),
  Jumlah_Cluster = c(c_fcm_awal, best_fcm$c),
  m = c(m_fcm_awal, best_fcm$m),
  Silhouette = c(
    round(sil_fcm_awal_mean, 4),
    round(sil_fcm_opt_mean, 4)
  ),
  Dunn_Index = c(
    round(dunn_fcm_awal, 4),
    round(dunn_fcm_opt, 4)
  ),
  Davies_Bouldin_Index = c(
    round(dbi_fcm_awal, 4),
    round(dbi_fcm_opt, 4)
  ),
  Partition_Coefficient = c(
    round(pc_pe_awal$partition_coefficient, 4),
    round(pc_pe_opt$partition_coefficient, 4)
  ),
  Partition_Entropy = c(
    round(pc_pe_awal$partition_entropy, 4),
    round(pc_pe_opt$partition_entropy, 4)
  ),
  Modified_Partition_Coefficient = c(
    round(pc_pe_awal$modified_partition_coefficient, 4),
    round(pc_pe_opt$modified_partition_coefficient, 4)
  ),
  Xie_Beni_Index = c(
    round(xb_fcm_awal, 4),
    round(xb_fcm_opt, 4)
  ),
  Ambiguous = c(
    sum(df$status_memb_awal == "Ambiguous"),
    sum(df$status_memb_opt == "Ambiguous")
  ),
  Purity = c(
    round(fcm_awal_external$purity, 4),
    round(fcm_opt_external$purity, 4)
  ),
  ARI = c(
    round(fcm_awal_external$ari, 4),
    round(fcm_opt_external$ari, 4)
  ),
  NMI = c(
    round(fcm_awal_external$nmi, 4),
    round(fcm_opt_external$nmi, 4)
  ),
  Keterangan = c(
    "Fitur numerik lengkap dengan c 2 dan m 2",
    paste0(
      best_fcm$feature_set,
      " dengan c ",
      best_fcm$c,
      " dan m ",
      best_fcm$m
    )
  )
)

cat("\nPerbandingan Hierarchical Clustering\n")
print(perbandingan_hc)

cat("\nPerbandingan Apriori\n")
print(perbandingan_apriori)

cat("\nPerbandingan Fuzzy C Means\n")
print(perbandingan_fcm)

ringkasan_clustering <- data.frame(
  Metode = c(
    "HC Awal",
    "HC Optimal",
    "FCM Awal",
    "FCM Optimal"
  ),
  Silhouette = c(
    round(sil_hc_awal_mean, 4),
    round(sil_hc_opt_mean, 4),
    round(sil_fcm_awal_mean, 4),
    round(sil_fcm_opt_mean, 4)
  ),
  Purity = c(
    round(hc_awal_external$purity, 4),
    round(hc_opt_external$purity, 4),
    round(fcm_awal_external$purity, 4),
    round(fcm_opt_external$purity, 4)
  ),
  ARI = c(
    round(hc_awal_external$ari, 4),
    round(hc_opt_external$ari, 4),
    round(fcm_awal_external$ari, 4),
    round(fcm_opt_external$ari, 4)
  ),
  NMI = c(
    round(hc_awal_external$nmi, 4),
    round(hc_opt_external$nmi, 4),
    round(fcm_awal_external$nmi, 4),
    round(fcm_opt_external$nmi, 4)
  )
)

cat("\nRingkasan metrik clustering utama\n")
print(ringkasan_clustering)

# Visualisasi hasil akhir
pca_data <- df %>%
  dplyr::select(age, trestbps, chol, thalch, oldpeak)

pca_scaled <- scale(pca_data)
pca_result <- prcomp(pca_scaled)

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Disease = df$disease_bin,
  HC_Opt = df$cluster_hc_opt,
  FCM_Opt = df$cluster_fcm_opt,
  Max_Membership = df$max_memb_opt,
  Membership_Status = df$status_memb_opt
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = HC_Opt)) +
  geom_point(alpha = 0.75) +
  labs(
    title = "Visualisasi HC optimal menggunakan PCA",
    x = "PC1",
    y = "PC2",
    color = "Cluster HC"
  ) +
  theme_minimal()

ggplot(pca_df, aes(x = PC1, y = PC2, color = FCM_Opt)) +
  geom_point(aes(size = Max_Membership), alpha = 0.75) +
  labs(
    title = "Visualisasi FCM optimal menggunakan PCA",
    x = "PC1",
    y = "PC2",
    color = "Cluster FCM",
    size = "Membership maksimum"
  ) +
  theme_minimal()

ggplot(pca_df, aes(x = PC1, y = PC2, color = Disease)) +
  geom_point(alpha = 0.75) +
  labs(
    title = "Sebaran status penyakit pada PCA",
    x = "PC1",
    y = "PC2",
    color = "Status penyakit"
  ) +
  theme_minimal()

ggplot(perbandingan_hc, aes(x = Metode, y = Silhouette, fill = Metode)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = Silhouette), vjust = -0.4, size = 4) +
  labs(
    title = "Perbandingan silhouette Hierarchical Clustering",
    x = "Model",
    y = "Silhouette"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(perbandingan_fcm, aes(x = Metode, y = Silhouette, fill = Metode)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = Silhouette), vjust = -0.4, size = 4) +
  labs(
    title = "Perbandingan silhouette Fuzzy C Means",
    x = "Model",
    y = "Silhouette"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(ringkasan_clustering, aes(x = Metode, y = Silhouette, fill = Metode)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = Silhouette), vjust = -0.4, size = 4) +
  labs(
    title = "Perbandingan silhouette seluruh model clustering",
    x = "Model",
    y = "Silhouette"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

if(length(rules_opt) > 0){
  rules_opt_df <- rules_to_df(rules_opt, trans_opt)
  
  ggplot(rules_opt_df, aes(x = support, y = confidence, size = lift)) +
    geom_point(alpha = 0.7) +
    labs(
      title = "Sebaran rules Apriori optimal",
      x = "Support",
      y = "Confidence",
      size = "Lift"
    ) +
    theme_minimal()
}

df_membership_status <- as.data.frame(table(df$status_memb_opt))
colnames(df_membership_status) <- c("Status", "Jumlah")

ggplot(df_membership_status, aes(x = Status, y = Jumlah, fill = Status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = Jumlah), vjust = -0.4, size = 4) +
  labs(
    title = "Distribusi status membership FCM optimal",
    x = "Status membership",
    y = "Jumlah pasien"
  ) +
  theme_minimal() +
  theme(legend.position = "none")