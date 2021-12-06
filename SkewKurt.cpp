# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
uint64_t choose_rcpp(uint64_t n, uint64_t k){
  if(k == 0) return 1;
  return (n * choose_rcpp(n - 1, k - 1)) / k;
}

// [[Rcpp::export]]
double asymskewness(Rcpp::NumericMatrix Data){
  int n = Data.nrow();
  int p = Data.ncol();
  
  uint64_t i;
  int j, k, l;
  
  uint64_t allsamplenum = (uint64_t)pow(2, p) - 1;
  Rcpp::LogicalMatrix allsampleindices(allsamplenum, p);
  uint64_t row_position = 0;
  for (k = 1; k <= p; ++k){
    std::string bitmask(k, 1);
    bitmask.resize(p, 0);
    do{
      for (j = 0; j < p; ++j){
        if (bitmask[j]){
          allsampleindices(row_position, j) = TRUE;
        }else{
          allsampleindices(row_position, j) = FALSE;
        }
      }
      row_position = row_position + 1;
    }while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  }
  
  arma::vec q_vector(allsamplenum);
  for (i = 0; i < allsamplenum; ++i){
    q_vector[i] = 0;
    for (j = 0; j < p; ++j){
      if (allsampleindices(i, j)){
        q_vector[i] = q_vector[i] + 1;
      }
    }
  }
  
  arma::vec DataMean(p);
  arma::mat DataVar(p, p);
  
  for (k = 0; k < p; ++k){
    DataMean[k] = 0;
    for (j = 0; j < n; ++j){
      DataMean[k] = DataMean[k] + Data(j, k);
    }
    DataMean[k] = DataMean[k] / n;
  }
  for (k = 0; k < p; ++k){
    for (j = 0; j < p; ++j){
      DataVar(k, j) = 0;
      for (l = 0; l < n; ++l){
        DataVar(k, j) = DataVar(k, j) + ((Data(l, k) - DataMean[k]) * (Data(l, j) - DataMean[j]));
      }
      DataVar(k, j) = DataVar(k, j) / (n - 1);
    }
  }
  
  uint64_t colnum_Z = 0;
  for (k = 1; k <= p; ++k){
    if (k == 1){
      colnum_Z = colnum_Z + (uint64_t)p;
    }else{
      colnum_Z = colnum_Z + choose_rcpp(p, k) * ((uint64_t)(k + ((k * (k - 1) * (k + 4)) / 6)));
    }
  }
  
  arma::vec m1(allsamplenum);
  arma::vec m1_n_mean(allsamplenum);
  arma::vec m1_n_sd(allsamplenum);
  arma::mat Z(n, colnum_Z);
  arma::vec groupsizes(allsamplenum);
  Rcpp::LogicalVector currentsampleindices(p);
  int q, currentrow, currentcolumn;
  uint64_t colindex_Z = 0;
  double b1_qi, temp_var, temp_scale;
  for (i = 0; i < allsamplenum; ++i){
    q = q_vector[i];
    for (j = 0; j < p; ++j){
      currentsampleindices[j] = allsampleindices(i, j);
    }
    
    arma::mat Data_i(n, q);
    arma::vec Data_i_mean(q);
    currentcolumn = 0;
    for (j = 0; j < p; ++j){
      if (currentsampleindices[j]){
        Data_i_mean[currentcolumn] = DataMean[j];
        for (k = 0; k < n; ++k){
          Data_i(k, currentcolumn) = Data(k, j);
        }
        currentcolumn = currentcolumn + 1;
      }
    }
    arma::mat Data_i_centered(n, q);
    for (k = 0; k < n; ++k){
      for (j = 0; j < q; ++j){
        Data_i_centered(k, j) = Data_i(k, j) - Data_i_mean[j];
      }
    }
    arma::mat Data_i_var(q, q);
    currentrow = 0;
    for (j = 0; j < p; ++j){
      if (currentsampleindices[j]){
        currentcolumn = 0;
        for (k = 0; k < p; ++k){
          if (currentsampleindices[k]){
            Data_i_var(currentrow, currentcolumn) = DataVar(j, k);
            currentcolumn = currentcolumn + 1;
          }
        }
        currentrow = currentrow + 1;
      }
    }
    
    arma::mat sqrt_Data_i_var(q, q);
    if (q > 1){
      arma::vec eigenvalues_vector;
      arma::mat eigenvectors_matrix;
      eig_sym(eigenvalues_vector, eigenvectors_matrix, Data_i_var);
      
      arma::vec sqrt_eigenvalues_vector(q);
      for (j = 0; j < q; ++j){
        sqrt_eigenvalues_vector[j] = sqrt(abs(eigenvalues_vector[j]));
      }
      for (j = 0; j < q; ++j){
        for (k = 0; k < q; ++k){
          sqrt_Data_i_var(j, k) = 0;
          for (l = 0; l < q; ++l){
            sqrt_Data_i_var(j, k) = sqrt_Data_i_var(j, k) +
              (eigenvectors_matrix(j, l) * sqrt_eigenvalues_vector[l] * eigenvectors_matrix(k, l));
          }
        }
      }
    }else{
      sqrt_Data_i_var(0, 0) = sqrt(Data_i_var(0, 0));
    }
    
    arma::mat Data_i_centeredscaled(n, q);
    arma::mat t_Data_i_centered(q, n), t_Data_i_centeredscaled(q, n);
    for (j = 0; j < q; ++j){
      for (k = 0; k < n; ++k){
        t_Data_i_centered(j, k) = Data_i_centered(k, j);
      }
    }
    t_Data_i_centeredscaled = arma::solve(sqrt_Data_i_var, t_Data_i_centered);
    for (k = 0; k < n; ++k){
      for (j = 0; j < q; ++j){
        Data_i_centeredscaled(k, j) = t_Data_i_centeredscaled(j, k);
      }
    }
    
    arma::mat X_Sinv_Y(n, n);
    for (j = 0; j < n; ++j){
      for (k = 0; k < n; ++k){
        X_Sinv_Y(j, k) = 0;
        for (l = 0; l < q; ++l){
          X_Sinv_Y(j, k) = X_Sinv_Y(j, k) + (Data_i_centeredscaled(j, l) * Data_i_centeredscaled(k, l));
        }
      }
    }
    
    b1_qi = 0;
    for (j = 0; j < n; ++j){
      for (k = 0; k < n; ++k){
        b1_qi = b1_qi + pow(X_Sinv_Y(j, k), 3);
      }
    }
    b1_qi = b1_qi / (n * n);
    
    m1[i] = b1_qi;
    
    m1_n_mean[i] = q * (q + 1) * (q + 2);
    m1_n_sd[i] = sqrt(12 * m1_n_mean[i]);
    
    if (q == 1){
      for (k = 0; k < n; ++k){
        Z(k, colindex_Z) = (Data_i_centered(k, 0) * (pow(Data_i_centered(k, 0), 2) - (3 * Data_i_var(0, 0)))) /
          pow(Data_i_var(0, 0), 1.5);
      }
      temp_var = 0;
      for (k = 0; k < n; ++k){
        temp_var = temp_var + (Z(k, colindex_Z) * Z(k, colindex_Z));
      }
      temp_var = temp_var / (n - 1);
      temp_scale = sqrt(6) / sqrt(temp_var);
      for (k = 0; k < n; ++k){
        Z(k, colindex_Z) = Z(k, colindex_Z) * temp_scale;
      }
      
      colindex_Z = colindex_Z + 1;
      
      groupsizes[i] = 1;
    }else{
      arma::vec normsq(n);
      for (k = 0; k < n; ++k){
        normsq[k] = X_Sinv_Y(k, k);
      }
      
      double alpha_1 = 6;
      double alpha_2 = 6;
      
      int multiplicity_1 = q;
      int multiplicity_2 = (int)(((q * (q - 1) * (q + 4)) / 6) + 0.2);
      int eigen_number = multiplicity_1 + multiplicity_2;
      
      groupsizes[i] = eigen_number;
      
      arma::vec h_eigenvalues_nonzero(eigen_number);
      if (alpha_1 >= alpha_2){
        for (j = 0; j < multiplicity_1; ++j){
          h_eigenvalues_nonzero[j] = alpha_1;
        }
        for (j = multiplicity_1; j < eigen_number; ++j){
          h_eigenvalues_nonzero[j] = alpha_2;
        }
      }else{
        for (j = 0; j < multiplicity_2; ++j){
          h_eigenvalues_nonzero[j] = alpha_2;
        }
        for (j = multiplicity_2; j < eigen_number; ++j){
          h_eigenvalues_nonzero[j] = alpha_1;
        }
      }
      
      Rcpp::NumericMatrix h_matrix(n, n);
      for (j = 0; j < n; ++j){
        for (k = 0; k < n; ++k){
          h_matrix(j, k) = pow(X_Sinv_Y(j, k), 3) - (3 * (normsq[j] + normsq[k]) * X_Sinv_Y(j, k)) +
            (3 * (q + 2) * X_Sinv_Y(j, k));
        }
      }
      
      Rcpp::Environment pkg = Rcpp::Environment::namespace_env("RSpectra");
      Rcpp::Function eigen_RSpectra = pkg["eigs_sym"];
      
      Rcpp::List h_eigen = eigen_RSpectra(h_matrix, Rcpp::Named("k", eigen_number), Rcpp::Named("which", "LM"));
      
      Rcpp::NumericVector h_eigenvalues = h_eigen["values"];
      Rcpp::NumericMatrix h_eigenvectors = h_eigen["vectors"];
      
      for (j = 0; j < eigen_number; ++j){
        for (k = 0; k < n; ++k){
          h_eigenvectors(k, j) = h_eigenvectors(k, j) * sqrt(n);
        }
      }
      
      Rcpp::IntegerVector indices(eigen_number);
      for (int j = 0; j < eigen_number; ++j){
        indices[j] = j;
      }
      std::sort(indices.begin(), indices.end(), [&](int i, int j){return h_eigenvalues[i] > h_eigenvalues[j];});
      
      Rcpp::NumericMatrix h_eigenvectors_sorted(n, eigen_number);
      for (j = 0; j < eigen_number; ++j){
        if (h_eigenvectors(0, indices[j]) >= 0){
          for (k = 0; k < n; ++k){
            h_eigenvectors_sorted(k, j) = h_eigenvectors(k, indices[j]);
          }
        }else{
          for (k = 0; k < n; ++k){
            h_eigenvectors_sorted(k, j) = - h_eigenvectors(k, indices[j]);
          }
        }
      }
      
      arma::mat sqrtLambda_fk_matrix(n, eigen_number);
      for (j = 0; j < eigen_number; ++j){
        for (k = 0; k < n; ++k){
          sqrtLambda_fk_matrix(k, j) = sqrt(h_eigenvalues_nonzero[j]) * h_eigenvectors_sorted(k, j);
        }
      }
      
      for (j = 0; j < eigen_number; ++j){
        for (k = 0; k < n; ++k){
          Z(k, colindex_Z) = sqrtLambda_fk_matrix(k, j);
        }
        colindex_Z = colindex_Z + 1;
      }
    }
  }
  
  uint64_t i1, i2;
  arma::vec mean_Z(colnum_Z);
  for (i1 = 0; i1 < colnum_Z; ++i1){
    mean_Z[i1] = 0;
    for (k = 0; k < n; ++k){
      mean_Z[i1] = mean_Z[i1] + Z(k, i1);
    }
    mean_Z[i1] = mean_Z[i1] / n;
  }
  arma::mat cov_Z(colnum_Z, colnum_Z);
  for (i1 = 0; i1 < colnum_Z; ++i1){
    for (i2 = 0; i2 < colnum_Z; ++i2){
      cov_Z(i1, i2) = 0;
      for (k = 0; k < n; ++k){
        cov_Z(i1, i2) = cov_Z(i1, i2) + ((Z(k, i1) - mean_Z[i1]) * (Z(k, i2) - mean_Z[i2]));
      }
      cov_Z(i1, i2) = cov_Z(i1, i2) / (n - 1);
    }
  }
  
  arma::vec mean_MVN(colnum_Z);
  for (i = 0; i < colnum_Z; ++i){
    mean_MVN[i] = 0;
  }
  int n_bootstrap = 1000;
  arma::mat simulated_temp = arma::mvnrnd(mean_MVN, cov_Z, n_bootstrap);
  
  arma::vec simulatedvalues(n_bootstrap);
  uint64_t current_position, position_end;
  double currentgroupsum, currentvalue;
  for (k = 0; k < n_bootstrap; ++k){
    current_position = 0;
    for (i = 0; i < allsamplenum; ++i){
      position_end = current_position + groupsizes[i];
      currentgroupsum = 0;
      while (current_position < position_end){
        currentgroupsum = currentgroupsum + pow(simulated_temp(current_position, k), 2);
        current_position = current_position + 1;
      }
      currentvalue = (currentgroupsum - m1_n_mean[i]) / m1_n_sd[i];
      if (i == 0){
        simulatedvalues[k] = currentvalue;
      }
      if (simulatedvalues[k] < currentvalue){
        simulatedvalues[k] = currentvalue;
      }
    }
  }
  
  double test_statistic = ((n * m1[0]) - m1_n_mean[0]) / m1_n_sd[0];
  for (i = 1; i < allsamplenum; ++i){
    currentvalue = ((n * m1[i]) - m1_n_mean[i]) / m1_n_sd[i];
    if (test_statistic < currentvalue){
      test_statistic = currentvalue;
    }
  }
  
  double pvalue = 0;
  for (k = 0; k < n_bootstrap; ++k){
    if (simulatedvalues[k] >= test_statistic){
      pvalue = pvalue + 1;
    }
  }
  pvalue = pvalue / n_bootstrap;
  
  return pvalue;
}

// [[Rcpp::export]]
double asymskewness_q(Rcpp::NumericMatrix Data, int q){
  int n = Data.nrow();
  int p = Data.ncol();
  
  uint64_t i;
  int j, k, l;
  
  uint64_t allsamplenum = choose_rcpp(p, q);
  Rcpp::LogicalMatrix allsampleindices(allsamplenum, p);
  uint64_t row_position = 0;
  std::string bitmask(q, 1);
  bitmask.resize(p, 0);
  do{
    for (j = 0; j < p; ++j){
      if (bitmask[j]){
        allsampleindices(row_position, j) = TRUE;
      }else{
        allsampleindices(row_position, j) = FALSE;
      }
    }
    row_position = row_position + 1;
  }while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  
  arma::vec DataMean(p);
  arma::mat DataVar(p, p);
  
  for (k = 0; k < p; ++k){
    DataMean[k] = 0;
    for (j = 0; j < n; ++j){
      DataMean[k] = DataMean[k] + Data(j, k);
    }
    DataMean[k] = DataMean[k] / n;
  }
  for (k = 0; k < p; ++k){
    for (j = 0; j < p; ++j){
      DataVar(k, j) = 0;
      for (l = 0; l < n; ++l){
        DataVar(k, j) = DataVar(k, j) + ((Data(l, k) - DataMean[k]) * (Data(l, j) - DataMean[j]));
      }
      DataVar(k, j) = DataVar(k, j) / (n - 1);
    }
  }
  
  uint64_t colnum_Z;
  if (q == 1){
    colnum_Z = (uint64_t)p;
  }else{
    colnum_Z = choose_rcpp(p, q) * ((uint64_t)(q + ((q * (q - 1) * (q + 4)) / 6)));
  }
  
  arma::vec m1(allsamplenum);
  arma::vec m1_n_mean(allsamplenum);
  arma::vec m1_n_sd(allsamplenum);
  arma::mat Z(n, colnum_Z);
  arma::vec groupsizes(allsamplenum);
  Rcpp::LogicalVector currentsampleindices(p);
  int currentrow, currentcolumn;
  uint64_t colindex_Z = 0;
  double b1_qi, temp_var, temp_scale;
  for (i = 0; i < allsamplenum; ++i){
    for (j = 0; j < p; ++j){
      currentsampleindices[j] = allsampleindices(i, j);
    }
    
    arma::mat Data_i(n, q);
    arma::vec Data_i_mean(q);
    currentcolumn = 0;
    for (j = 0; j < p; ++j){
      if (currentsampleindices[j]){
        Data_i_mean[currentcolumn] = DataMean[j];
        for (k = 0; k < n; ++k){
          Data_i(k, currentcolumn) = Data(k, j);
        }
        currentcolumn = currentcolumn + 1;
      }
    }
    arma::mat Data_i_centered(n, q);
    for (k = 0; k < n; ++k){
      for (j = 0; j < q; ++j){
        Data_i_centered(k, j) = Data_i(k, j) - Data_i_mean[j];
      }
    }
    arma::mat Data_i_var(q, q);
    currentrow = 0;
    for (j = 0; j < p; ++j){
      if (currentsampleindices[j]){
        currentcolumn = 0;
        for (k = 0; k < p; ++k){
          if (currentsampleindices[k]){
            Data_i_var(currentrow, currentcolumn) = DataVar(j, k);
            currentcolumn = currentcolumn + 1;
          }
        }
        currentrow = currentrow + 1;
      }
    }
    
    arma::mat sqrt_Data_i_var(q, q);
    if (q > 1){
      arma::vec eigenvalues_vector;
      arma::mat eigenvectors_matrix;
      eig_sym(eigenvalues_vector, eigenvectors_matrix, Data_i_var);
      
      arma::vec sqrt_eigenvalues_vector(q);
      for (j = 0; j < q; ++j){
        sqrt_eigenvalues_vector[j] = sqrt(abs(eigenvalues_vector[j]));
      }
      for (j = 0; j < q; ++j){
        for (k = 0; k < q; ++k){
          sqrt_Data_i_var(j, k) = 0;
          for (l = 0; l < q; ++l){
            sqrt_Data_i_var(j, k) = sqrt_Data_i_var(j, k) +
              (eigenvectors_matrix(j, l) * sqrt_eigenvalues_vector[l] * eigenvectors_matrix(k, l));
          }
        }
      }
    }else{
      sqrt_Data_i_var(0, 0) = sqrt(Data_i_var(0, 0));
    }
    
    arma::mat Data_i_centeredscaled(n, q);
    arma::mat t_Data_i_centered(q, n), t_Data_i_centeredscaled(q, n);
    for (j = 0; j < q; ++j){
      for (k = 0; k < n; ++k){
        t_Data_i_centered(j, k) = Data_i_centered(k, j);
      }
    }
    t_Data_i_centeredscaled = arma::solve(sqrt_Data_i_var, t_Data_i_centered);
    for (k = 0; k < n; ++k){
      for (j = 0; j < q; ++j){
        Data_i_centeredscaled(k, j) = t_Data_i_centeredscaled(j, k);
      }
    }
    
    arma::mat X_Sinv_Y(n, n);
    for (j = 0; j < n; ++j){
      for (k = 0; k < n; ++k){
        X_Sinv_Y(j, k) = 0;
        for (l = 0; l < q; ++l){
          X_Sinv_Y(j, k) = X_Sinv_Y(j, k) + (Data_i_centeredscaled(j, l) * Data_i_centeredscaled(k, l));
        }
      }
    }
    
    b1_qi = 0;
    for (j = 0; j < n; ++j){
      for (k = 0; k < n; ++k){
        b1_qi = b1_qi + pow(X_Sinv_Y(j, k), 3);
      }
    }
    b1_qi = b1_qi / (n * n);
    
    m1[i] = b1_qi;
    
    m1_n_mean[i] = q * (q + 1) * (q + 2);
    m1_n_sd[i] = sqrt(12 * m1_n_mean[i]);
    
    if (q == 1){
      for (k = 0; k < n; ++k){
        Z(k, colindex_Z) = (Data_i_centered(k, 0) * (pow(Data_i_centered(k, 0), 2) - (3 * Data_i_var(0, 0)))) /
          pow(Data_i_var(0, 0), 1.5);
      }
      temp_var = 0;
      for (k = 0; k < n; ++k){
        temp_var = temp_var + (Z(k, colindex_Z) * Z(k, colindex_Z));
      }
      temp_var = temp_var / (n - 1);
      temp_scale = sqrt(6) / sqrt(temp_var);
      for (k = 0; k < n; ++k){
        Z(k, colindex_Z) = Z(k, colindex_Z) * temp_scale;
      }
      
      colindex_Z = colindex_Z + 1;
      
      groupsizes[i] = 1;
    }else{
      arma::vec normsq(n);
      for (k = 0; k < n; ++k){
        normsq[k] = X_Sinv_Y(k, k);
      }
      
      double alpha_1 = 6;
      double alpha_2 = 6;
      
      int multiplicity_1 = q;
      int multiplicity_2 = (int)(((q * (q - 1) * (q + 4)) / 6) + 0.2);
      int eigen_number = multiplicity_1 + multiplicity_2;
      
      groupsizes[i] = eigen_number;
      
      arma::vec h_eigenvalues_nonzero(eigen_number);
      if (alpha_1 >= alpha_2){
        for (j = 0; j < multiplicity_1; ++j){
          h_eigenvalues_nonzero[j] = alpha_1;
        }
        for (j = multiplicity_1; j < eigen_number; ++j){
          h_eigenvalues_nonzero[j] = alpha_2;
        }
      }else{
        for (j = 0; j < multiplicity_2; ++j){
          h_eigenvalues_nonzero[j] = alpha_2;
        }
        for (j = multiplicity_2; j < eigen_number; ++j){
          h_eigenvalues_nonzero[j] = alpha_1;
        }
      }
      
      Rcpp::NumericMatrix h_matrix(n, n);
      for (j = 0; j < n; ++j){
        for (k = 0; k < n; ++k){
          h_matrix(j, k) = pow(X_Sinv_Y(j, k), 3) - (3 * (normsq[j] + normsq[k]) * X_Sinv_Y(j, k)) +
            (3 * (q + 2) * X_Sinv_Y(j, k));
        }
      }
      
      Rcpp::Environment pkg = Rcpp::Environment::namespace_env("RSpectra");
      Rcpp::Function eigen_RSpectra = pkg["eigs_sym"];
      
      Rcpp::List h_eigen = eigen_RSpectra(h_matrix, Rcpp::Named("k", eigen_number), Rcpp::Named("which", "LM"));
      
      Rcpp::NumericVector h_eigenvalues = h_eigen["values"];
      Rcpp::NumericMatrix h_eigenvectors = h_eigen["vectors"];
      
      for (j = 0; j < eigen_number; ++j){
        for (k = 0; k < n; ++k){
          h_eigenvectors(k, j) = h_eigenvectors(k, j) * sqrt(n);
        }
      }
      
      Rcpp::IntegerVector indices(eigen_number);
      for (int j = 0; j < eigen_number; ++j){
        indices[j] = j;
      }
      std::sort(indices.begin(), indices.end(), [&](int i, int j){return h_eigenvalues[i] > h_eigenvalues[j];});
      
      Rcpp::NumericMatrix h_eigenvectors_sorted(n, eigen_number);
      for (j = 0; j < eigen_number; ++j){
        if (h_eigenvectors(0, indices[j]) >= 0){
          for (k = 0; k < n; ++k){
            h_eigenvectors_sorted(k, j) = h_eigenvectors(k, indices[j]);
          }
        }else{
          for (k = 0; k < n; ++k){
            h_eigenvectors_sorted(k, j) = - h_eigenvectors(k, indices[j]);
          }
        }
      }
      
      arma::mat sqrtLambda_fk_matrix(n, eigen_number);
      for (j = 0; j < eigen_number; ++j){
        for (k = 0; k < n; ++k){
          sqrtLambda_fk_matrix(k, j) = sqrt(h_eigenvalues_nonzero[j]) * h_eigenvectors_sorted(k, j);
        }
      }
      
      for (j = 0; j < eigen_number; ++j){
        for (k = 0; k < n; ++k){
          Z(k, colindex_Z) = sqrtLambda_fk_matrix(k, j);
        }
        colindex_Z = colindex_Z + 1;
      }
    }
  }
  
  uint64_t i1, i2;
  arma::vec mean_Z(colnum_Z);
  for (i1 = 0; i1 < colnum_Z; ++i1){
    mean_Z[i1] = 0;
    for (k = 0; k < n; ++k){
      mean_Z[i1] = mean_Z[i1] + Z(k, i1);
    }
    mean_Z[i1] = mean_Z[i1] / n;
  }
  arma::mat cov_Z(colnum_Z, colnum_Z);
  for (i1 = 0; i1 < colnum_Z; ++i1){
    for (i2 = 0; i2 < colnum_Z; ++i2){
      cov_Z(i1, i2) = 0;
      for (k = 0; k < n; ++k){
        cov_Z(i1, i2) = cov_Z(i1, i2) + ((Z(k, i1) - mean_Z[i1]) * (Z(k, i2) - mean_Z[i2]));
      }
      cov_Z(i1, i2) = cov_Z(i1, i2) / (n - 1);
    }
  }
  
  arma::vec mean_MVN(colnum_Z);
  for (i = 0; i < colnum_Z; ++i){
    mean_MVN[i] = 0;
  }
  int n_bootstrap = 1000;
  arma::mat simulated_temp = arma::mvnrnd(mean_MVN, cov_Z, n_bootstrap);
  
  arma::vec simulatedvalues(n_bootstrap);
  uint64_t current_position, position_end;
  double currentgroupsum, currentvalue;
  for (k = 0; k < n_bootstrap; ++k){
    current_position = 0;
    for (i = 0; i < allsamplenum; ++i){
      position_end = current_position + groupsizes[i];
      currentgroupsum = 0;
      while (current_position < position_end){
        currentgroupsum = currentgroupsum + pow(simulated_temp(current_position, k), 2);
        current_position = current_position + 1;
      }
      currentvalue = (currentgroupsum - m1_n_mean[i]) / m1_n_sd[i];
      if (i == 0){
        simulatedvalues[k] = currentvalue;
      }
      if (simulatedvalues[k] < currentvalue){
        simulatedvalues[k] = currentvalue;
      }
    }
  }
  
  double test_statistic = ((n * m1[0]) - m1_n_mean[0]) / m1_n_sd[0];
  for (i = 1; i < allsamplenum; ++i){
    currentvalue = ((n * m1[i]) - m1_n_mean[i]) / m1_n_sd[i];
    if (test_statistic < currentvalue){
      test_statistic = currentvalue;
    }
  }
  
  double pvalue = 0;
  for (k = 0; k < n_bootstrap; ++k){
    if (simulatedvalues[k] >= test_statistic){
      pvalue = pvalue + 1;
    }
  }
  pvalue = pvalue / n_bootstrap;
  
  return pvalue;
}

// [[Rcpp::export]]
double asymkurtosis(Rcpp::NumericMatrix Data){
  int n = Data.nrow();
  int p = Data.ncol();
  
  int j, k, l;
  uint64_t i, i1, i2;
  
  uint64_t allsamplenum = (uint64_t)pow(2, p) - 1;
  Rcpp::LogicalMatrix allsampleindices(allsamplenum, p);
  uint64_t row_position = 0;
  for (k = 1; k <= p; ++k){
    std::string bitmask(k, 1);
    bitmask.resize(p, 0);
    do{
      for (j = 0; j < p; ++j){
        if (bitmask[j]){
          allsampleindices(row_position, j) = TRUE;
        }else{
          allsampleindices(row_position, j) = FALSE;
        }
      }
      row_position = row_position + 1;
    }while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  }
  
  arma::vec q_vector(allsamplenum);
  for (i = 0; i < allsamplenum; ++i){
    q_vector[i] = 0;
    for (j = 0; j < p; ++j){
      if (allsampleindices(i, j)){
        q_vector[i] = q_vector[i] + 1;
      }
    }
  }
  
  arma::vec m2(allsamplenum);
  arma::mat Data_W(n, allsamplenum);
  Rcpp::LogicalVector currentsampleindices(p);
  int q, currentcolumn;
  double b2_qi;
  for (i = 0; i < allsamplenum; ++i){
    q = q_vector[i];
    for (j = 0; j < p; ++j){
      currentsampleindices[j] = allsampleindices(i, j);
    }
    
    arma::mat Data_i(n, q);
    arma::vec Data_i_mean(q);
    currentcolumn = 0;
    for (j = 0; j < p; ++j){
      if (currentsampleindices[j]){
        Data_i_mean[currentcolumn] = 0;
        for (k = 0; k < n; ++k){
          Data_i(k, currentcolumn) = Data(k, j);
          Data_i_mean[currentcolumn] = Data_i_mean[currentcolumn] + Data_i(k, currentcolumn);
        }
        Data_i_mean[currentcolumn] = Data_i_mean[currentcolumn] / n;
        currentcolumn = currentcolumn + 1;
      }
    }
    arma::mat Data_i_centered(n, q);
    for (k = 0; k < n; ++k){
      for (j = 0; j < q; ++j){
        Data_i_centered(k, j) = Data_i(k, j) - Data_i_mean[j];
      }
    }
    arma::mat S_i(q, q);
    for (k = 0; k < q; ++k){
      for (j = 0; j < q; ++j){
        S_i(k, j) = 0;
        for (l = 0; l < n; ++l){
          S_i(k, j) = S_i(k, j) + (Data_i_centered(l, k) * Data_i_centered(l, j));
        }
        S_i(k, j) = S_i(k, j) / (n - 1);
      }
    }
    
    arma::mat t_Data_i_centered(q, n);
    for (j = 0; j < q; ++j){
      for (k = 0; k < n; ++k){
        t_Data_i_centered(j, k) = Data_i_centered(k, j);
      }
    }
    
    arma::mat S_i_inv_Data_i_centered = arma::solve(S_i, t_Data_i_centered);
    
    arma::vec vector_i(n);
    for (k = 0; k < n; ++k){
      vector_i[k] = 0;
      for (j = 0; j < q; ++j){
        vector_i[k] = vector_i[k] + (Data_i_centered(k, j) * S_i_inv_Data_i_centered(j, k));
      }
    }
    
    b2_qi = 0;
    for (k = 0; k < n; ++k){
      b2_qi = b2_qi + pow(vector_i[k], 2);
    }
    b2_qi = b2_qi / n;
    
    m2[i] = b2_qi;
    
    for (k = 0; k < n; ++k){
      Data_W(k, i) = pow(vector_i[k], 2) - (2 * (q + 2) * vector_i[k]);
    }
  }
  
  arma::mat Omega_hat(allsamplenum, allsamplenum);
  arma::vec current_vector1(n), current_vector2(n);
  double currentmean1, currentmean2, currentvar1, currentvar2;
  double currentcov, currentcorr;
  for (i1 = 0; i1 < allsamplenum; ++i1){
    for (i2 = 0; i2 < allsamplenum; ++i2){
      if (i1 == i2){
        Omega_hat(i1, i2) = 1;
      }else if (i1 < i2){
        currentmean1 = 0;
        currentmean2 = 0;
        for (k = 0; k < n; ++k){
          current_vector1[k] = Data_W(k, i1);
          current_vector2[k] = Data_W(k, i2);
          
          currentmean1 = currentmean1 + current_vector1[k];
          currentmean2 = currentmean2 + current_vector2[k];
        }
        currentmean1 = currentmean1 / n;
        currentmean2 = currentmean2 / n;
        
        currentvar1 = 0;
        currentvar2 = 0;
        currentcov = 0;
        for (k = 0; k < n; ++k){
          currentvar1 = currentvar1 + pow(current_vector1[k] - currentmean1, 2);
          currentvar2 = currentvar2 + pow(current_vector2[k] - currentmean2, 2);
          currentcov = currentcov + ((current_vector1[k] - currentmean1) * (current_vector2[k] - currentmean2));
        }
        currentvar1 = currentvar1 / (n - 1);
        currentvar2 = currentvar2 / (n - 1);
        currentcov = currentcov / (n - 1);
        currentcorr = currentcov / (sqrtl(currentvar1) * sqrtl(currentvar2));
        
        Omega_hat(i1, i2) = currentcorr;
      }else{
        Omega_hat(i1, i2) = Omega_hat(i2, i1);
      }
    }
  }
  
  arma::vec mean_MVN(allsamplenum);
  for (i = 0; i < allsamplenum; ++i){
    mean_MVN[i] = 0;
  }
  int n_bootstrap = 1000;
  arma::mat simulated_temp = arma::mvnrnd(mean_MVN, Omega_hat, n_bootstrap);
  arma::vec simulatedvalues(n_bootstrap);
  for (k = 0; k < n_bootstrap; ++k){
    simulatedvalues[k] = abs(simulated_temp(0, k));
    for (i = 1; i < allsamplenum; ++i){
      if (simulatedvalues[k] < abs(simulated_temp(i, k))){
        simulatedvalues[k] = abs(simulated_temp(i, k));
      }
    }
  }
  
  int q_current;
  double value_current;
  q_current = q_vector[0];
  double test_statistic = abs(m2[0] - (q_current * (q_current + 2))) / sqrtl(8 * q_current * (q_current + 2));
  for (i = 1; i < allsamplenum; ++i){
    q_current = q_vector[i];
    value_current = abs(m2[i] - (q_current * (q_current + 2))) / sqrtl(8 * q_current * (q_current + 2));
    if (test_statistic < value_current){
      test_statistic = value_current;
    }
  }
  test_statistic = sqrt(n) * test_statistic;
  
  double pvalue = 0;
  for (k = 0; k < n_bootstrap; ++k){
    if (simulatedvalues[k] >= test_statistic){
      pvalue = pvalue + 1;
    }
  }
  pvalue = pvalue / n_bootstrap;
  
  return pvalue;
}

// [[Rcpp::export]]
double asymkurtosis_q(Rcpp::NumericMatrix Data, int q){
  int n = Data.nrow();
  int p = Data.ncol();
  
  int j, k, l;
  uint64_t i, i1, i2;
  
  uint64_t allsamplenum = choose_rcpp(p, q);
  Rcpp::LogicalMatrix allsampleindices(allsamplenum, p);
  uint64_t row_position = 0;
  std::string bitmask(q, 1);
  bitmask.resize(p, 0);
  do{
    for (j = 0; j < p; ++j){
      if (bitmask[j]){
        allsampleindices(row_position, j) = TRUE;
      }else{
        allsampleindices(row_position, j) = FALSE;
      }
    }
    row_position = row_position + 1;
  }while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  
  arma::vec m2(allsamplenum);
  arma::mat Data_W(n, allsamplenum);
  Rcpp::LogicalVector currentsampleindices(p);
  int currentcolumn;
  double b2_qi;
  for (i = 0; i < allsamplenum; ++i){
    for (j = 0; j < p; ++j){
      currentsampleindices[j] = allsampleindices(i, j);
    }
    
    arma::mat Data_i(n, q);
    arma::vec Data_i_mean(q);
    currentcolumn = 0;
    for (j = 0; j < p; ++j){
      if (currentsampleindices[j]){
        Data_i_mean[currentcolumn] = 0;
        for (k = 0; k < n; ++k){
          Data_i(k, currentcolumn) = Data(k, j);
          Data_i_mean[currentcolumn] = Data_i_mean[currentcolumn] + Data_i(k, currentcolumn);
        }
        Data_i_mean[currentcolumn] = Data_i_mean[currentcolumn] / n;
        currentcolumn = currentcolumn + 1;
      }
    }
    arma::mat Data_i_centered(n, q);
    for (k = 0; k < n; ++k){
      for (j = 0; j < q; ++j){
        Data_i_centered(k, j) = Data_i(k, j) - Data_i_mean[j];
      }
    }
    arma::mat S_i(q, q);
    for (k = 0; k < q; ++k){
      for (j = 0; j < q; ++j){
        S_i(k, j) = 0;
        for (l = 0; l < n; ++l){
          S_i(k, j) = S_i(k, j) + (Data_i_centered(l, k) * Data_i_centered(l, j));
        }
        S_i(k, j) = S_i(k, j) / (n - 1);
      }
    }
    
    arma::mat t_Data_i_centered(q, n);
    for (j = 0; j < q; ++j){
      for (k = 0; k < n; ++k){
        t_Data_i_centered(j, k) = Data_i_centered(k, j);
      }
    }
    
    arma::mat S_i_inv_Data_i_centered = arma::solve(S_i, t_Data_i_centered);
    
    arma::vec vector_i(n);
    for (k = 0; k < n; ++k){
      vector_i[k] = 0;
      for (j = 0; j < q; ++j){
        vector_i[k] = vector_i[k] + (Data_i_centered(k, j) * S_i_inv_Data_i_centered(j, k));
      }
    }
    
    b2_qi = 0;
    for (k = 0; k < n; ++k){
      b2_qi = b2_qi + pow(vector_i[k], 2);
    }
    b2_qi = b2_qi / n;
    
    m2[i] = b2_qi;
    
    for (k = 0; k < n; ++k){
      Data_W(k, i) = pow(vector_i[k], 2) - (2 * (q + 2) * vector_i[k]);
    }
  }
  
  arma::mat Omega_hat(allsamplenum, allsamplenum);
  arma::vec current_vector1(n), current_vector2(n);
  double currentmean1, currentmean2, currentvar1, currentvar2;
  double currentcov, currentcorr;
  for (i1 = 0; i1 < allsamplenum; ++i1){
    for (i2 = 0; i2 < allsamplenum; ++i2){
      if (i1 == i2){
        Omega_hat(i1, i2) = 1;
      }else if (i1 < i2){
        currentmean1 = 0;
        currentmean2 = 0;
        for (k = 0; k < n; ++k){
          current_vector1[k] = Data_W(k, i1);
          current_vector2[k] = Data_W(k, i2);
          
          currentmean1 = currentmean1 + current_vector1[k];
          currentmean2 = currentmean2 + current_vector2[k];
        }
        currentmean1 = currentmean1 / n;
        currentmean2 = currentmean2 / n;
        
        currentvar1 = 0;
        currentvar2 = 0;
        currentcov = 0;
        for (k = 0; k < n; ++k){
          currentvar1 = currentvar1 + pow(current_vector1[k] - currentmean1, 2);
          currentvar2 = currentvar2 + pow(current_vector2[k] - currentmean2, 2);
          currentcov = currentcov + ((current_vector1[k] - currentmean1) * (current_vector2[k] - currentmean2));
        }
        currentvar1 = currentvar1 / (n - 1);
        currentvar2 = currentvar2 / (n - 1);
        currentcov = currentcov / (n - 1);
        currentcorr = currentcov / (sqrtl(currentvar1) * sqrtl(currentvar2));
        
        Omega_hat(i1, i2) = currentcorr;
      }else{
        Omega_hat(i1, i2) = Omega_hat(i2, i1);
      }
    }
  }
  
  arma::vec mean_MVN(allsamplenum);
  for (i = 0; i < allsamplenum; ++i){
    mean_MVN[i] = 0;
  }
  int n_bootstrap = 1000;
  arma::mat simulated_temp = arma::mvnrnd(mean_MVN, Omega_hat, n_bootstrap);
  arma::vec simulatedvalues(n_bootstrap);
  for (k = 0; k < n_bootstrap; ++k){
    simulatedvalues[k] = abs(simulated_temp(0, k));
    for (i = 1; i < allsamplenum; ++i){
      if (simulatedvalues[k] < abs(simulated_temp(i, k))){
        simulatedvalues[k] = abs(simulated_temp(i, k));
      }
    }
  }
  
  double value_current;
  double test_statistic = abs(m2[0] - (q * (q + 2))) / sqrtl(8 * q * (q + 2));
  for (i = 1; i < allsamplenum; ++i){
    value_current = abs(m2[i] - (q * (q + 2))) / sqrtl(8 * q * (q + 2));
    if (test_statistic < value_current){
      test_statistic = value_current;
    }
  }
  test_statistic = sqrt(n) * test_statistic;
  
  double pvalue = 0;
  for (k = 0; k < n_bootstrap; ++k){
    if (simulatedvalues[k] >= test_statistic){
      pvalue = pvalue + 1;
    }
  }
  pvalue = pvalue / n_bootstrap;
  
  return pvalue;
}
