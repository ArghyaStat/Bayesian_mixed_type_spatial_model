get_score_cnt = function(prediction_cnt, obs, u_cnt, weights_cnt){
  distr_obs = c()
  for(k in 1:length(u_cnt)){
    distr_obs = cbind(distr_obs, ifelse(u_cnt[k] < obs, 0, 1))
  }
  weights_mat = matrix(weights_cnt, ncol = length(weights_cnt), nrow = length(obs), byrow = TRUE)
  score_cnt = sum(weights_mat * (distr_obs - prediction_cnt)^2)
  score_cnt
}

get_score_ba = function(prediction_ba, obs, u_ba, weights_ba){
  distr_obs = c()
  for(k in 1:length(u_ba)){
    distr_obs = cbind(distr_obs, ifelse(u_ba[k] < obs, 0, 1))
  }
  weights_mat = matrix(weights_ba, ncol = length(weights_ba), nrow = length(obs), byrow = TRUE)
  score_ba = sum(weights_mat * (distr_obs - prediction_ba)^2)
  score_ba
} 

