# ==============================================================================
# Projeto: Regressão Birnbaum-Saunders de 3 Parâmetros (BS3)
# Passo 1: Gerador e Estrutura do Modelo de Regressão
# ==============================================================================

#' Gerador de Variáveis Aleatórias BS de 3 Parâmetros (alpha, beta, gamma)
rbs3 <- function(n, alpha, beta, gamma) {
  if (alpha <= 0 || any(beta <= 0)) stop("alpha e beta devem ser positivos.")
  
  z <- rnorm(n)
  # Transformação exata para BS3
  y <- gamma + beta * (1 + (alpha^2 * z^2)/2 + alpha * z * sqrt(1 + (alpha^2 * z^2)/4))
  return(y)
}

#' Densidade da BS de 3 Parâmetros
dbs3 <- function(y, alpha, beta, gamma, log = FALSE) {
  pdf_vals <- numeric(length(y))
  pos <- y > gamma
  
  if (any(pos)) {
    y_pos <- y[pos]
    termo_x <- y_pos - gamma
    term1 <- (sqrt(termo_x / beta) + sqrt(beta / termo_x)) / (2 * alpha * termo_x)
    z <- (1 / alpha) * (sqrt(termo_x / beta) - sqrt(beta / termo_x))
    pdf_vals[pos] <- term1 * dnorm(z)
  }
  
  if (log) return(log(pdf_vals))
  return(pdf_vals)
}

# ==============================================================================
# Exemplo de Simulação de Dados de Regressão sob H0
# ==============================================================================

set.seed(123)

n <- 40 # Tamanho amostral

# Parâmetros Reais
alpha_real <- 0.5
gamma_real <- 1.5 # Parâmetro de localização
theta0_real <- 0.8 # Intercepto

# Covariável
x1 <- runif(n, min = 0, max = 2)

# Teste de Hipóteses: H0: theta1 = 0 vs H1: theta1 != 0
theta1_h0 <- 0 # Sob H0, o efeito da covariável é nulo

# Vetor de betas dependente do modelo linear log(beta_i) = theta0 + theta1 * x1_i
log_beta <- theta0_real + theta1_h0 * x1
beta_vetor <- exp(log_beta)

# Gerar variável resposta Y
y_simulado <- rbs3(n = n, alpha = alpha_real, beta = beta_vetor, gamma = gamma_real)

# Visualização inicial
plot(x1, y_simulado, pch = 19, col = "darkblue",
     main = "Dados Simulados sob H0: theta1 = 0",
     xlab = "Covariável (x1)", ylab = "Resposta (Y)")


# ==============================================================================
# Projeto: Regressão Birnbaum-Saunders de 3 Parâmetros (BS3)
# Passo 2: Função de Log-Verossimilhança e Teste de Wald para H0: theta1 = 0
# ==============================================================================

# 1. Função de Log-Verossimilhança Negativa (para minimização via optim)
log_lik_bs3_reg <- function(par, y, x) {
  alpha  <- par[1]
  gamma  <- par[2]
  theta0 <- par[3]
  theta1 <- par[4]
  
  # Restrições do domínio
  if (alpha <= 0 || any(y <= gamma)) return(1e10)
  
  # Estrutura de regressão para beta_i
  beta <- exp(theta0 + theta1 * x)
  
  termo_x <- y - gamma
  z <- (1 / alpha) * (sqrt(termo_x / beta) - sqrt(beta / termo_x))
  
  # Densidade BS3 em escala logarítmica
  log_pdf <- log(sqrt(termo_x / beta) + sqrt(beta / termo_x)) - 
             log(2 * alpha * termo_x) + dnorm(z, log = TRUE)
  
  # Retorna a log-verossimilhança negativa
  return(-sum(log_pdf))
}

# 2. Função para Ajustar o Modelo e Calcular a Estatística de Wald
teste_wald_bs3 <- function(y, x, chute_inicial = NULL) {
  if (is.null(chute_inicial)) {
    # Chutes iniciais razoáveis
    gamma_init <- min(y) * 0.9
    alpha_init <- 0.5
    theta0_init <- log(mean(y - gamma_init))
    theta1_init <- 0
    chute_inicial <- c(alpha_init, gamma_init, theta0_init, theta1_init)
  }
  
  # Optimização numérica da verossimilhança
  fit <- optim(
    par = chute_inicial, 
    fn = log_lik_bs3_reg, 
    y = y, 
    x = x, 
    method = "L-BFGS-B",
    lower = c(0.001, -Inf, -Inf, -Inf),
    upper = c(Inf, min(y) - 1e-4, Inf, Inf),
    hessian = TRUE
  )
  
  # Inversão da Hessiana para obter a matriz de variância-covariância
  var_cov <- tryCatch(solve(fit$hessian), error = function(e) return(NULL))
  
  if (is.null(var_cov)) {
    return(c(stat_wald = NA, p_valor = NA, erro = TRUE))
  }
  
  # Estimativa de theta1 e seu Erro Padrão
  hat_theta1 <- fit$par[4]
  se_theta1  <- sqrt(diag(var_cov)[4])
  
  # Estatística de Wald (Normal Padrão)
  w_stat  <- hat_theta1 / se_theta1
  p_valor <- 2 * (1 - pnorm(abs(w_stat))) # Teste bicaudal
  
  return(c(
    hat_alpha  = fit$par[1],
    hat_gamma  = fit$par[2],
    hat_theta0 = fit$par[3],
    hat_theta1 = hat_theta1,
    se_theta1  = se_theta1,
    stat_wald  = w_stat,
    p_valor    = p_valor,
    erro       = FALSE
  ))
}
set.seed(42)

n <- 40
alpha_real  <- 0.5
gamma_real  <- 1.0
theta0_real <- 0.8
theta1_real <- 0 # Sob H0

x <- runif(n, 0, 2)
beta_vetor <- exp(theta0_real + theta1_real * x)

# Gerar Y ~ BS3
y <- rbs3(n, alpha = alpha_real, beta = beta_vetor, gamma = gamma_real)

# Executar o Teste de Wald
resultado_teste <- teste_wald_bs3(y, x)
print(resultado_teste)


# ==============================================================================
# Projeto: Regressão Birnbaum-Saunders de 3 Parâmetros (BS3)
# Passo 3: Simulação de Monte Carlo (25 Cenários: N x n)
# ==============================================================================

library(parallel) # Para execução em paralelo se necessário

set.seed(2024)

# Vetores da grade de simulação
vetor_N <- c(1000, 5000, 10000)
vetor_n <- c(10, 20, 30, 40)

# Parâmetros reais do modelo sob H0
alpha_real  <- 0.5
gamma_real  <- 1.0
theta0_real <- 0.8
theta1_real <- 0.0 # H0 verdadeira: sem efeito da covariável
alpha_nominal <- 0.05

# Função para executar 1 cenário (N réplicas com amostra n)
simular_cenario_wald <- function(N, n) {
  p_valores <- numeric(N)
  erros_convergencia <- 0
  
  for (i in 1:N) {
    # 1. Gerar covariável x
    x <- runif(n, min = 0, max = 2)
    
    # 2. Gerar resposta Y sob H0 (beta_i = exp(theta0))
    beta_vetor <- exp(theta0_real + theta1_real * x)
    y <- rbs3(n = n, alpha = alpha_real, beta = beta_vetor, gamma = gamma_real)
    
    # 3. Executar o Teste de Wald
    res <- teste_wald_bs3(y, x)
    
    # Controlar possíveis falhas na inversão da matriz Hessiana
    if (is.na(res["p_valor"]) || res["erro"] == TRUE) {
      erros_convergencia <- erros_convergencia + 1
      p_valores[i] <- NA
    } else {
      p_valores[i] <- res["p_valor"]
    }
  }
  
  # Remover réplicas com falha numérica
  p_validos <- na.omit(p_valores)
  N_efetivo <- length(p_validos)
  
  # Métricas do Erro do Tipo I
  taxa_erro <- mean(p_validos <= alpha_nominal)
  erro_padrao <- sqrt(taxa_erro * (1 - taxa_erro) / N_efetivo)
  
  return(data.frame(
    N = N,
    n = n,
    N_efetivo = N_efetivo,
    falhas = erros_convergencia,
    taxa_erro_tipo1 = round(taxa_erro, 5),
    erro_padrao = round(erro_padrao, 5),
    ic_inferior = round(taxa_erro - 1.96 * erro_padrao, 5),
    ic_superior = round(taxa_erro + 1.96 * erro_padrao, 5)
  ))
}

# Criar a grade dos 25 cenários
grade_25 <- expand.grid(N = vetor_N, n = vetor_n)

cat("Iniciando a simulação de Monte Carlo para os 25 cenários...\n")
t_inicio <- Sys.time()

resultados_wald_lista <- vector("list", nrow(grade_25))

for (k in 1:nrow(grade_25)) {
  cat(sprintf("[%d/25] Executando N = %d, n = %d...\n", k, grade_25$N[k], grade_25$n[k]))
  resultados_wald_lista[[k]] <- simular_cenario_wald(
    N = grade_25$N[k], 
    n = grade_25$n[k]
  )
}

t_fim <- Sys.time()
df_resultados_wald <- do.call(rbind, resultados_wald_lista)

cat("\nSimulação concluída em:", round(difftime(t_fim, t_inicio, units = "mins"), 2), "minutos.\n")

# Matriz de Resultados (Proporção do Erro Tipo I)
matriz_wald <- xtabs(taxa_erro_tipo1 ~ N + n, data = df_resultados_wald)
print(matriz_wald)
