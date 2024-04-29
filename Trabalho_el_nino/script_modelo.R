# Carregando os pacotes necessários

library("GGally")     # grafico em matriz
library("MASS")		    # stepwise (aic)
library("mixlm")		  # stepwise (valor-p)
#library("glmulti")	  # todas as regressoes
library("tidyverse")	# manipulacao de dados
library("rvest")
library("reshape2")
library("dplyr")
library("ggplot2")


# Importando os dados a serem utilizados
url <- "https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php"
html <- read_html(url)
tabelas <- html |> html_nodes("table")

tabela_desejada <- tabelas[[9]]

dados <- tabela_desejada|> 
  html_table(header=T) |>
  filter(Year !="Year" & Year != 2024)

dados <-  pivot_longer(data = dados, cols = colnames(dados)[2:ncol(dados)],
               names_to = "meses") |>
  filter(meses %in% c("JFM", "AMJ", "JAS", "OND"))

dados <- dados |>
  summarise(
    ano = as.numeric(Year),
    trimestre = case_when(
    meses == "JFM" ~ 1,
    meses == "AMJ" ~ 2,
    meses == "JAS" ~ 3,
    meses == "OND" ~ 4
  ),
  
  fenomeno = case_when(
    value >= 0.5 ~ "nino",
    value <= -0.5 ~ "nina",
    TRUE ~ "neutro"
  )
  )

pesca <- read.delim("~/Analise_de_regressao/Trabalho_el_nino/pesca.txt")

dados2 <- left_join(pesca, dados, by = c("ano", "trimestre")) |>
  mutate(
    cpue2 = log(cpue)
    ) |>
  select(-c(dias_pesca, captura))


#### Análise Descritiva

#Histograma de cpue

ggplot(dados2, aes(x = cpue)) +
  geom_density(alpha = 0.5, fill = "violet", col = "violet") +  
  labs(
    #title = "Histograma de Densidade",
    x = "Captura por unidade de pesca",
    y = "Densidade"
  ) + theme_minimal() + 
  geom_vline(xintercept = mean(dados2$cpue), linetype = "dashed", color = "black")

# Box plot fenomeno x cpue

ggplot(dados2, aes(x = fenomeno, y = cpue, fill = fenomeno)) +
  geom_boxplot() +
  labs(
   # title = " ",
    x = "Fenômeno natural",
    y = "Captura por unidade de esforço"
  ) + scale_x_discrete(labels = c("Nenhum", "La Niña", "El Niño")) +
  theme_minimal()

# Densidade de cpue agrupado por ano

ggplot(dados2, aes(x = cpue, fill = as.factor(ano))) +
  geom_density(alpha = 0.3) +
  labs(
    title = " ",
    x = "Captura por unidade de esforço",
    y = "Densidade",
    fill = 'Ano'
  ) + theme_minimal() 
 

# Densidade de cpue agrupado por frota

ggplot(dados2, aes(x = cpue, fill = frota)) +
  geom_density(alpha = 0.3) +
  labs(
    title = " ",
    x = "Captura por unidade de esforço",
    y = "Densidade",
    fill = 'Frota'
  ) + theme_minimal() 

# Box plot trimestre x cpue

ggplot(dados2, aes(x = trimestre, y = cpue, fill = as.factor(trimestre))) +
  geom_boxplot() +
  labs(
    # title = " ",
    x = "Trimestre",
    y = "Captura por unidade de esforço",
    fill = "Trimestre"
  ) + 
  theme_minimal()

#### Correlação

ggcorr(select(dados2, -c(cpue2)), label = T)

ggcorr(select(dados2, -c(cpue2)), geom = "blank", label = TRUE, hjust = 0.75) +
  geom_point(size = 10, aes(color = coefficient >= 0, alpha = abs(coefficient) >= 0.05)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = FALSE, alpha = FALSE)

#### Regressão

dados2 <- dados2 |> select(-c(cpue))

modelo <- stats::lm(cpue2 ~ ., data=dados2)
summary(modelo)

opt_model_step_aic<- stepAIC(modelo, direction="both") 
summary(opt_model_step_aic)
