# Carregando os pacotes necessários

library("GGally")     # grafico em matriz
library("MASS")		    # stepwise (aic)
library("mixlm")		  # stepwise (valor-p)
library("glmulti")	  # todas as regressoes
library("tidyverse")	# manipulacao de dados
library("rvest")
library("reshape2")


# Importando os dados a serem utilizados
url <- "https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php"
html <- read_html(url)
tabelas <- html |> html_nodes("table")

tabela_desejada <- tabelas[[9]]

dados <- tabela_desejada|> 
  html_table(header=T) |>
  filter(Year !="Year" & Year != 2024) |>
  pivot_longer(cols = colnames(dados)[2:ncol(dados)],
               names_to = "trimestre")
