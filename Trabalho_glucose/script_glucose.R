library("GGally")     # grafico em matriz
library("MASS")		    # stepwise (aic)
library("mixlm")		  # stepwise (valor-p)
library("tidyverse")	# manipulacao de dados
library("rvest")
library("reshape2")
library("dplyr")
library("ggplot2")
library("readr")
library("janitor")
library("plotly")
library(car)    	  # vif - multicolinearidade
library(nortest)  	# normalidade
library(lmtest)		  # homocedasticidade e auto-correlação
library(gamlss)		  # incorporando heterocedasticidade
library(nlme)		    # incorporando auto-correlação


# IMPORTANDO DADOS

df <- read_csv("Trabalho_glucose/glicose.csv") |>
  clean_names()

# ANALISANDO DADOS FALTANTES

unique(df$pregnancies) #certo
unique(df$glucose) #nível de glicose zerado?
unique(df$blood_pressure) #pressão zerada
unique(df$skin_thickness) #espessura zerada
unique(df$insulin) #nível de insulina zerado
unique(df$bmi) #imc zerado
unique(df$diabetes_pedigree_function) #certo
unique(df$age) #certo
unique(df$outcome) #certo

# Verificando o número de dados faltantes para cada coluna

glucose0 <- nrow(filter(df, glucose == 0)) #5 dados faltantes
pressure0 <- nrow(filter(df, blood_pressure == 0)) # 35 dados faltantes
skin0 <- nrow(filter(df, skin_thickness == 0)) #227 dados faltantes
insulin0 <- nrow(filter(df, insulin == 0)) #374 dados faltantes
bmi0 <- nrow(filter(df, bmi==0)) #11 dados faltantes

# Retirando os dados faltantes

df <- df |>
  select(-c(insulin, skin_thickness)) |>
  filter(glucose != 0 & blood_pressure != 0 & bmi != 0)

# ANÁLISE DESCRITIVA


# Histograma da variável 'glucose'

ggplot(df, aes(x = glucose)) +
  geom_density(alpha = 0.5, fill = "violet", col = "violet") +  
  labs(
    #title = "Histograma de Densidade",
    x = "Captura por unidade de pesca",
    y = "Densidade"
  ) + theme_minimal() + 
  geom_vline(xintercept = mean(df$glucose), linetype = "dashed", color = "black")

# Gráfico de dispersão glucose x blood_pressure

ggplot(data = df, aes(x = glucose, y = blood_pressure)) +
  geom_point(col = "violet") +  
  labs(x = "Concentração de glicose", y = "Pressão arterial", title = " ")  + theme_minimal()

# Gráfico de dispersão glucose x bmi

ggplot(data = df, aes(x = glucose, y = bmi)) +
  geom_point(col = "violet") +  
  labs(x = "Concentração de glicose", y = "Índice de massa corporal", title = " ")  + theme_minimal()

# Gráfico de dispersão glucose x diabates_pedigree_function

ggplot(data = df, aes(x = glucose, y = diabetes_pedigree_function)) +
  geom_point(col = "violet") +  
  labs(x = "Concentração de glicose", y = "Diabetes função da genealogia", title = " ")  + theme_minimal()

# Gráfico densidade glucose para os diabéticos e não diabéticos

ggplot(df, aes(x = glucose, fill = as.factor(outcome))) +
  geom_density(alpha = 0.3) +
  labs(
    title = " ",
    x = "Concentração de glicose",
    y = "Teste de diabetes",
    fill = 'Teste de diabetes'
  ) + theme_minimal() + 
  scale_fill_manual(values = c("0" = "violet", "1" = "grey"),
                     labels = c("Saudável", "Diabético"))

# Box plot glucose por idade discretizada

df2 <- df
df2 <- df2|>
  mutate(
    age_interval = case_when(
      age <= 29 & age >= 21 ~ "21 a 29",
      age <= 35 & age >= 30 ~ "30 a 35",
      age <= 43 & age >= 36 ~ "36 a 43",
      age <= 50 & age >= 44 ~ "44 a 50",
      age <= 56 & age >= 51 ~ "51 a 56",
      age <= 63 & age >= 57 ~ "57 a 63",
      age <= 70 & age >= 64 ~ "64 a 70",
      age <= 75 & age >= 71 ~ "71 a 75",
      age <= 81 & age >= 76 ~ "76 a 81",
    ),
    
    preg_interval = case_when(
      pregnancies == 0 ~ "Sem gestações",
      pregnancies >= 1 & pregnancies <= 2 ~ "1 a 2",
      pregnancies >= 3 & pregnancies <= 5 ~ "3 a 5",
      pregnancies >= 10 & pregnancies <= 13 ~ "10 a 13",
      pregnancies >= 6 & pregnancies <= 9 ~ "6 a 9",
      pregnancies >= 14 & pregnancies <= 17 ~ "14 a 17"
      )
  )

ggplot(df2, aes(x = age_interval, y = glucose, fill = age_interval)) +
  geom_boxplot() +
  labs(
    x = "Idade em anos",
    y = "Concentração de glicose"
  ) + theme_minimal()


# Gráfico de violino glucose por número de gestações discretizada

df2$preg_interval <- factor(df2$preg_interval, levels = c("Sem gestações", "1 a 2", "3 a 5", "6 a 9", "10 a 13", "14 a 17"))


ggplot(df2, aes(x = preg_interval, y = glucose, fill = preg_interval)) +
  geom_violin() +
  labs(
    title = " ",
    x = "Número de gestações",
    y = "Concentração de glicose",
    fill = 'Número de gestações'
  ) +
  theme_minimal()

# ANÁLISE DE CORRELAÇÃO

ggcorr(df, geom = "blank", label = TRUE, hjust = 0.75) +
  geom_point(size = 10, aes(color = coefficient >= 0, alpha = abs(coefficient) >= 0.05)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = FALSE, alpha = FALSE)

# MODELO 

modelo <- stats::lm(glucose ~ ., data=df)
summary(modelo)

opt_model_step_aic<- stepAIC(modelo, direction="both") 
summary(opt_model_step_aic)

# PONTOS INFLUENTES/DE ALAVANCA/OUTLIERS

fit<- opt_model_step_aic

n<- nrow(df)    		        # número de observações
k<- length(fit$coef) 		    # k=p+1 (número de coeficientes)

corte.hii<- 2*k/n		        # corte para elementos da diagonal de H
corte.cook<- qf(0.5,k,n-k)	# corte para Distância de Cook
corte.stu<- 2			          # corte para resíduos estudentizados

rst<- rstudent(fit)		      # resíduos estudentizados
hii<- hatvalues(fit) 		    # valores da diagonal da matriz H
dcook<- cooks.distance(fit)	# distância de Cook

obs<- 1:n

df.fit<- data.frame(obs,rst,hii,dcook)



# GRÁFICO - RESÍDUOS ESTUDENTIZADOS

df.fit %>% ggplot(aes(x=obs,y=rst)) + 
  geom_point() + 
  geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Resíduo Estudentizado") + 
  theme_bw()



# GRÁFICO - ALAVANCAGEM

df.fit %>% ggplot(aes(x=obs,y=hii,ymin=0,ymax=hii)) + 
  geom_point() + 
  geom_linerange() + 
  geom_hline(yintercept = corte.hii, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Alavancagem") + 
  theme_bw()



# GRÁFICO - DISTÂNCIA DE COOK

df.fit %>% ggplot(aes(x=obs,y=dcook,ymin=0,ymax=dcook)) + 
  geom_point() + 
  geom_linerange() +
  geom_hline(yintercept = corte.cook, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Distância de Cook") + 
  theme_bw()



# GRÁFICO GERAL DE DIAGNÓSTICO

texto<- paste("observ:",obs,"\n",
              "resid_stud:",round(rst,2),"\n",
              "alavancagem:",round(hii,2),"\n",
              "D_Cook:",round(dcook,2),"\n",
              "corte_Cook:",round(corte.cook,2))

ggplotly(
  
  df.fit %>% ggplot(aes(x=hii,y=rst,text=texto)) +
    geom_point(aes(size=dcook)) + 
    xlim(0, max(max(hii),corte.hii)) + 
    ylim(-max(abs(rst),corte.stu), max(abs(rst),corte.stu)) + 
    geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
    geom_vline(xintercept = corte.hii, color="red", linetype="dashed") + 
    theme_bw() +
    theme(legend.position="none") + 
    xlab("Alavancagem") + 
    ylab("Resíduo Estudentizado"),
  tooltip = c("text")
  
)


# ANALISANDO OS PRESSUPOSTOS

# NORMALIDADE 
# HISTOGRAMA

rst %>% data.frame() %>% ggplot(aes(x=rst)) + 
  geom_histogram(aes(y=..density..)) + 
  geom_density(alpha=.1, fill="blue") +
  theme_bw()



# GRÁFICO QUANTIL-QUANTIL

rst %>% data.frame() %>% ggplot(aes(sample=rst)) + 
  stat_qq() + 
  stat_qq_line() +
  theme_bw()

# TESTE DE NORMALIDADE

t1 <- ks.test(rst,"pnorm")	#KS
t2 <- lillie.test(rst)		  # Lilliefors
t3 <- cvm.test(rst)		      # Cramér-von Mises
t4 <- shapiro.test(rst)		  # Shapiro-Wilk
t5 <- sf.test(rst)		      # Shapiro-Francia
t6 <- ad.test(rst)		      # Anderson-Darling

# Tabela de resultados
testes <- c(t1$method, t2$method, t3$method, t4$method, t5$method,t6$method)
estt <- as.numeric(c(t1$statistic,
                     t2$statistic,
                     t3$statistic,
                     t4$statistic,
                     t5$statistic,
                     t6$statistic))
valorp <- c(t1$p.value, t2$p.value, t3$p.value, t4$p.value, t5$p.value,t6$p.value)
resultados <- cbind(estt, valorp)
rownames(resultados) <- testes
colnames(resultados) <- c("Estatística", "p")
round(resultados, digits = 4)



# TRANSFORMAÇÃO BOXCOX

transf_boxcox<- fit %>% boxcox(data=df)
lambda<- transf_boxcox$x[which.max(transf_boxcox$y)]
lambda


# HOMOSCEDASTICIDADE #

# VERIFICANDO A HETEROSCEDASTICIDADE DOS ERROS

bptest(fit)			      # teste de Breusch-Pagan

# H0: homoscesticidade
# H1: heteroscedasticidade

# CASO HAJA HETEROSCEDASTICIDADE DEVEMOS:
# (1) TRANSFORMAR OS DADOS PARA ELIMINA-LA; OU
# (2) INCORPORA-LA AO MODELO COM OS COMANDOS ABAIXO

fit.gamlss<- gamlss(glucose ~ ., sigma.formula = ~ ., data=df)
fit.gamlss %>% summary()

# ERROS NÃO-CORRELACIONADOS

# VERIFICANDO SE HÁ CORRELAÇÃO SERIAL

bgtest(fit)			      # teste de Breusch-Godfrey

# H0: erros não correlacionados
# H1: erros correlacionados

# CASO HAJA AUTOCORRELAÇÃO DEVEMOS:
# INCORPORA-LA AO MODELO COM OS COMANDOS ABAIXO

fit.gls<- gls(y ~ ., data=df, correlation=corARMA(p=1))	# AR(1)
fit.gls %>% summary()

