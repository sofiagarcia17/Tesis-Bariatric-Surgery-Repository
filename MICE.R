library(mice)
library(VIM)  
library(parallel)
library(doParallel)
library(naniar)
library(nnet)
#==============================preprocesado====================================#

library(readr)
datos <- read_csv("BARCO_base_def.csv")

# Ver primeras filas
head(datos)

summary(datos)

# Detectar variables de tipo character y convertir a factor
datos[] <- lapply(datos, function(x) {
  if (is.character(x)) as.factor(x) else x
})

summary(datos)


# Resumen de valores faltantes
summary(datos)
md.pattern(datos)  # patrón de NA
aggr(datos, numbers = TRUE, sortVars = TRUE,
     labels = names(datos), cex.axis = 0.7, gap = 3,
     ylab = c("Histograma de NA", "Patrón de NA"))

# Gráfico de vacíos
vis_miss(datos, cluster = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # mapa de missing (naniar)

# Variables con >50% faltantes (CORREGIDO)
datos %>% 
  select(where(~ mean(is.na(.x)) > 0.5)) %>% 
  names()

# 2. Porcentaje exacto de faltantes por variable
miss_var_summary(datos) %>% 
  filter(pct_miss > 50) %>% 
  arrange(desc(pct_miss)) %>% 
  mutate(pct_miss = round(pct_miss, 1)) %>% 
  knitr::kable(col.names = c("Variable", "Faltantes", "% Faltantes"))


#verificar el patron de datos faltantes
# Little's MCAR
mcar_test(datos)
# No hay evidencia estadísticamente significativa para decir que los datos faltantes siguen un patrón distinto de completamente al azar.
# MCAR es el escenario más favorable para imputar.

#==============================Metodos MICE====================================#

# Detectar tipos
str(datos)

str(datos$GÉNERO)

metodos <- make.method(datos)

# Asignar métodos de imputación según tipo de variable
for (v in names(datos)) {
  if (is.factor(datos[[v]]) && length(levels(datos[[v]])) == 2) {
    # Variable categórica binaria
    metodos[v] <- "logreg"
  } else if (is.factor(datos[[v]]) && length(levels(datos[[v]])) > 2) {
    # Variable categórica con más de 2 niveles
    metodos[v] <- "polyreg"
  } else if (is.numeric(datos[[v]])) {
    # Variable numérica
    metodos[v] <- "pmm"
  }
}

metodos

#====utilizar el anterior o este, deberia ser lo mismo====#

# Configurar método de imputación por tipo de variable:
#metodo_imputacion <- c(
#  "numeric" = "pmm",          # Predictive Mean Matching (continuas)
#  "factor" = "polyreg",       # Regresión politómica (categóricas)
#  "binary" = "logreg"         # Regresión logística (binarias)
#)



# Paso 1: Calcular porcentaje de datos faltantes por variable
missing_pct <- sapply(datos, function(x) mean(is.na(x)))

# Paso 3: Desactivar imputación para variables con más de 50% datos faltantes
metodos[missing_pct > 0.5] <- ""

# Variables con <50% faltantes (CORREGIDO)
datos %>% 
  select(where(~ mean(is.na(.x)) < 0.5)) %>% 
  names()

# Variables con menos del 50% de datos faltantes

# Porcentaje exacto de faltantes por variable
miss_var_summary(datos) %>% 
  filter(pct_miss < 50) %>% 
  arrange(desc(pct_miss)) %>% 
  mutate(pct_miss = round(pct_miss, 1)) %>% 
  knitr::kable(col.names = c("Variable", "Faltantes", "% Faltantes"))
#==============================MICE====================================#

set.seed(123) # para reproducibilidad

# 3. Preparación para alta dimensionalidad (reducir predictores)
pred_matrix <- quickpred(datos, 
                         minpuc = 0.5,   # Mínimo 50% de casos completos
                         mincor = 0.3 )  # Correlación mínima de 0.3
                         #include = c("tabaquismo", "EDAD", "GÉNERO"))  # Variables clave


# Configuración de procesamiento paralelo (acelera cálculo)
cl <- makeCluster(detectCores() - 1)  # Usa todos los núcleos menos uno
registerDoParallel(cl)


imp <- mice(datos, 
            method = metodos, # tu vector de métodos
            predictorMatrix = pred_matrix, 
            m = 5, # número de datasets imputados
            maxit = 50, # número de iteraciones
            seed = 123)# reproducibilidad

# Detener cluster paralelo
stopCluster(cl)


imp
summary(imp)

# Gráfico de convergencia
plot(imp)

# Ver imputaciones para una variable específica
imp$imp$TABACO   # ejemplo: ver imputaciones para tabaco

stripplot(imp, pch = 20, cex = 1.2)  # visualiza distribuciones imputadas vs observadas

#densityplot(imp)  # densidades por variable imputada
# Tarda mucho

#====================== Guardar conjunto ======================================#

# Extraer dataset imputado (primer conjunto)
datos_imputados <- complete(imp, 1)

# Guardar resultados
write.csv(datos_imputados, "datos_imputados.csv", row.names = FALSE)

# Gráfico de vacíos
vis_miss(datos_imputados, cluster = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # mapa de missing (naniar)

#==================== Verificar ===============================================#



library(ggplot2)

# Crear dataframe con ambas versiones (conserva longitud original)
df_comparacion <- data.frame(
  CC_original = datos$AC_URICO,  # Conserva NAs
  CC_imputado = datos_imputados$AC_URICO
)

# Gráfico corregido
ggplot(df_comparacion) +
  geom_density(aes(x = CC_original, fill = "Original"), alpha = 0.5, na.rm = TRUE) +
  geom_density(aes(x = CC_imputado, fill = "Imputado"), alpha = 0.5) +
  scale_fill_manual(values = c("Original" = "blue", "Imputado" = "red")) +
  labs(title = "Distribución de ácido úrico",
       x = "mg/dL",
       y = "Densidad",
       fill = "Grupo") +
  theme_minimal()

# Test KS para distribuciones (variables continuas)
ks.test(na.omit(datos$AC_URICO), datos_imputados$AC_URICO)


#==================== Graficos ===============================================#

# Gráfico de vacíos antes
vis_miss(datos, cluster = TRUE) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) # mapa de missing (naniar)

# Gráfico de vacíos despues
vis_miss(datos_imputados, cluster = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # mapa de missing (naniar)

library(ggplot2)

# Función para graficar y hacer KS test para una variable
comparar_imputacion <- function(var_name) {
  df_comp <- data.frame(
    original = datos[[var_name]],
    imputado = datos_imputados[[var_name]]
  )
  
  # Gráfico
  p <- ggplot(df_comp) +
    geom_density(aes(x = original, fill = "Original"), alpha = 0.5, na.rm = TRUE) +
    geom_density(aes(x = imputado, fill = "Imputado"), alpha = 0.5) +
    scale_fill_manual(values = c("Original" = "blue", "Imputado" = "red")) +
    labs(title = paste("Distribución de", var_name, "Original vs imputado"),
         x = var_name,
         y = "Densidad",
         fill = "Grupo") +
    theme_minimal()
  
  print(p)
  
  # Test KS
  ks <- ks.test(na.omit(datos[[var_name]]), datos_imputados[[var_name]])
  print(ks)
}

# Variables a comparar
variables <- c("AC_URICO","Ca", "LDL", "PROT_TOTAL", "ALBUMINA")

# Ejecutar para cada variable
for (v in variables) {
  comparar_imputacion(v)
}

############
# 1. AC_URICO
df_AC_URICO <- data.frame(
  original = datos$AC_URICO,
  imputado = datos_imputados$AC_URICO
)

ks_AC_URICO <- ks.test(na.omit(datos$AC_URICO), datos_imputados$AC_URICO)
pval_AC_URICO <- signif(ks_AC_URICO$p.value, 3)

ggplot(df_AC_URICO) +
  geom_density(aes(x = original, fill = "Original"), alpha = 0.5, na.rm = TRUE) +
  geom_density(aes(x = imputado, fill = "Imputado"), alpha = 0.5) +
  scale_fill_manual(values = c("Original" = "blue", "Imputado" = "red")) +
  labs(title = "Distribución de ácido úrico Original vs Imputado",
       x = "mg/dL",
       y = "Densidad",
       fill = "Grupo") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste0("KS test p = ", pval_AC_URICO),
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")

print(ks_AC_URICO)


# 1. Ca
df_Ca <- data.frame(
  original = datos$Ca,
  imputado = datos_imputados$Ca
)

ks_Ca <- ks.test(na.omit(datos$Ca), datos_imputados$Ca)
pval_Ca <- signif(ks_Ca$p.value, 3)

ggplot(df_Ca) +
  geom_density(aes(x = original, fill = "Original"), alpha = 0.5, na.rm = TRUE) +
  geom_density(aes(x = imputado, fill = "Imputado"), alpha = 0.5) +
  scale_fill_manual(values = c("Original" = "blue", "Imputado" = "red")) +
  labs(title = "Distribución de Calcio Original vs Imputado",
       x = "mg/dL",
       y = "Densidad",
       fill = "Grupo") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste0("KS test p = ", pval_Ca),
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")

print(ks_Ca)

# 2. LDL
df_LDL <- data.frame(
  original = datos$LDL,
  imputado = datos_imputados$LDL
)

ks_LDL <- ks.test(na.omit(datos$LDL), datos_imputados$LDL)
pval_LDL <- signif(ks_LDL$p.value, 3)

ggplot(df_LDL) +
  geom_density(aes(x = original, fill = "Original"), alpha = 0.5, na.rm = TRUE) +
  geom_density(aes(x = imputado, fill = "Imputado"), alpha = 0.5) +
  scale_fill_manual(values = c("Original" = "blue", "Imputado" = "red")) +
  labs(title = "Distribución de LDL Original vs Imputado",
       x = "mg/dL",
       y = "Densidad",
       fill = "Grupo") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste0("KS test p = ", pval_LDL),
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")

print(ks_LDL)

# 3. PROT_TOTAL
df_PROT_TOTAL <- data.frame(
  original = datos$PROT_TOTAL,
  imputado = datos_imputados$PROT_TOTAL
)

ks_PROT_TOTAL <- ks.test(na.omit(datos$PROT_TOTAL), datos_imputados$PROT_TOTAL)
pval_PROT_TOTAL <- signif(ks_PROT_TOTAL$p.value, 3)

ggplot(df_PROT_TOTAL) +
  geom_density(aes(x = original, fill = "Original"), alpha = 0.5, na.rm = TRUE) +
  geom_density(aes(x = imputado, fill = "Imputado"), alpha = 0.5) +
  scale_fill_manual(values = c("Original" = "blue", "Imputado" = "red")) +
  labs(title = "Distribución de proteina Original vs Imputado",
       x = "g/dL",
       y = "Densidad",
       fill = "Grupo") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste0("KS test p = ", pval_PROT_TOTAL),
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")

print(ks_PROT_TOTAL)

###===============================####Experimento#####===============================#####

# Modelo predictivo multinomial en imputaciones
fit_multinom <- with(
  imp,
  nnet::multinom(COMPLICACIONES ~ EDAD + IMC + PESO + LDL + ALBUMINA)
)

# Combinar resultados con pooling de Rubin
resumen <- pool(fit_multinom)

summary(resumen)

# Extraer una imputación y calcular R² aproximado
library(DescTools)
modelo_ejemplo <- multinom(COMPLICACIONES ~ COLECISTECTOMIA+CIELO+HTA+SAOS+REF_SUT_MANGA+
                           HELYCOBACTER+TRANSITO+OCUPACIÓN+ESTADO_CIVIL+OTRA+EXTRAOPERATORIAS+
                           ECOABD_NUEVA+ESPIROMETRÍA+EDAD+TALLA+ERITROCITOS+RCTO_LEUCOCITOS+RCTO_PLAQUETAS+
                           ALBUMINA+PÉPTIDO_C+BUN+CREATINEMIA+Ca+K+Cl,
                           data = complete(imp, 1))

PseudoR2(modelo_ejemplo, which = "Nagelkerke")

table(complete(imp, 1)$COMPLICACIONES, useNA = "ifany")


datos1 <- complete(imp, 1)

# Recodificación
datos1$COMPLICACIONES_BIN <- ifelse(datos1$COMPLICACIONES == "NO", "NO", "SI")
datos1$COMPLICACIONES_BIN <- factor(datos1$COMPLICACIONES_BIN)

table(datos1$COMPLICACIONES_BIN)

# Modelo logístico binario
modelo_bin <- glm(COMPLICACIONES_BIN ~ COLECISTECTOMIA + CIELO + HTA + SAOS + REF_SUT_MANGA +
                    HELYCOBACTER + TRANSITO + OCUPACIÓN + ESTADO_CIVIL + OTRA + EXTRAOPERATORIAS +
                    ECOABD_NUEVA + ESPIROMETRÍA + EDAD + TALLA + ERITROCITOS + RCTO_LEUCOCITOS + 
                    RCTO_PLAQUETAS + ALBUMINA + PÉPTIDO_C + BUN + CREATINEMIA + Ca + K + Cl,
                  data = datos1,
                  family = binomial)   # <- clave para modelo logístico

# Predicciones en una imputación
pred <- predict(modelo_ejemplo, newdata = complete(imp, 1))

table(Predicho = pred, Real = complete(imp, 1)$COMPLICACIONES)
