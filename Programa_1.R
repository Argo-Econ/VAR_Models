# ------------------------------------------------------------------------------#
# Introducción a modelos multivariados ----
## Vectores autoregresivos VAR(p)
## Arturo Yesid Gonzalez ----
# ------------------------------------------------------------------------------#

# carga de librerias (funciones para el desarrollo practico)

pacman::p_load(readxl, sandwich, car, lmtest, TSstudio, lmtest, forecast
               , vars, urca, tsDyn, tidyverse, glue, xts, tseries)


# ------------------------------------------------------------------------------#
# Carga de paths ----

entradas <- "Datos_ent/"
salidas <- "Datos_sal/"


# ------------------------------------------------------------------------------#
# Lectura datos ----
Datos_ent  <- read_xlsx(glue("{entradas}Base.xlsx")
                        ,sheet = "Base",range = "a3:e321")

Datos_exo  <- read_xlsx(glue("{entradas}Base.xlsx")
                        ,sheet = "Exogenas",range = "a3:c435")

Datos_ent_xts <- xts(Datos_ent[,-1],order.by = as.Date(Datos_ent$fecha))
Datos_exo_xts <- xts(Datos_ent[,-1],order.by = as.Date(Datos_ent$fecha))

# Gráficos ----
#------------------------------------------------------------------------------#
ts_plot(Datos_ent_xts
        ,type = "multiple"
        ,slider = F)

# Test estacionariedad ----

## TRM
adf.test(Datos_ent_xts[,1])
kpss.test(Datos_ent_xts[,1]) # H0: serie estacionaria
pp.test(Datos_ent_xts[,1])   # H0: serie no estacionaria

## Brent
adf.test(Datos_ent_xts[,2])
kpss.test(Datos_ent_xts[,2]) # H0: serie estacionaria
pp.test(Datos_ent_xts[,2])   # H0: serie no estacionaria

## DXY
adf.test(Datos_ent_xts[,3])
kpss.test(Datos_ent_xts[,3]) # H0: serie estacionaria
pp.test(Datos_ent_xts[,3])   # H0: serie no estacionaria

## EMBI
adf.test(Datos_ent_xts[,4])
kpss.test(Datos_ent_xts[,4]) # H0: serie estacionaria
pp.test(Datos_ent_xts[,4])   # H0: serie no estacionaria


# diff aplicar ----
ndiffs(Datos_ent_xts[,1])
ndiffs(Datos_ent_xts[,2])
ndiffs(Datos_ent_xts[,3])
ndiffs(Datos_ent_xts[,4])


# Calculo retornos log. ----
Datos_ent_xts_dlx <- Datos_ent_xts %>% log() %>% diff(.,lag=1, differences = 1 ) %>% na.omit()
Datos_ent_xts_slx <- Datos_ent_xts %>% log() %>% diff(.,lag=12, differences = 1 ) %>% na.omit()

Datos_exo_xts_dlx <- Datos_exo_xts %>% log() %>% diff(.,lag=1, differences = 1 ) %>% na.omit()
Datos_exo_xts_slx <- Datos_exo_xts %>% log() %>% diff(.,lag=12, differences = 1 ) %>% na.omit()


# Gráficos ----
#------------------------------------------------------------------------------#
ts_plot(Datos_ent_xts_slx
        ,type = "single"
        ,slider = F)


# Las series de retornos log mensuales son I(0)
## TRM
adf.test(Datos_ent_xts_dlx[,1])
kpss.test(Datos_ent_xts_dlx[,1]) # H0: serie estacionaria
pp.test(Datos_ent_xts_dlx[,1])   # H0: serie no estacionaria

# Las series de retornos log anuales son I(0)
## TRM
adf.test(Datos_ent_xts_slx[,1])
kpss.test(Datos_ent_xts_slx[,1]) # H0: serie estacionaria
pp.test(Datos_ent_xts_slx[,1])   # H0: serie no estacionaria


# Correlaciones dinámicas ----
#------------------------------------------------------------------------------#

test <- ccf(x = as.numeric(Datos_ent_xts_slx$EMBI) ,y = as.numeric(Datos_ent_xts_slx$TRM)
                      ,lag.max = 18,plot = T,type = "correlation")



# Identificación ----
#------------------------------------------------------------------------------#

## retornos log. mensuales
VARselect(Datos_ent_xts_dlx,lag.max = 10,type = "both")

## retornos log. anuales
VARselect(Datos_ent_xts_slx,lag.max = 10,type = "const")



# Estimación ----
#------------------------------------------------------------------------------#

mod1 <- VAR(y = Datos_ent_xts_dlx, p = 1, type = "both")
summary(mod1)

summary(mod1, equation = "TRM")
windows()
plot(mod1, names = "TRM")

mod2 <- VAR(y = Datos_ent_xts_slx, lag.max = 3, type = "both", ic = "AIC")
summary(mod2)

summary(mod2, equation = "TRM")
windows()
plot(mod2, names = "TRM")


# Diagnostico ----
#------------------------------------------------------------------------------#

roots(mod1) # los valores propios deben ser inferiores a la unidad para que el
            # sistema VAR estimado sea estable

## correlación serial ----
serial_mod1 <- serial.test(mod1,lags.pt = 16, type = "PT.asymptotic")

windows()
plot(serial_mod1, names= "TRM")


serial_mod2 <- serial.test(mod2,lags.pt = 16, type = "PT.asymptotic")

windows()
plot(serial_mod2, names= "TRM")

## Normalidad ----
normal_mod1 <- normality.test(mod1, multivariate.only = "TRUE")

normal_mod2 <- normality.test(mod2, multivariate.only = "TRUE")


## Heterocedasticidad ----
hetero_mod1 <- arch.test(mod1, multivariate.only = "TRUE")

hetero_mod2 <- arch.test(mod2, multivariate.only = "TRUE")


## Estabilidad ----
windows()
plot(stability(mod1, type = "fluctuation"))

windows()
plot(stability(mod2, type = "fluctuation"))

windows()
plot(stability(mod2, type = "OLS-CUSUM"))



#------------------------------------------------------------------------------#
# Prueba alterna Ljung Box Q* ----

#------------------------------------------------------------------------------#
# Prueba de Ljung Box sobre cada residual de las ecuaciones del modelo
# ningun p-value de cada rezago debe estar por debajo de la linea roja
#------------------------------------------------------------------------------#
names(Datos_ent_xts_slx)

resi <- residuals(mod2)
i=1
k=1

for(k in 1:ncol(resi)){
  Ljung_Box <- NULL
  
  for(i in 1:18){
    Ljung_Box1 <- Box.test(resi[,k], lag = i, type="Ljung")$p.value
    Ljung_Box <- rbind(Ljung_Box, cbind(i,Ljung_Box1))
  }
  colnames(Ljung_Box) <- c("Rezago","p-valor");
  cat(" \n \n")
  cat("Prueba Ljung_Box \n \n")
  Ljung_Box
  windows()
  plot(Ljung_Box, main=paste("Prueba Ljung_Box \n H0: residuales ecuación "
                             ,names(Datos_ent_xts_slx)[k]," son iid",sep=""))
  abline(h=0.05,col="red")
}


#------------------------------------------------------------------------------#
# Uso de modelos ----

## Pronóstico ----
pronos_mod1 <- predict(mod1, n.ahead = 15, ci=0.95)
windows()
plot(pronos_mod1)

windows()
fanchart(pronos_mod1)

pronos_mod2 <- predict(mod2, n.ahead = 15, ci=0.95)
windows()
plot(pronos_mod2,names="TRM")

windows()
fanchart(pronos_mod2,names="TRM")


## Causalidad Granger ---- "Si una variable aporta sobre el pronóstico de mi variable de interés"
granger_mod1 <- causality(mod1,cause = "DXY")
granger_mod1

granger_mod1 <- causality(mod1,cause = "Brent")
granger_mod1

granger_mod2 <- causality(mod2,cause = "DXY")
granger_mod2

granger_mod2 <- causality(mod2,cause = "Brent")
granger_mod2


## Impulso respuesta IRF ----
irf_mod1 <- irf(mod1,impulse = c("TRM","Brent","DXY","EMBI")
                ,response = "TRM", cumulative = T,
                ,n.ahead = 24)
windows()
plot(irf_mod1)

irf_mod2 <- irf(mod2,impulse = c("TRM","Brent","DXY","EMBI")
                ,response = "TRM", cumulative = F,ortho=T
                ,boot=T,ci=0.95,n.ahead = 24)
windows()
plot(irf_mod2)


## FEVD analisis ----
## Descomposición de la varianza de pronóstico

fevd_mod1 <- vars::fevd(mod1, n.ahead = 12)

windows()
plot(fevd_mod1)

fevd_mod2 <- vars::fevd(mod2, n.ahead = 18)

windows()
plot(fevd_mod2,names="TRM")


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#



