# Proyecto académico: Modelos de mezclas gaussianas (GMM) y Expectation-Maximization (EM)

## 1. Contexto y propósito del trabajo

El trabajo estudia los **Modelos de Mezclas Gaussianas (Gaussian Mixture
Models, GMM)** y el algoritmo **Expectation-Maximization (EM)** mediante
una aplicación práctica en RStudio.

Se consideran dos fuentes de datos:

1.  **USGS Earthquake Catalog**:
    https://earthquake.usgs.gov/data/comcat/
2.  **Inside Airbnb -- Get the Data**:
    https://insideairbnb.com/es/get-the-data/

La recomendación preliminar es utilizar **Inside Airbnb**,
preferiblemente los datos detallados (*Detailed Listings*) de Bogotá,
porque permiten construir perfiles de alojamientos a partir de variables
cuantitativas con interpretación económica y operativa. USGS sigue
siendo una alternativa válida, aunque exige más cuidado con la
estructura espacial, temporal y estadística de los eventos sísmicos.

El objetivo no es asumir de antemano cuántos grupos existen, sino
comprobar si una mezcla de distribuciones gaussianas representa
razonablemente la heterogeneidad de los datos y si EM permite estimar
componentes interpretables.

------------------------------------------------------------------------

## 2. Fundamentos conceptuales

### 2.1. ¿Qué es un modelo GMM?

Un modelo de mezclas gaussianas supone que los datos observados proceden
de una combinación de varias distribuciones normales (gaussianas).

La densidad del modelo es:

\[
p(x)=`\sum`{=tex}\_{k=1}\^{K}`\pi`{=tex}\_k`\mathcal{N}`{=tex}(x`\mid`{=tex}`\mu`{=tex}\_k,`\Sigma`{=tex}\_k)
\]

Donde:

-   (x): observación.
-   (K): número de componentes gaussianos.
-   (`\pi`{=tex}\_k): peso del componente (k), con (`\pi`{=tex}\_k
    `\geq 0`{=tex}) y (`\sum`{=tex}\_k `\pi`{=tex}\_k=1).
-   (`\mu`{=tex}\_k): vector de medias del componente.
-   (`\Sigma`{=tex}\_k): matriz de covarianzas del componente.
-   (`\mathcal{N}`{=tex}): densidad normal multivariada.

Cada componente representa una distribución gaussiana. Los pesos
expresan la contribución relativa de los componentes a la mezcla; no
tienen por qué coincidir exactamente con proporciones de grupos reales
conocidos.

### 2.2. Distribución normal multivariada

Para varias variables, se representa:

\[ X`\sim `{=tex}N(`\mu`{=tex},`\Sigma`{=tex}) \]

-   El vector (`\mu`{=tex}) contiene la media de cada variable.
-   La matriz (`\Sigma`{=tex}) contiene varianzas en la diagonal y
    covarianzas fuera de ella.

La covarianza permite modelar relaciones entre variables. Por ello, los
componentes pueden tener formas elípticas y diferentes orientaciones,
según la estructura de covarianza.

### 2.3. Clustering suave

GMM puede calcular la probabilidad posterior de que cada observación
pertenezca a cada componente. No obliga a que una observación tenga una
pertenencia exclusiva durante la estimación.

Por ejemplo, una observación podría tener estas probabilidades:

  Observación     Componente 1   Componente 2   Componente 3
  ------------- -------------- -------------- --------------
  A                       0.90           0.08           0.02
  B                       0.15           0.75           0.10
  C                       0.45           0.40           0.15

La clasificación final más probable puede asignarse al componente con
mayor probabilidad, pero conviene conservar las probabilidades para
analizar la ambigüedad.

### 2.4. ¿Qué es Expectation-Maximization?

EM es un algoritmo iterativo de estimación de parámetros que se utiliza
cuando hay variables latentes o información incompleta. En GMM, la
variable latente indica a qué componente pertenece cada observación.

EM alterna dos pasos:

**E-step (Expectation):** calcula las responsabilidades o probabilidades
posteriores de pertenencia:

\[ `\gamma`{=tex}*{ik}=
`\frac{\pi_k\mathcal{N}(x_i\mid\mu_k,\Sigma_k)}`{=tex}
{`\sum`{=tex}*{j=1}\^{K}`\pi`{=tex}\_j`\mathcal{N}`{=tex}(x_i`\mid`{=tex}`\mu`{=tex}\_j,`\Sigma`{=tex}\_j)}
\]

**M-step (Maximization):** actualiza los parámetros usando esas
responsabilidades.

Primero:

\[ N_k=`\sum`{=tex}*{i=1}\^{n}`\gamma`{=tex}*{ik} \]

Pesos:

\[ `\pi`{=tex}\_k=`\frac{N_k}{n}`{=tex} \]

Medias:

\[
`\mu`{=tex}*k=`\frac{1}{N_k}`{=tex}`\sum`{=tex}*{i=1}\^{n}`\gamma`{=tex}\_{ik}x_i
\]

Covarianzas:

\[
`\Sigma`{=tex}*k=`\frac{1}{N_k}`{=tex}`\sum`{=tex}*{i=1}\^{n}`\gamma`{=tex}\_{ik}(x_i-`\mu`{=tex}\_k)(x_i-`\mu`{=tex}\_k)\^T
\]

Después se repiten E-step y M-step hasta alcanzar un criterio de
convergencia, normalmente relacionado con el cambio de la
log-verosimilitud.

La log-verosimilitud observada es:

\[
`\ell`{=tex}(`\theta`{=tex})=`\sum`{=tex}*{i=1}\^{n}`\log`{=tex}`\left`{=tex}(`\sum`{=tex}*{k=1}\^{K}`\pi`{=tex}\_k`\mathcal{N}`{=tex}(x_i`\mid`{=tex}`\mu`{=tex}\_k,`\Sigma`{=tex}\_k)`\right`{=tex})
\]

EM normalmente no disminuye la verosimilitud en cada iteración, pero
puede converger a un máximo local o punto estacionario. La
inicialización importa.

### 2.5. Diferencias entre GMM y K-Means

  -----------------------------------------------------------------------
  Característica          K-Means                 GMM
  ----------------------- ----------------------- -----------------------
  Tipo                    Aprendizaje no          Aprendizaje no
                          supervisado             supervisado

  Asignación              Dura                    Probabilística durante
                                                  la estimación

  Representación          Centros/centroides      Distribuciones
                                                  gaussianas

  Covarianzas             No las modela           Puede estimarlas
                          explícitamente          

  Grupos                  Favorece grupos         Puede modelar
                          compactos y             componentes elípticos
                          aproximadamente         
                          esféricos               

  Pesos de mezcla         No son parámetros del   Se estiman
                          modelo                  

  Complejidad             Generalmente menor      Generalmente mayor
  -----------------------------------------------------------------------

K-Means también puede producir grupos de tamaños diferentes, pero no
estima explícitamente pesos probabilísticos ni matrices de covarianza.

------------------------------------------------------------------------

## 3. Cuándo utilizar GMM y qué datos son adecuados

GMM puede ser apropiado cuando:

-   No hay etiquetas previas y se busca descubrir estructura en los
    datos.
-   Las observaciones pueden representarse razonablemente mediante una
    mezcla de gaussianas.
-   Los grupos tienen dispersiones, tamaños u orientaciones diferentes.
-   Hay solapamiento entre grupos.
-   Interesa obtener probabilidades de pertenencia.
-   También interesa modelar densidades y detectar observaciones de baja
    densidad.

Características recomendables:

-   Variables numéricas continuas.
-   Suficientes observaciones para estimar los parámetros.
-   Tratamiento de datos faltantes y valores no finitos.
-   Atención a escalas diferentes.
-   Dimensionalidad razonable en relación con el tamaño de muestra.
-   Evaluación de valores extremos y distribuciones.

GMM no exige que la distribución global de cada variable sea normal: la
mezcla completa puede ser multimodal o asimétrica. El supuesto relevante
es que cada componente pueda aproximarse mediante una distribución
normal multivariada.

### 3.1. Tipos de covarianza

En `mclust` se pueden comparar diferentes estructuras de covarianza. En
términos generales:

-   **Spherical:** componentes esféricos; una varianza por componente.
-   **Diagonal:** varianzas diferentes por variable, sin covarianzas
    dentro del componente.
-   **Tied/compartida:** una estructura de covarianza común para los
    componentes.
-   **Full:** cada componente tiene su propia matriz de covarianza
    completa.

La nomenclatura exacta y las parametrizaciones disponibles dependen de
la implementación. `mclust` utiliza una familia de modelos geométricos
de covarianza, identificados mediante códigos de modelo.

### 3.2. Cuándo no es recomendable

Hay que considerar alternativas o adaptaciones cuando:

-   Los grupos tienen formas muy irregulares o no gaussianas.
-   Los datos son exclusivamente categóricos.
-   Hay muchas variables y pocas observaciones.
-   Existen valores extremos muy influyentes.
-   Las matrices de covarianza son singulares o inestables.
-   Los datos tienen dependencia temporal que el modelo estático no
    representa.

Alternativas posibles: DBSCAN/HDBSCAN, clustering espectral, modelos de
mezclas categóricas, mezclas t robustas, modelos dinámicos o HMM, según
el problema.

------------------------------------------------------------------------

## 4. Selección del número de componentes

El número de componentes (K) no debe fijarse sin exploración. Se pueden
comparar varios valores y estructuras de covarianza.

### 4.1. BIC y AIC

Una convención habitual del BIC es:

\[ BIC=-2`\ell`{=tex}(`\hat{\theta}`{=tex})+p`\log`{=tex}(n) \]

Con esta convención se prefiere el BIC más bajo.

Sin embargo, algunas implementaciones, entre ellas `mclust`, presentan
el BIC con signo opuesto:

\[ BIC\_{mclust}=2`\ell`{=tex}(`\hat{\theta}`{=tex})-p`\log`{=tex}(n) \]

En `mclust` se prefiere el BIC más alto. Es importante indicar la
convención al reportar resultados.

El AIC convencional es:

\[ AIC=-2`\ell`{=tex}(`\hat{\theta}`{=tex})+2p \]

AIC y BIC equilibran ajuste y complejidad mediante penalizaciones
distintas. Ninguno garantiza recuperar los grupos "verdaderos"; son
criterios para comparar modelos candidatos.

### 4.2. Evaluación complementaria

-   **Silhouette:** evalúa separación y cohesión después de asignar
    etiquetas duras. No aprovecha las probabilidades posteriores.
-   **Davies-Bouldin:** evalúa similitud entre agrupaciones duras.
-   **ARI:** compara con etiquetas reales si existen.
-   **Log-likelihood:** evalúa el ajuste probabilístico.
-   **Incertidumbre posterior:** permite estudiar ambigüedad en las
    asignaciones.
-   **Estabilidad:** evalúa si los componentes se mantienen ante
    distintas inicializaciones o muestras.

No se debe elegir un modelo únicamente por una métrica.

------------------------------------------------------------------------

## 5. Comparación de las fuentes de datos

### 5.1. USGS Earthquake Catalog

Fuente: https://earthquake.usgs.gov/data/comcat/

Contiene registros de eventos sísmicos. Entre las variables que pueden
aparecer:

  Variable      Descripción                Tipo
  ------------- -------------------------- ------------
  `time`        Fecha y hora               Fecha
  `latitude`    Latitud                    Numérica
  `longitude`   Longitud                   Numérica
  `depth`       Profundidad en km          Numérica
  `mag`         Magnitud                   Numérica
  `magType`     Tipo de magnitud           Categórica
  `place`       Descripción de ubicación   Texto
  `nst`         Número de estaciones       Numérica
  `gap`         Separación azimutal        Numérica
  `rms`         Error cuadrático medio     Numérica

**Posible tema:** identificación de patrones de eventos sísmicos según
magnitud y profundidad mediante GMM y EM.

Ventajas: - Gran cantidad de observaciones. - Variables físicas
cuantitativas. - Aplicación científica clara. - Posibilidad de
visualizar los componentes geográficamente.

Dificultades: - La distribución espacial está condicionada por placas
tectónicas y estructuras geológicas. - Puede haber secuencias temporales
y réplicas dependientes. - Magnitudes y profundidades pueden presentar
asimetrías y valores extremos. - Latitud y longitud requieren
tratamiento geográfico adecuado; no conviene aplicar sin más una
distancia euclidiana sobre coordenadas angulares.

Recomendación si se elige USGS: comenzar con magnitud y profundidad,
analizar distribuciones y ajustar GMM sin incorporar inicialmente
coordenadas geográficas. Después se puede mapear la localización de los
componentes como análisis complementario.

### 5.2. Inside Airbnb

Fuente: https://insideairbnb.com/es/get-the-data/

Proyecto que publica datos de alojamientos de Airbnb en diferentes
ciudades. Los archivos pueden incluir:

  Archivo                    Contenido                               Uso inicial
  -------------------------- --------------------------------------- -----------------
  `listings.csv.gz`          Información detallada de alojamientos   Sí
  `calendar.csv.gz`          Disponibilidad diaria                   No inicialmente
  `reviews.csv.gz`           Reseñas                                 No inicialmente
  `neighbourhoods.csv`       Barrios                                 Opcional
  `neighbourhoods.geojson`   Límites geográficos                     Opcional

Variables habituales de `listings`:

  Variable                   Descripción               Tipo
  -------------------------- ------------------------- ------------
  `price`                    Precio por noche          Numérica
  `latitude`                 Latitud                   Numérica
  `longitude`                Longitud                  Numérica
  `minimum_nights`           Mínimo de noches          Numérica
  `number_of_reviews`        Número de reseñas         Numérica
  `reviews_per_month`        Reseñas por mes           Numérica
  `availability_365`         Días disponibles al año   Numérica
  `accommodates`             Capacidad de huéspedes    Numérica
  `bedrooms`                 Habitaciones              Numérica
  `beds`                     Camas                     Numérica
  `room_type`                Tipo de alojamiento       Categórica
  `neighbourhood_cleansed`   Barrio                    Categórica

La disponibilidad exacta de variables depende del archivo y la ciudad.
Debe revisarse el diccionario de datos de la descarga concreta.

**Posible tema recomendado:** identificación de perfiles de alojamientos
turísticos en Bogotá mediante modelos de mezclas gaussianas estimados
con EM.

Ventajas: - Variables cuantitativas con interpretación económica y
operativa. - Resultados fáciles de comunicar. - Posibilidad de
caracterizar componentes mediante variables que no entrenaron el
modelo. - Posibilidad de mapear posteriormente los perfiles.

Dificultades: - Los precios suelen ser asimétricos y contener valores
extremos. - Puede haber datos faltantes. - Las variables tienen escalas
distintas. - Variables como `room_type` y `property_type` son
categóricas y no deben introducirse directamente en un GMM gaussiano
convencional. - La geografía puede introducir patrones espaciales que
requieren análisis específico.

### 5.3. Recomendación comparativa

La recomendación preliminar es **Inside Airbnb de Bogotá**. Permite
formular una pregunta de investigación sencilla de interpretar y
trabajar con variables cuantitativas de características de los
alojamientos.

No se debe afirmar antes del análisis que sus datos son más gaussianos o
que GMM necesariamente ajustará bien. Esto debe comprobarse
empíricamente.

------------------------------------------------------------------------

## 6. Propuesta de investigación con Inside Airbnb

### Título provisional

**Identificación de perfiles de alojamientos turísticos en Bogotá
mediante modelos de mezclas gaussianas y el algoritmo
Expectation-Maximization.**

### Pregunta de investigación

¿Es posible identificar perfiles diferenciados de alojamientos
turísticos en Bogotá a partir de sus características económicas,
capacidad y disponibilidad, utilizando modelos de mezclas gaussianas
estimados mediante Expectation-Maximization?

### Objetivo general

Aplicar un modelo de mezclas gaussianas, estimado mediante el algoritmo
Expectation-Maximization, para identificar y caracterizar perfiles de
alojamientos turísticos publicados en Airbnb en Bogotá.

### Objetivos específicos

1.  Realizar un análisis exploratorio de las características de los
    alojamientos.
2.  Identificar variables cuantitativas relevantes para aplicar GMM.
3.  Estimar diferentes modelos de mezclas gaussianas mediante EM.
4.  Determinar el número de componentes y la estructura de covarianza
    apropiados.
5.  Caracterizar los perfiles identificados mediante probabilidades de
    pertenencia.
6.  Evaluar la calidad y estabilidad de los resultados y comparar el
    modelo con K-Means.

### Variables iniciales propuestas

  Dimensión   Variable                  Nombre en la base
  ----------- ------------------------- --------------------
  Económica   Precio por noche          `price`
  Física      Capacidad de huéspedes    `accommodates`
  Física      Número de habitaciones    `bedrooms`
  Operativa   Días disponibles al año   `availability_365`
  Operativa   Estancia mínima           `minimum_nights`

Se propone transformar el precio:

\[ log_price=`\log`{=tex}(1+price) \]

y estandarizar las variables antes de ajustar el modelo.

No se recomienda incluir inicialmente latitud y longitud en el modelo si
el objetivo principal es caracterizar perfiles económicos y operativos.
Se pueden conservar para mapear los componentes posteriormente.

Tampoco se recomienda incluir directamente `room_type` en un GMM
gaussiano convencional. Puede utilizarse después para describir los
componentes.

------------------------------------------------------------------------

## 7. Procedimiento en RStudio

### 7.1. Instalar y cargar paquetes

``` r
install.packages(c(
  "tidyverse",
  "mclust",
  "janitor",
  "skimr",
  "cluster",
  "factoextra"
))

library(tidyverse)
library(mclust)
library(janitor)
library(skimr)
library(cluster)
library(factoextra)
```

Paquetes: - `tidyverse`: manipulación y visualización. - `mclust`:
modelos de mezclas gaussianas, EM y BIC. - `janitor`: limpieza de
nombres. - `skimr`: exploración descriptiva. - `cluster`: métricas de
agrupación. - `factoextra`: visualización de agrupaciones.

### 7.2. Descargar e importar los datos

Desde https://insideairbnb.com/es/get-the-data/:

1.  Buscar Bogotá.
2.  Descargar el archivo de alojamientos detallados (`listings.csv.gz` o
    equivalente).
3.  Guardar el archivo en una carpeta del proyecto.

Importación:

``` r
ruta <- file.choose()

airbnb <- readr::read_csv(
  ruta,
  show_col_types = FALSE
)

head(airbnb)
dim(airbnb)
names(airbnb)
glimpse(airbnb)
skimr::skim(airbnb)
```

`read_csv()` puede leer archivos `.csv.gz` directamente.

### 7.3. Limpiar nombres y seleccionar variables

``` r
airbnb <- airbnb %>%
  janitor::clean_names()

airbnb_modelo <- airbnb %>%
  select(
    id,
    price,
    accommodates,
    bedrooms,
    availability_365,
    minimum_nights,
    latitude,
    longitude,
    room_type,
    neighbourhood_cleansed
  )
```

Se conservan identificador, ubicación y variables categóricas para
interpretación posterior, aunque no entren todas al modelo.

### 7.4. Convertir precio y revisar datos

``` r
airbnb_modelo <- airbnb_modelo %>%
  mutate(
    price = readr::parse_number(price)
  )

glimpse(airbnb_modelo)

airbnb_modelo %>%
  summarise(
    across(everything(), ~sum(is.na(.)))
  )
```

Revisar valores extremos y rangos:

``` r
airbnb_modelo %>%
  summarise(
    precio_min = min(price, na.rm = TRUE),
    precio_max = max(price, na.rm = TRUE),
    capacidad_max = max(accommodates, na.rm = TRUE),
    habitaciones_max = max(bedrooms, na.rm = TRUE)
  )
```

No eliminar valores extremos automáticamente. Primero determinar si son
errores o valores válidos poco frecuentes.

### 7.5. Transformar el precio

``` r
airbnb_modelo <- airbnb_modelo %>%
  mutate(
    log_price = log1p(price)
  )
```

Visualizar:

``` r
ggplot(airbnb_modelo, aes(x = log_price)) +
  geom_histogram(
    bins = 40,
    fill = "steelblue",
    color = "white"
  ) +
  labs(
    title = "Distribución del precio transformado",
    x = "Logaritmo del precio",
    y = "Frecuencia"
  ) +
  theme_minimal()
```

La transformación puede reducir asimetría, pero no garantiza normalidad.

### 7.6. Construir una base analítica consistente

Es importante que la base auxiliar y la matriz del modelo conserven
exactamente las mismas filas y el mismo orden.

``` r
datos_analisis <- airbnb_modelo %>%
  mutate(
    log_price = log1p(price)
  ) %>%
  filter(
    is.finite(log_price),
    is.finite(accommodates),
    is.finite(bedrooms),
    is.finite(availability_365),
    is.finite(minimum_nights),
    !is.na(latitude),
    !is.na(longitude),
    accommodates >= 1,
    bedrooms >= 0,
    availability_365 >= 0,
    availability_365 <= 365,
    minimum_nights >= 1
  ) %>%
  drop_na(
    log_price,
    accommodates,
    bedrooms,
    availability_365,
    minimum_nights
  )
```

Construir matriz:

``` r
datos_gmm <- datos_analisis %>%
  select(
    log_price,
    accommodates,
    bedrooms,
    availability_365,
    minimum_nights
  ) %>%
  mutate(
    across(everything(), as.numeric)
  )

X <- scale(datos_gmm)
X <- as.matrix(X)

nrow(datos_analisis)
nrow(X)
```

Las dos cantidades deben coincidir.

Conviene revisar variables constantes o casi constantes y correlaciones
fuertes antes de estimar covarianzas completas.

### 7.7. Explorar distribuciones y relaciones

``` r
datos_gmm %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Valor"
  ) %>%
  ggplot(aes(x = Valor)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "white"
  ) +
  facet_wrap(
    ~Variable,
    scales = "free"
  ) +
  theme_minimal()
```

``` r
pairs(
  as.data.frame(X),
  pch = 19,
  cex = 0.3
)
```

Buscar asimetrías, valores extremos, relaciones entre variables y
patrones potenciales. La normalidad global no es requisito del GMM; se
evalúa la plausibilidad de una mezcla de componentes gaussianos.

### 7.8. Ajustar modelos GMM con `mclust`

``` r
set.seed(123)

modelo_gmm <- Mclust(
  X,
  G = 1:8
)

summary(modelo_gmm)
```

`Mclust()` estima modelos mediante EM y compara números de componentes y
estructuras de covarianza según los modelos disponibles.

Consultar resultados:

``` r
modelo_gmm$G
modelo_gmm$modelName
modelo_gmm$loglik
modelo_gmm$bic
```

En `mclust`, el BIC se expresa como:

\[ BIC=2`\ell`{=tex}(`\hat{\theta}`{=tex})-p`\log`{=tex}(n) \]

Por tanto, se prefiere el BIC más alto.

Visualización:

``` r
plot(
  modelo_gmm,
  what = "BIC"
)
```

### 7.9. Obtener probabilidades y clasificación

``` r
probabilidades <- modelo_gmm$z

head(probabilidades)

clasificacion <- modelo_gmm$classification

table(clasificacion)

incertidumbre <- modelo_gmm$uncertainty

summary(incertidumbre)
```

La incertidumbre individual es:

\[ U_i=1-`\max`{=tex}\_k P(Z_i=k`\mid `{=tex}x_i) \]

Una incertidumbre alta indica que la observación no tiene una
pertenencia claramente dominante.

### 7.10. Visualizar clasificación e incertidumbre

``` r
plot(
  modelo_gmm,
  what = "classification"
)

plot(
  modelo_gmm,
  what = "uncertainty"
)
```

### 7.11. Visualizar con PCA

PCA puede ayudar a representar en dos dimensiones los resultados de un
GMM ajustado sobre cinco variables. En este caso PCA es para
visualización, no necesariamente para entrenar el GMM.

``` r
pca <- prcomp(X)

datos_pca <- as.data.frame(
  pca$x[, 1:2]
)

datos_pca$Grupo <- factor(
  modelo_gmm$classification
)

ggplot(
  datos_pca,
  aes(
    x = PC1,
    y = PC2,
    color = Grupo
  )
) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Perfiles de alojamientos identificados mediante GMM",
    x = "Componente principal 1",
    y = "Componente principal 2"
  )
```

### 7.12. Caracterizar los componentes

Añadir clasificación a la base auxiliar, manteniendo correspondencia de
filas:

``` r
datos_analisis$Grupo <- factor(
  modelo_gmm$classification
)
```

Resumen de perfiles:

``` r
perfiles <- datos_analisis %>%
  group_by(Grupo) %>%
  summarise(
    Cantidad = n(),
    Precio_mediano = median(price, na.rm = TRUE),
    Capacidad_mediana = median(accommodates, na.rm = TRUE),
    Habitaciones_mediana = median(bedrooms, na.rm = TRUE),
    Disponibilidad_mediana = median(availability_365, na.rm = TRUE),
    Estancia_minima_mediana = median(minimum_nights, na.rm = TRUE),
    .groups = "drop"
  )

print(perfiles)
```

Los nombres de los perfiles se asignan después de examinar los
resultados. No asumir de antemano que habrá grupos económicos, medios y
costosos.

### 7.13. Describir los componentes con variables no utilizadas en el ajuste

Por ejemplo, comparar tipo de alojamiento:

``` r
datos_analisis %>%
  count(
    Grupo,
    room_type
  ) %>%
  group_by(Grupo) %>%
  mutate(
    Porcentaje = n / sum(n) * 100
  ) %>%
  ungroup()
```

También se pueden describir los componentes por barrio, sin interpretar
automáticamente una asociación como causal.

### 7.14. Evaluar con Silhouette

``` r
sil <- silhouette(
  modelo_gmm$classification,
  dist(X)
)

mean(sil[, 3])
```

Esta evaluación usa etiquetas duras y distancia euclidiana. Es
complementaria, no sustituye BIC ni la evaluación probabilística.

Incertidumbre promedio:

``` r
mean(modelo_gmm$uncertainty)
```

También se recomienda repetir el ajuste con diferentes semillas o
muestras para evaluar estabilidad.

------------------------------------------------------------------------

## 8. Comparación con K-Means

Se recomienda comparar ambos métodos usando las mismas variables y el
mismo número de grupos, tomando como referencia el número de componentes
elegido por GMM.

``` r
set.seed(123)

k <- modelo_gmm$G

modelo_kmeans <- kmeans(
  X,
  centers = k,
  nstart = 25
)

sil_gmm <- silhouette(
  modelo_gmm$classification,
  dist(X)
)

sil_kmeans <- silhouette(
  modelo_kmeans$cluster,
  dist(X)
)

comparacion <- data.frame(
  Modelo = c("GMM", "K-Means"),
  Silhouette = c(
    mean(sil_gmm[, 3]),
    mean(sil_kmeans[, 3])
  )
)

print(comparacion)
```

Comparar: 1. Estructura de grupos encontrada. 2. Separación según
Silhouette. 3. Diferencias en dispersiones y covarianzas. 4.
Probabilidades e incertidumbre disponibles en GMM. 5. Estabilidad de las
soluciones. 6. Interpretabilidad de los componentes.

No concluir que un algoritmo es superior únicamente por una métrica.

------------------------------------------------------------------------

## 9. Limitaciones y precauciones metodológicas

1.  **Inicialización:** EM puede alcanzar máximos locales. Se deben usar
    múltiples inicializaciones cuando la implementación lo permita y
    revisar estabilidad.
2.  **Número de componentes:** debe justificarse comparando modelos y
    considerando interpretación sustantiva.
3.  **Valores extremos:** pueden distorsionar medias y covarianzas. No
    eliminarlos sin justificación.
4.  **Covarianzas singulares:** pueden aparecer con pocos datos,
    variables redundantes o componentes pequeños. Considerar
    regularización o modelos más parsimoniosos.
5.  **Supuesto gaussiano:** la adecuación de la mezcla debe evaluarse,
    no darse por sentada.
6.  **Interpretación:** los componentes son patrones estadísticos, no
    necesariamente categorías naturales ni grupos causales.
7.  **Datos de Airbnb:** disponibilidad anual no equivale necesariamente
    a ocupación; precio publicado no equivale a precio efectivamente
    pagado; reseñas no equivalen al total de reservas.
8.  **Sesgo y cobertura:** los datos representan anuncios publicados y
    disponibles en la fuente, no todo el mercado de alojamiento
    turístico.
9.  **Causalidad:** el clustering identifica perfiles, no demuestra qué
    factores causan precios o disponibilidad.
10. **Datos temporales:** una descarga es una fotografía del mercado en
    una fecha determinada. Debe registrarse la fecha de descarga y la
    versión de los datos.

------------------------------------------------------------------------

## 10. Aplicación alternativa con USGS

Si se decide utilizar USGS, un tema posible es:

**Identificación de patrones de eventos sísmicos según magnitud y
profundidad mediante modelos de mezclas gaussianas.**

Variables iniciales: - `mag`: magnitud. - `depth`: profundidad.

Se recomienda empezar con estas variables y analizar sus distribuciones.
Las coordenadas geográficas pueden utilizarse para mapear los
componentes después del ajuste, pero no conviene introducirlas sin
tratamiento geográfico.

La distribución espacial y temporal de los terremotos puede presentar
dependencia, secuencias de réplicas y estructuras no gaussianas. Es
necesario considerar esas limitaciones antes de interpretar los
componentes como regiones o tipos de actividad sísmica.

Fuente y documentación: - https://earthquake.usgs.gov/data/comcat/ -
Servicio FDSN Event Web Service:
https://earthquake.usgs.gov/fdsnws/event/1/

------------------------------------------------------------------------

## 11. Investigaciones y bibliografía de referencia

Estas referencias permiten fundamentar el método y mostrar aplicaciones
en otras disciplinas:

1.  Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum
    likelihood from incomplete data via the EM algorithm. *Journal of
    the Royal Statistical Society: Series B (Methodological), 39*(1),
    1--22. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x

2.  Reynolds, D. A., & Rose, R. C. (1995). Robust text-independent
    speaker identification using Gaussian mixture speaker models. *IEEE
    Transactions on Speech and Audio Processing, 3*(1), 72--83.

3.  Reynolds, D. A., Quatieri, T. F., & Dunn, R. B. (2000). Speaker
    verification using adapted Gaussian mixture models. *Digital Signal
    Processing, 10*(1--3), 19--41.
    https://doi.org/10.1006/dspr.1999.0361

4.  Fraley, C., & Raftery, A. E. (2002). Model-based clustering,
    discriminant analysis, and density estimation. *Journal of the
    American Statistical Association, 97*(458), 611--631.
    https://doi.org/10.1198/016214502760047131

5.  McNicholas, P. D., & Murphy, T. B. (2010). Model-based clustering of
    microarray expression data via latent Gaussian mixture models.
    *Bioinformatics, 26*(21), 2705--2712.
    https://doi.org/10.1093/bioinformatics/btq498

6.  Liu, T.-C., et al. (2022). GMMchi: Gene expression clustering using
    Gaussian mixture modeling. *BMC Bioinformatics, 23*, 457.
    https://doi.org/10.1186/s12859-022-05006-0

### Documentación técnica

-   `mclust` en CRAN: https://cran.r-project.org/package=mclust
-   `mclust` y documentación: https://mclust-org.github.io/mclust/
-   Inside Airbnb: https://insideairbnb.com/es/get-the-data/
-   USGS ComCat: https://earthquake.usgs.gov/data/comcat/

------------------------------------------------------------------------

## 12. Estructura sugerida del informe académico

1.  **Introducción:** aprendizaje no supervisado y motivación del
    problema.
2.  **Marco teórico:** distribución normal multivariada, GMM, parámetros
    y EM.
3.  **Metodología:** fuente de datos, variables, limpieza,
    transformaciones y estandarización.
4.  **Estimación:** modelos candidatos, número de componentes,
    estructuras de covarianza y BIC.
5.  **Resultados:** probabilidades, clasificación, incertidumbre y
    perfiles.
6.  **Evaluación:** Silhouette como medida complementaria, estabilidad y
    comparación con K-Means.
7.  **Discusión:** interpretación, limitaciones y alcance.
8.  **Conclusiones:** respuesta a la pregunta de investigación.

### Productos esperados

-   Base analítica depurada y documentada.
-   Código reproducible en RStudio.
-   Análisis exploratorio.
-   Comparación de modelos mediante BIC.
-   Modelo GMM seleccionado.
-   Probabilidades posteriores e incertidumbre.
-   Visualizaciones de componentes.
-   Tabla de caracterización de perfiles.
-   Comparación con K-Means.
-   Conclusiones sustentadas en resultados reales.

------------------------------------------------------------------------

## 13. Decisiones pendientes antes de ejecutar el análisis

La recomendación actual es:

-   **Fuente:** Inside Airbnb.
-   **Ciudad:** Bogotá.
-   **Archivo:** Detailed Listings (`listings.csv.gz` o versión
    equivalente).
-   **Enfoque:** perfiles económicos, físicos y operativos.
-   **Variables iniciales:** precio transformado, capacidad,
    habitaciones, disponibilidad y estancia mínima.
-   **Lenguaje:** R.
-   **Paquete principal:** `mclust`.
-   **Método:** GMM estimado mediante EM.
-   **Selección:** BIC, comparando componentes y estructuras de
    covarianza.
-   **Evaluación:** incertidumbre, estabilidad, Silhouette
    complementario y comparación con K-Means.

Estas son decisiones preliminares. El número de componentes, la
estructura de covarianza y la interpretación de los perfiles solo se
definirán después de descargar, explorar y analizar los datos reales.
