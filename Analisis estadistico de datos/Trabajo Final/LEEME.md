# Trabajo final: perfiles de alojamiento en Bogotá con GMM

El notebook corre de principio a fin con los datos de Bogotá. Ya están la preparación de los datos y el ajuste de los modelos (componentes principales, K-means, jerárquico y mezcla gaussiana), que es la parte más delicada. Faltan las gráficas, las tablas de resumen y lo que dicen los resultados. Lo repartí en tres bloques, y quién toma cada uno lo definimos en el grupo.

## Qué hay en esta carpeta

- `Trabajo Final.Rmd`: el notebook. Se abre en RStudio y se teje con Knit, que tarda unos dos minutos.
- `datos/airbnb_bogota.csv`: los 19.505 anuncios de la copia de Inside Airbnb del 14 de septiembre de 2026, con las 14 columnas que usamos y sin datos de los anfitriones. El bloque `crear-datos` del notebook muestra cómo se sacó del archivo original.
- `Presentación.pptx`: la presentación para el video, con la plantilla de la universidad. Tiene 29 diapositivas y la recortamos cuando sepamos cuánto puede durar el video.
- `logo_ur.png`: el logo del encabezado del PDF.
- `contexto_proyecto_GMM_EM_Inside_Airbnb_USGS.md`: la guía del proyecto, sin cambios.

Solo se necesitan dos paquetes, tidyverse y mclust. Si no tienen mclust, `install.packages("mclust")`.

## Qué falta

Todo lo pendiente dice `PENDIENTE`, así que basta con buscar esa palabra. En los bloques de código, el comentario dice qué gráfica o tabla hacer y con qué objeto del notebook. En el texto dice qué escribir, y debajo hay una línea con tres puntos que se reemplaza por lo que escribamos. Cuando algo esté hecho se borra su PENDIENTE.

En la presentación, cada diapositiva de resultados tiene un recuadro gris donde va la figura o la tabla y uno punteado que dice qué escribir.

| Quién | Secciones del notebook | Diapositivas |
|---|---|---|
| Bloque A | Datos (el tope de precio), Descripción y Componentes principales | 11 a 13 |
| Bloque B | K-means y Agrupamiento jerárquico | 14 a 16 |
| Bloque C | Modelo de mezcla gaussiana, Mezcla gaussiana en los datos de Bogotá y Perfiles | 17 a 22 |
| Bayron | Introducción, Comparación de los métodos, Conclusiones y revisión final | 1 a 10 y 23 a 29 |

Los colores de los grupos (`colores_grupos`) y el tema de las gráficas (`tema`) ya están en el bloque `colores`, para que todas las gráficas se vean iguales. La comparación de los métodos y las limitaciones ya tienen texto en el notebook. Las conclusiones las escribo cuando estén los tres bloques.

## Antes de cambiar algo

Lo que se puede cambiar está en el bloque `configuracion`: semilla 123, 4 grupos, tope de 5.000.000 COP por noche y muestra de 2.000 anuncios para el jerárquico. Si cambia la semilla, el número de grupos o el código de los modelos, cambia también la numeración de los grupos y hay que revisar lo que ya esté escrito sobre ellos.

Las diapositivas de datos y metodología ya tienen los números del notebook como está hoy. Si alguien cambia la configuración, me avisa para actualizarlas.

mclust se carga antes que tidyverse a propósito. Los dos tienen funciones `count()` y `map()`, y así quedan activas las de tidyverse.

Los resultados salen de R 4.6.1 y mclust 6.1.3. Con otras versiones pueden cambiar un poco.

## Cambios frente a la guía

La guía proponía precio, capacidad, habitaciones, días disponibles y noches mínimas. El número de habitaciones salió porque el 72,8 % de los anuncios tiene una sola y su correlación con la capacidad es 0,53, así que aporta poco. En su lugar entró el número de reseñas, como medida de la actividad del anuncio. Las noches mínimas y las reseñas entran en logaritmo, igual que el precio, porque tienen colas muy largas.

Muchos anuncios repiten exactamente el mismo valor (una noche mínima, treinta noches, dos personas). Con eso el modelo VVV no se puede estimar desde tres grupos, así que el GMM se ajusta con una distribución previa, `priorControl()`. Además el BIC sigue subiendo cuando se agregan grupos, así que los 4 grupos se escogen donde la curva se dobla y no en el máximo del BIC.
