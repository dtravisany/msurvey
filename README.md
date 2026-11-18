# Receta para el análisis de Microbial Surveys en el tiempo

Práctico de amplicones 16S con QIIME 2: de las lecturas crudas a los cambios de la comunidad en el tiempo y el espacio

## Objetivos:

En el práctico de clases recorrimos con QIIME 2 el tutorial _Gut-to-soil axis_. Ahí los datos venían
importados, ya sin primers y con una profundidad bajísima. En este práctico cada grupo hará **el análisis completo
de un estudio real**, desde los archivos FASTQ tal como están en las bases de datos públicas, y responderá una
pregunta biológica sobre **cómo cambia una comunidad microbiana en el tiempo o a lo largo de un gradiente**.

Los tres estudios usan la misma región del gen 16S (**V4**, primers **515F/806R** del _Earth Microbiome
Project_) y el mismo secuenciador (**Illumina MiSeq 2×250**), así que los tres grupos corren el mismo flujo. Lo
que cambia es el ambiente, el diseño y la pregunta.

Se realizarán las siguientes tareas:

- Importación de lecturas paired-end a QIIME 2 con un archivo _manifest_.
- Revisión de la calidad de las lecturas con `qiime demux summarize`.
- (Grupo 1) Recorte de los primers con [__cutadapt__](https://cutadapt.readthedocs.io/).
- Corrección de errores y construcción de la tabla de ASVs con [__DADA2__](https://benjjneb.github.io/dada2/).
- Clasificación taxonómica con un clasificador _naive Bayes_ entrenado con [__SILVA 144__](https://www.arb-silva.de/) recortado a la región V4.
- Eliminación de cloroplastos y mitocondrias.
- Diversidad alfa y beta con [__q2-boots__](https://github.com/caporaso-lab/q2-boots) y [__q2-kmerizer__](https://github.com/bokulich-lab/q2-kmerizer).
- Análisis de series temporales con [__q2-longitudinal__](https://github.com/qiime2/q2-longitudinal).
- Pruebas estadísticas: Kruskal-Wallis, correlación de Spearman, PERMANOVA, PERMDISP y adonis.
- Abundancia diferencial con [__ANCOM-BC2__](https://github.com/FrederickHuangLin/ANCOMBC).
- Reproducibilidad con la procedencia de QIIME 2 (`replay-provenance`).


## Materiales:

Cada grupo analizará **todas las muestras** de un estudio publicado. Los datos son públicos y vienen del
[European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home); unos pocos archivos que en ENA están
dañados se obtuvieron desde el [SRA](https://www.ncbi.nlm.nih.gov/sra) de NCBI.

|Grupo| Tema | Estudio | Accesión | Muestras | Diseño |
| ----- | ----- | ----- | ----- | ----- | ----- |
|1| Océano (fiordo con salmonicultura) | Gutiérrez et al. 2018, _Frontiers in Marine Science_ 5:277. [doi:10.3389/fmars.2018.00277](https://doi.org/10.3389/fmars.2018.00277) | [`PRJNA359487`](https://www.ebi.ac.uk/ena/browser/view/PRJNA359487) | 50 | Fiordo Puyuhuapi (Aysén): una boya muestreada cada mes entre junio de 2013 y septiembre de 2014 a 2 y 20 m, más un perfil de 5 profundidades (1 a 100 m) en julio de 2015 |
|2| One Health | Lopatto et al. 2019, _PLoS ONE_ 14:e0220770. [doi:10.1371/journal.pone.0220770](https://doi.org/10.1371/journal.pone.0220770) | [`PRJNA485542`](https://www.ebi.ac.uk/ena/browser/view/PRJNA485542) | 130 | Purín de un plantel de 6.500 cerdos inyectado en un campo de maíz (Iowa): suelo antes, a los 3 días (en la línea de inyección y entre líneas) y a los 97 y 130 días, más el purín |
|4| Tratamiento de aguas | Brown et al. 2019, _Water Research_ 152:264–273. [doi:10.1016/j.watres.2019.01.003](https://doi.org/10.1016/j.watres.2019.01.003) | [`PRJEB48758`](https://www.ebi.ac.uk/ena/browser/view/PRJEB48758) | 198 | Planta de lodos activados (noreste de Inglaterra): lodo activado y agua servida decantada, muestreados cada semana durante dos años (2011–2013) |

Los datos ya se encuentran en el servidor y no es necesario descargarlos. Se ponen los links para que puedan
revisar de dónde vienen y descargar el artículo asociado a su estudio: **léanlo antes de empezar**, porque su
informe debe comparar sus resultados con los del artículo.

### Ubicación de los datos en el servidor:

```
/mnt/biostore/dipBG/microbial_surveys/
├── grupo1_puyuhuapi/        fastq/  manifest.tsv  metadata.tsv
├── grupo2_purin_suelo/      fastq/  manifest.tsv  metadata.tsv
├── grupo4_aguas_servidas/   fastq/  manifest.tsv  metadata.tsv
└── clasificador/            silva-144-ssu-nr99-515f-806r-classifier.qza
```

- `fastq/`: un par de archivos por muestra, `<muestra>_R1.fastq.gz` (forward) y `<muestra>_R2.fastq.gz` (reverse).
  Son las lecturas completas, tal como están en ENA.
- `manifest.tsv`: le dice a QIIME 2 qué archivos corresponden a cada muestra.
- `metadata.tsv`: la información de cada muestra (fecha, profundidad, tratamiento, variables ambientales). Ábranlo
  y revisen sus columnas antes de empezar; las van a necesitar en casi todos los pasos.

Columnas principales de cada `metadata.tsv`:

|Grupo| Columnas |
| ---- | ---- |
|1| `fecha`, `dias` (desde el 1 de junio de 2013), `mes`, `estacion` (austral), `serie` (`mensual_2013_2014` o `perfil_julio_2015`), `profundidad_m`, `capa` (`superficie` 1–5 m, `subsuperficie` 8–20 m, `profunda` 50–100 m), `nitrato_uM`, `fosfato_uM`, `silicato_uM`, `clorofila_a_ug_L`, `amonio_uM` (de la Tabla 1 del artículo) |
|2| `fecha`, `tipo_muestra` (`purin`, `suelo`, `suelo_linea_purin`), `grupo` (`purin`, `suelo_pre`, `linea_d3`, `suelo_d3`, `suelo_d97`, `suelo_d130`), `momento`, `dias_desde_aplicacion` (el purín se aplicó el 12 de noviembre de 2016), `replicado` |
|4| `fecha`, `semana`, `dias` (desde el 1 de junio de 2011), `anio_muestreo`, `estacion` (boreal), `tipo_muestra` (`lodo_activado` o `agua_servida_decantada`), y variables semanales de la planta: `temperatura_C`, `pH`, `oxigeno_disuelto_mg_L`, `caudal_afluente_L_s`, `TRH_h` (tiempo de retención hidráulica), `TRS_d` (edad del lodo), `DQO_afluente_mg_L`, `DQO_efluente_mg_L`, `NH4_afluente_mg_L`, `NH4_efluente_mg_L`, `nitrato_efluente_mg_L` y conteos de bacterias, bacterias oxidantes de amonio (AOB) y virus en el lodo |

### Clasificador taxonómico:

El clasificador del tutorial (GTDB r202) no tiene secuencias de cloroplastos ni de mitocondrias, así que no
permite eliminarlas, y en muestras de agua con fitoplancton los cloroplastos pueden ser varios por ciento de las
lecturas. Por eso en este práctico usaremos un clasificador entrenado con **SILVA 144** (versión NR99) y
recortado a la región que amplifican sus primers. Ya está entrenado en el servidor (tarda 18 minutos y usa 9 GB
de RAM). Así se construyó, con [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt), por si quieren hacer uno
para otra región:

```bash
qiime rescript get-silva-data --p-version 144 --p-target SSURef_NR99 --p-seq-format unaligned ...
qiime rescript reverse-transcribe ...                     # ARN -> ADN
qiime rescript cull-seqs ...                              # elimina secuencias con muchas bases ambiguas
qiime rescript filter-seqs-length-by-taxon --p-labels Archaea Bacteria Eukaryota --p-min-lens 900 1200 1400 ...
qiime rescript dereplicate --p-mode uniq ...
qiime feature-classifier extract-reads \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT ...                   # deja solo la V4, sin primers
qiime rescript filter-seqs-length --p-global-min 200 --p-global-max 450 ...
qiime rescript dereplicate --p-mode uniq ...
qiime feature-classifier fit-classifier-naive-bayes ...
```

`extract-reads` recuperó la región V4 del 99,5% de las referencias de bacterias, el 99,2% de las de arqueas, el
94% de los cloroplastos y el 88% de las mitocondrias.

### Software:

Todo el práctico se hace con **QIIME 2 2026.10**, la misma versión del tutorial de clases. Ya está instalado en
el servidor en el ambiente `rachis-qiime2-2026.10`.

### Activar el ambiente

Actívelo al inicio de cada sesión:

```bash
conda activate rachis-qiime2-2026.10
```

Sabrá que está activo porque el prompt cambia a `(rachis-qiime2-2026.10) usuario@servidor:$`. Recuerde que cada
comando `qiime` tarda unos 20 segundos en partir, aunque después haga muy poco.


# Inicio del Práctico
---

## 0. Carpeta de trabajo y variables

Cree una carpeta para el práctico y defina dos variables: `D`, con la carpeta de **su** grupo, y `C`, con el
clasificador.

```bash
mkdir -p ~/series_temporales
cd ~/series_temporales

D=/mnt/biostore/dipBG/microbial_surveys/grupo1_puyuhuapi        # grupo 1
# D=/mnt/biostore/dipBG/microbial_surveys/grupo2_purin_suelo    # grupo 2
# D=/mnt/biostore/dipBG/microbial_surveys/grupo4_aguas_servidas # grupo 4
C=/mnt/biostore/dipBG/microbial_surveys/clasificador/silva-144-ssu-nr99-515f-806r-classifier.qza

cp $D/metadata.tsv .
```

> Las variables se pierden si se desconecta. Cada vez que vuelva a entrar al servidor: `conda activate
> rachis-qiime2-2026.10`, `cd ~/series_temporales` y defina de nuevo `D` y `C`.

Algunos pasos tardan más de media hora. Córralos dentro de `screen`, para que sigan aunque se corte la conexión:

```bash
screen -S practico     # abre una sesión; adentro active el ambiente y defina D y C
# Ctrl-a d             # se sale de la sesión sin cerrarla
screen -r practico     # vuelve a la sesión
```

## 1. Los metadatos

```bash
qiime metadata tabulate \
  --m-input-file metadata.tsv \
  --o-visualization metadata.qzv
```

Como en el tutorial, los `.qzv` se copian al computador y se abren en [view.qiime2.org](https://view.qiime2.org):

```bash
scp <usuario>@<servidor>:series_temporales/metadata.qzv .
```

> Reemplace `<usuario>` por su usuario/grupo y `<servidor>` por la IP o el nombre del servidor del curso.

### En su informe deberá describir el diseño de su estudio con los metadatos: cuántas muestras hay de cada tipo, de qué fechas y en qué condiciones, y cuál es la variable de tiempo o de gradiente.

## 2. Importar las lecturas

En el tutorial los datos venían ya importados en `demux.qza`. Aquí los importamos nosotros a partir del
`manifest.tsv`, una tabla que asocia cada muestra con sus dos archivos FASTQ:

```bash
head -3 $D/manifest.tsv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path $D/manifest.tsv \
  --output-path demux.qza

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
```

| Parámetro | Descripción |
| ---- | ---- |
| `--type` | qué tipo de dato es: lecturas paired-end con calidad, separadas por muestra |
| `--input-format` | el formato de entrada: un _manifest_ con calidades Phred+33 |

La importación tarda alrededor de un minuto. Revise en `demux.qzv` cuántas lecturas tiene cada muestra (pestaña
_Overview_) y la calidad por posición (pestaña _Interactive Quality Plot_).

### (Solo grupo 1) Recortar los primers

Las lecturas del grupo 1 todavía traen los primers: las primeras 19 bases de `R1` son el 515F y las primeras 20
de `R2` son el 806R. Los primers se deben eliminar, porque son secuencia sintética: no informan sobre la bacteria y,
como son degenerados, pueden crear variantes falsas. Los grupos 2 y 4 se saltan este paso, porque sus lecturas
vienen sin primers.

```bash
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --p-cores 4 \
  --o-trimmed-sequences demux-sin-primers.qza \
  --o-stats cutadapt-stats.qza

qiime demux summarize \
  --i-data demux-sin-primers.qza \
  --o-visualization demux-sin-primers.qzv
```

| Parámetro | Descripción |
| ---- | ---- |
| `--p-front-f` / `--p-front-r` | primer que se busca al inicio de `R1` / `R2` (con las bases degeneradas del EMP: Y = C/T, M = A/C, N = cualquiera, V = A/C/G, W = A/T) |
| `--p-discard-untrimmed` | descarta los pares en los que no se encontró el primer |
| `--o-stats` | estadísticas por muestra de lo recortado |

Tarda un par de minutos. Compare `demux.qzv` con `demux-sin-primers.qzv`: el largo de las lecturas y cuántas se
perdieron.

### En su informe deberá mostrar el gráfico de calidad de sus lecturas y describir el número de lecturas por muestra (mínimo, mediana y máximo). El grupo 1 debe indicar cuántas lecturas se perdieron al recortar los primers y por qué es necesario recortarlos.

## 3. DADA2: de lecturas a ASVs

DADA2 corrige los errores de secuenciación y entrega la tabla de ASVs. Su decisión principal es dónde **truncar**
cada lectura: cortar donde la calidad cae, pero dejando suficiente largo para que `R1` y `R2` se solapen. El
amplicón V4 sin primers mide ~253 pb, y DADA2 exige al menos 12 pb de solapamiento, así que el largo truncado de
`R1` más el de `R2` debe superar los ~265 pb (con margen).

Corra DADA2 **dentro de `screen`**. Use los valores de su grupo:

|Grupo| `--i-demultiplexed-seqs` | `--p-trunc-len-f` | `--p-trunc-len-r` | Tiempo con 8 hilos |
| ---- | ---- | ---- | ---- | ---- |
|1| `demux-sin-primers.qza` | 220 | 180 | ~36 min |
|2| `demux.qza` | 240 | 200 | ~39 min |
|4| `demux.qza` | 180 | 140 | ~45 min |

> **Grupo 4:** los autores subieron las lecturas ya recortadas por calidad, así que no todas miden 250 pb (en
> `demux.qzv` se ve que el número de lecturas baja hacia el final). DADA2 descarta las lecturas más cortas que el
> largo de truncado, por eso aquí se trunca antes: con 180 y 140 se conserva cerca del 90% de las lecturas y
> todavía quedan ~67 pb de solapamiento.

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-sin-primers.qza \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 180 \
  --p-n-threads 8 \
  --o-representative-sequences asv-seqs.qza \
  --o-table asv-table.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza
```

> Por favor no use más de 8 hilos (`--p-n-threads`): el servidor lo comparten los tres grupos.

Revise cuántas lecturas pasaron cada etapa de DADA2:

```bash
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```

### En su informe deberá justificar los valores de truncado a partir del gráfico de calidad y del largo del amplicón, y reportar qué porcentaje de las lecturas sobrevivió a cada etapa (filtrado, corrección, unión de pares y quimeras). Comente qué muestras quedaron con pocas lecturas.

## 4. La tabla de ASVs y sus secuencias

```bash
qiime feature-table summarize \
  --i-table asv-table.qza \
  --m-metadata-file metadata.tsv \
  --o-summary asv-table.qzv \
  --o-sample-frequencies frecuencias-muestras.qza \
  --o-feature-frequencies frecuencias-asvs.qza

qiime feature-table tabulate-seqs \
  --i-data asv-seqs.qza \
  --m-metadata-file frecuencias-asvs.qza \
  --o-visualization asv-seqs.qzv
```

En `asv-seqs.qzv` revise el largo de las ASVs: la gran mayoría debería medir ~253 pb, el largo de la V4. Las que
miden mucho menos o mucho más suelen ser artefactos o secuencias que no son 16S bacteriano.

## 5. Taxonomía

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier $C \
  --i-reads asv-seqs.qza \
  --p-n-jobs 4 \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

Tarda uno o dos minutos. Cada ASV recibe un linaje (`d__` dominio, `p__` filo, `c__` clase, `o__` orden,
`f__` familia, `g__` género, `s__` especie) y una confianza.

> **Ojo con `Incertae_Sedis`.** SILVA usa ese nombre (en latín, «de posición incierta») para los grupos que no
> tienen un nombre formal en ese rango. Por eso pueden aparecer «géneros» llamados `g__Incertae_Sedis` que son
> linajes completamente distintos. Lea siempre el linaje completo: por ejemplo,
> `o__Chloroplast;f__Incertae_Sedis;g__Incertae_Sedis` es un cloroplasto, y `o__SAR324_clade;f__Incertae_Sedis` es
> una bacteria marina.

## 6. Eliminar cloroplastos, mitocondrias y ASVs raras

El 16S de los cloroplastos y de las mitocondrias también se amplifica, pero no son bacterias de la comunidad: son
organelos de algas, plantas o animales. Se eliminan junto con las ASVs que no llegaron a filo (casi siempre
artefactos). Después, como en el tutorial, se sacan las ASVs presentes en una sola muestra.

```bash
qiime taxa filter-table \
  --i-table asv-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table asv-table-sin-organelos.qza

qiime feature-table filter-features \
  --i-table asv-table-sin-organelos.qza \
  --p-min-samples 2 \
  --o-filtered-table asv-table-filtrada.qza

qiime feature-table filter-seqs \
  --i-data asv-seqs.qza \
  --i-table asv-table-filtrada.qza \
  --o-filtered-data asv-seqs-filtrada.qza

qiime feature-table summarize \
  --i-table asv-table-filtrada.qza \
  --m-metadata-file metadata.tsv \
  --o-summary asv-table-filtrada.qzv \
  --o-sample-frequencies frecuencias-muestras-filtrada.qza \
  --o-feature-frequencies frecuencias-asvs-filtrada.qza

qiime taxa barplot \
  --i-table asv-table-filtrada.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

En `taxa-bar-plots.qzv`, ordene las muestras por la variable de tiempo de su grupo (`dias` o
`dias_desde_aplicacion`) con el menú _Sort Samples By_, y mire los niveles 2 (filo) y 6 (género).

### En su informe deberá reportar qué porcentaje de las lecturas eran cloroplastos y mitocondrias, y mostrar el gráfico de barras ordenado en el tiempo (o a lo largo del gradiente), describiendo los filos y géneros dominantes y si tienen sentido para su ambiente.

## 7. La profundidad de muestreo

Como en el tutorial, antes de la diversidad hay que elegir una profundidad de rarefacción: las muestras con menos
lecturas quedan fuera.

```bash
qiime diversity alpha-rarefaction \
  --i-table asv-table-filtrada.qza \
  --p-max-depth 30000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```

Use la pestaña _Interactive Sample Detail_ de `asv-table-filtrada.qzv` para ver cuántas muestras pierde con cada
profundidad. Los valores que usamos en el ensayo fueron:

|Grupo| `--p-max-depth` | Profundidad elegida | Muestras que quedan |
| ---- | ---- | ---- | ---- |
|1| 30000 | 5000 | 50 de 50 |
|2| 20000 | 10000 | 126 de 129 |
|4| 40000 | 10000 | 196 de 197 |

### En su informe deberá mostrar la curva de rarefacción y justificar la profundidad elegida: cuántas muestras pierde y de qué grupos.

## 8. Diversidad alfa y beta con k-mers

Igual que en el tutorial, cambiando la profundidad por la de su grupo:

```bash
qiime boots kmer-diversity \
  --i-table asv-table-filtrada.qza \
  --i-sequences asv-seqs-filtrada.qza \
  --m-metadata-file metadata.tsv \
  --p-sampling-depth 5000 \
  --p-n 10 \
  --p-no-replacement \
  --output-dir boots-kmer-diversity
```

Tarda de 1 a 7 minutos según el grupo. Abra `boots-kmer-diversity/scatter_plot.qzv`. Además de la PCoA
(_braycurtis Axis 1_ contra _braycurtis Axis 2_, coloreando por su variable principal), aproveche que el gráfico
acepta columnas numéricas: ponga en el eje X la variable de tiempo de su grupo (`dias` o `dias_desde_aplicacion`)
y en el eje Y el _braycurtis Axis 1_ o la `shannon_entropy`. Así se ve la trayectoria de la comunidad en el tiempo.

### En su informe deberá mostrar la PCoA de Bray-Curtis y la de Jaccard, explicar la diferencia entre ambas métricas y describir qué variable de sus metadatos ordena mejor las muestras.

## 9. Cambios en el tiempo: volatilidad

`q2-longitudinal` está hecho para series temporales. El gráfico de volatilidad muestra cómo cambia una variable
(por ejemplo, la diversidad de Shannon) a lo largo del tiempo, con una línea por grupo.

|Grupo| `--p-state-column` | `--p-default-group-column` |
| ---- | ---- | ---- |
|1| `dias` | `capa` |
|2| `dias_desde_aplicacion` | `tipo_muestra` |
|4| `dias` | `tipo_muestra` |

```bash
qiime longitudinal volatility \
  --m-metadata-file metadata.tsv \
  --m-metadata-file boots-kmer-diversity/alpha_diversities/shannon.qza \
  --m-metadata-file boots-kmer-diversity/alpha_diversities/observed_features.qza \
  --p-state-column dias \
  --p-default-group-column capa \
  --p-default-metric shannon_entropy \
  --o-visualization volatilidad.qzv
```

En el menú del gráfico puede cambiar la métrica (`observed_features` o cualquier columna numérica de sus
metadatos, por ejemplo la temperatura) y la variable de agrupación.

### En su informe deberá mostrar el gráfico de volatilidad y describir cómo cambia la diversidad en el tiempo: ¿hay estacionalidad?, ¿vuelve a su estado inicial?

## 10. Pruebas estadísticas

Los gráficos sugieren diferencias; las pruebas dicen si son mayores de lo que daría el azar.

**Diversidad alfa.** Kruskal-Wallis compara la diversidad entre grupos de una variable categórica, y Spearman la
correlaciona con las variables numéricas (fecha, nutrientes, temperatura). Son iguales para todos los grupos:

```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity boots-kmer-diversity/alpha_diversities/shannon.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization shannon-grupos.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity boots-kmer-diversity/alpha_diversities/shannon.qza \
  --m-metadata-file metadata.tsv \
  --p-method spearman \
  --o-visualization shannon-correlacion.qzv
```

**Diversidad beta.** PERMANOVA prueba si la composición difiere entre grupos; PERMDISP, si los grupos difieren en
dispersión (su supuesto); y adonis entrega el R², la fracción de la variación que explica cada variable. Los
comandos dependen del diseño de cada grupo:

#### Grupo 1: superficie contra subsuperficie, y estaciones

La serie mensual y el perfil de 2015 son dos diseños distintos, así que se analizan por separado:

```bash
qiime diversity filter-distance-matrix \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[serie]='mensual_2013_2014'" \
  --o-filtered-distance-matrix bc-mensual.qza

qiime diversity beta-group-significance \
  --i-distance-matrix bc-mensual.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column capa \
  --p-pairwise \
  --o-visualization bc-mensual-capa-permanova.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix bc-mensual.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column capa \
  --p-method permdisp \
  --o-visualization bc-mensual-capa-permdisp.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix bc-mensual.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column estacion \
  --p-pairwise \
  --o-visualization bc-mensual-estacion-permanova.qzv

qiime diversity adonis \
  --i-distance-matrix bc-mensual.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "capa+estacion" \
  --o-visualization bc-mensual-adonis.qzv

qiime diversity filter-distance-matrix \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[serie]='perfil_julio_2015'" \
  --o-filtered-distance-matrix bc-perfil.qza

qiime diversity beta-group-significance \
  --i-distance-matrix bc-perfil.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column capa \
  --p-pairwise \
  --o-visualization bc-perfil-capa-permanova.qzv
```

#### Grupo 2: el suelo antes y después del purín

```bash
qiime diversity beta-group-significance \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column grupo \
  --p-pairwise \
  --o-visualization bc-grupo-permanova.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column grupo \
  --p-method permdisp \
  --o-visualization bc-grupo-permdisp.qzv

qiime diversity filter-distance-matrix \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[tipo_muestra]!='purin'" \
  --o-filtered-distance-matrix bc-suelo.qza

qiime diversity adonis \
  --i-distance-matrix bc-suelo.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "grupo" \
  --o-visualization bc-suelo-adonis.qzv
```

La PERMANOVA por pares del grupo 2 tarda unos 5 minutos (son 15 comparaciones).

#### Grupo 4: lodo contra afluente, y la estacionalidad del lodo

Primero se compara el lodo activado con el agua servida que entra al reactor; después, solo en el lodo, se
prueba el efecto de la estación y de las variables de operación de la planta:

```bash
qiime diversity beta-group-significance \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column tipo_muestra \
  --o-visualization bc-tipo-permanova.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column tipo_muestra \
  --p-method permdisp \
  --o-visualization bc-tipo-permdisp.qzv

qiime diversity filter-distance-matrix \
  --i-distance-matrix boots-kmer-diversity/distance_matrices/braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[tipo_muestra]='lodo_activado'" \
  --o-filtered-distance-matrix bc-lodo.qza

qiime diversity beta-group-significance \
  --i-distance-matrix bc-lodo.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column estacion \
  --p-pairwise \
  --o-visualization bc-lodo-estacion-permanova.qzv

qiime diversity adonis \
  --i-distance-matrix bc-lodo.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "temperatura_C+oxigeno_disuelto_mg_L+NH4_afluente_mg_L+anio_muestreo" \
  --o-visualization bc-lodo-adonis.qzv
```

> En adonis el orden de las variables importa: cada una explica la variación que dejaron las anteriores. Pruebe
> también otro orden y compare. No se incluye la edad del lodo (`TRS_d`) porque le faltan 13 semanas, y adonis no
> acepta valores faltantes.

### En su informe deberá presentar en una tabla los resultados de las pruebas (estadístico, p y q, y el R² de adonis), interpretarlos y discutir el resultado de PERMDISP: ¿la diferencia es de composición, de dispersión o de ambas?

## 11. Abundancia diferencial con ANCOM-BC2

Para saber **qué** géneros cambian, se agrupan las ASVs por género (nivel 6) y se usa ANCOM-BC2, como en el
tutorial.

#### Grupo 1: ¿qué géneros cambian con la profundidad y la estación?

```bash
qiime feature-table filter-samples \
  --i-table asv-table-filtrada.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[serie]='mensual_2013_2014'" \
  --o-filtered-table tabla-mensual.qza

qiime taxa collapse \
  --i-table tabla-mensual.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table tabla-mensual-genero.qza

qiime composition ancombc2 \
  --i-table tabla-mensual-genero.qza \
  --m-metadata-file metadata.tsv \
  --p-fixed-effects-formula "capa + estacion" \
  --p-reference-levels "capa::superficie" "estacion::invierno" \
  --o-ancombc2-output ancombc2-genero.qza

qiime composition da-barplot \
  --i-data ancombc2-genero.qza \
  --o-visualization ancombc2-genero.qzv
```

#### Grupo 2: ¿qué géneros del purín llegan al suelo y cuánto duran?

Se comparan las muestras de suelo con el suelo antes del purín (`suelo_pre`):

```bash
qiime feature-table filter-samples \
  --i-table asv-table-filtrada.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[tipo_muestra]!='purin'" \
  --o-filtered-table tabla-suelo.qza

qiime taxa collapse \
  --i-table tabla-suelo.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table tabla-suelo-genero.qza

qiime composition ancombc2 \
  --i-table tabla-suelo-genero.qza \
  --m-metadata-file metadata.tsv \
  --p-fixed-effects-formula grupo \
  --p-reference-levels "grupo::suelo_pre" \
  --o-ancombc2-output ancombc2-genero.qza

qiime composition da-barplot \
  --i-data ancombc2-genero.qza \
  --o-visualization ancombc2-genero.qzv
```

#### Grupo 4: ¿qué géneros del lodo cambian con las estaciones?

```bash
qiime feature-table filter-samples \
  --i-table asv-table-filtrada.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[tipo_muestra]='lodo_activado'" \
  --o-filtered-table tabla-lodo.qza

qiime taxa collapse \
  --i-table tabla-lodo.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table tabla-lodo-genero.qza

qiime composition ancombc2 \
  --i-table tabla-lodo-genero.qza \
  --m-metadata-file metadata.tsv \
  --p-fixed-effects-formula estacion \
  --p-reference-levels "estacion::invierno" \
  --o-ancombc2-output ancombc2-genero.qza

qiime composition da-barplot \
  --i-data ancombc2-genero.qza \
  --o-visualization ancombc2-genero.qzv
```

ANCOM-BC2 tarda entre 4 y 6 minutos. En `ancombc2-genero.qzv` revise, para cada comparación, el LFC (en
logaritmo natural: un LFC de 1 son 2,7 veces), el q y la columna `passed_ss`: si es falsa, el resultado cambia
según el pseudoconteo que se suma a los ceros, y es frágil.

### En su informe deberá presentar una tabla con los géneros significativos (q < 0,05) de cada comparación, con su LFC, su q y si pasan el análisis de sensibilidad, e interpretar su biología. Si ningún resultado es robusto, discuta por qué y qué haría para confirmarlo.

## 12. Reproducibilidad

Como en el tutorial, QIIME 2 escribe solo el script con todo lo que hicieron:

```bash
cd ~
qiime tools replay-provenance \
  --in-fp series_temporales/ \
  --recurse \
  --out-fp series-temporales-replayed.bash
```


# Resumen del flujo

```
FASTQ (R1/R2) + manifest.tsv
   │  qiime tools import  →  demux.qza  →  demux summarize (calidad)
   │  (grupo 1: cutadapt, recorte de primers)
   ▼
DADA2  →  tabla de ASVs + secuencias
   │  classify-sklearn (SILVA 144 V4)  →  ¿QUIÉNES están?
   │  filter-table (sin cloroplastos ni mitocondrias) + filter-features
   ▼
alpha-rarefaction  →  profundidad
   │
   ├──► kmer-diversity  →  alfa (Shannon, riqueza) y beta (Bray-Curtis, Jaccard, PCoA)
   │        ├──► volatility  →  ¿CÓMO cambia en el tiempo?
   │        └──► Kruskal-Wallis, Spearman, PERMANOVA, PERMDISP, adonis  →  ¿es real la diferencia?
   │
   └──► taxa collapse + ANCOM-BC2  →  ¿QUÉ géneros cambian?
```

# Entrega del informe

Su informe debe incluir, para su estudio:

1. **Diseño**: descripción de las muestras y de la variable de tiempo o gradiente, a partir de los metadatos.
2. **Calidad y DADA2**: lecturas por muestra, gráfico de calidad, valores de truncado justificados y porcentaje
   de lecturas que sobrevivió a cada etapa (y, en el grupo 1, el recorte de primers).
3. **Composición taxonómica**: porcentaje de cloroplastos y mitocondrias, y gráfico de barras ordenado en el
   tiempo o el gradiente, con los taxones dominantes.
4. **Diversidad**: curva de rarefacción y profundidad elegida, PCoA de Bray-Curtis y de Jaccard, y gráfico de
   volatilidad.
5. **Estadística**: tabla con las pruebas de diversidad alfa y beta (Kruskal-Wallis, Spearman, PERMANOVA,
   PERMDISP y adonis) y su interpretación.
6. **Abundancia diferencial**: tabla de géneros significativos de ANCOM-BC2 con su LFC, su q y `passed_ss`.
7. **Discusión**: ¿responden sus resultados la pregunta del estudio? Compárelos con los del artículo.
