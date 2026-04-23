# Parcial II Bioinformática

## Punto 1: Calidad y Ensamblaje Genómico De Novo

### Análisis de calidad

El análisis de calidad se realizó con FastQC en el clúster:

```bash
module load fastqc

fastqc SRR8528336_1.fastq.gz SRR8528336_2.fastq.gz
```
<img width="647" height="357" alt="tabla1" src="https://github.com/user-attachments/assets/11eb0673-8907-41d1-bc30-79ed2ce34ecf" />
<img width="1097" height="818" alt="calidad1" src="https://github.com/user-attachments/assets/e7a02e20-ed0d-4833-be03-6decd9b4615f" />
<img width="629" height="352" alt="tabla2" src="https://github.com/user-attachments/assets/3af4f5b2-f54f-45ce-8b02-a739fbf04bae" />
<img width="1148" height="821" alt="calidad2" src="https://github.com/user-attachments/assets/2644c073-7d70-4d14-ac63-a65b6df41319" />

En general, las secuencias de ambos archivos (SRR8528336_1.fastq.gz y SRR8528336_2.fastq.gz) presentan buena calidad. Las lecturas tienen una longitud uniforme de 101 pb y no se observan secuencias clasificadas como de baja calidad. Además, el contenido de GC es de aproximadamente 43%, lo cual es un valor esperado.

Al revisar la calidad por posición, se observa que las lecturas mantienen valores altos en la mayor parte de su longitud, especialmente al inicio y en la región central. Hacia el extremo final (3’) hay una disminución gradual en la calidad, lo cual es típico en secuenciación Illumina. Esta caída es un poco más evidente en el segundo archivo, pero en general la calidad sigue siendo adecuada.



### Ensamblaje con SPAdes

El ensamblaje de novo se realizó con SPAdes:

```bash
module load spades/3.15.0

spades.py --careful \
-1 SRR8528336_1.fastq.gz \
-2 SRR8528336_2.fastq.gz \
-o spades_output
```



### Evaluación del ensamblaje

Se utilizaron los resultados de SPAdes para correr QUAST:

```bash
quast.py scaffolds.fasta -o quast_results
```
<img width="577" height="608" alt="report" src="https://github.com/user-attachments/assets/06dc3353-2b74-48af-b1f5-e3d8a67a4b68" />

En general, el ensamblaje no quedó muy continuo, sino bastante fragmentado. El N50 es de 1178 pb, lo que indica que la mitad del genoma está en fragmentos relativamente pequeños. Además, el L50 es de 507 contigs, lo que significa que se necesitan muchos fragmentos para cubrir el 50% del ensamblaje.

El tamaño total ensamblado fue de 3,089,161 pb, pero este valor está distribuido en muchos contigs pequeños. El contig más largo tiene una longitud de 15,601 pb, lo que también refleja que la continuidad del ensamblaje es baja.



## Punto 2: Búsqueda del Gen Candidato en el Genoma

### Preparación de secuencias

Se descargaron las secuencias proteicas desde NCBI y se concatenaron en un solo archivo:

<img width="862" height="865" alt="NCBI" src="https://github.com/user-attachments/assets/dac3e027-73c8-41de-9ab1-0378df7e81fc" />

```bash
cat *.fasta > proteins_query.fasta
```

Luego se editaron los encabezados para dejar únicamente el código de acceso y el nombre de la especie.

<img width="1411" height="450" alt="Atom" src="https://github.com/user-attachments/assets/1c3bda89-c506-47c8-ac24-36f7756c9f18" />


### Creación de la base de datos

Se creó una base de datos local con los contigs ensamblados:

```bash
makeblastdb -in contigs.fasta -dbtype nucl -parse_seqids -out db_contigs -title "Gasteracantha_genome"
```



### Ejecución de BLAST

Se realizó la búsqueda usando tblastn:

```bash
tblastn -query proteins_query.fasta \
-db db_contigs \
-outfmt 7 \
-out blast_obp \
-num_threads 1
```

<img width="1692" height="866" alt="blast" src="https://github.com/user-attachments/assets/db09f14a-8ab3-4145-8a5e-6ce2b48f43b8" />

### Interpretación de resultados

El contig más cercano a la proteína OBP es NODE_3138_length_349_cov_1.095238, ya que presenta el mejor resultado en el BLAST. Tiene un porcentaje de identidad de 88.1%, lo que indica una similitud bastante alta con la proteína, y un E-value de 7.55e-33, que es extremadamente bajo y confirma que el alineamiento es significativo y no se debe al azar. Además, el bit score es alto (172) y no hay gaps, lo que sugiere que el alineamiento es de buena calidad. Sin embargo, la longitud del alineamiento es de solo 59 aminoácidos, por lo que probablemente este contig corresponde solo a una parte del gen y no a la proteína completa.

## Punto 3: Análisis Poblacional y Gráficas en R

### Carga de datos

El análisis se realizó en R cargando las librerías necesarias:

```r
library(ggplot2)
library(dplyr)

df_altitud <- read.csv("df_altitud.csv", header=TRUE)
df_temperatura <- read.csv("df_temperatura.csv", header=TRUE)
```

Revisar
```r
head(df_altitud)
head(df_temperatura)
```


### Gráfica de abundancia relativa

```r
plot_altitud <- ggplot(df_altitud, aes(x=Altitud, y=Abundancia, fill=Grupo)) +
  geom_area() +
  scale_fill_manual(values = c("Blanco"="white",
                               "Amarillo"="yellow",
                               "Rojo"="red")) +
  ggtitle("Abundancia relativa de morfos según altitud") +
  xlab("Altitud") +
  ylab("Abundancia relativa")

pdf("grafica_altitud.pdf")
plot_altitud
dev.off()
```

<img width="717" height="717" alt="abundancia_relativa" src="https://github.com/user-attachments/assets/e912e271-ab11-4313-8ec2-d09f129337ce" />

### Gráfica de temperatura y probabilidad de supervivencia

```r
plot_temp <- ggplot(df_temperatura, aes(x=Temperatura, y=Supervivencia)) +
  geom_point() +
  ggtitle("Probabilidad de supervivencia según temperatura") +
  xlab("Temperatura") +
  ylab("Probabilidad de supervivencia") +
  theme(panel.background = element_rect(fill="white", colour="black"))

pdf("grafica_temperatura.pdf")
plot_temp
dev.off()
```

<img width="857" height="857" alt="temperatura" src="https://github.com/user-attachments/assets/0d214a0c-6b20-4ed6-804b-131a54d661d7" />

### Interpretación

A partir de las gráficas, se observa que la distribución de los morfos cambia a lo largo de la altitud, ya que el morfo rojo predomina en altitudes más bajas, mientras que el amarillo aumenta en altitudes intermedias y altas, y el blanco se encuentra en las zonas más altas. Esto sugiere que diferentes condiciones ambientales asociadas a la altitud favorecen distintos fenotipos. Por otro lado, la gráfica de temperatura muestra que la probabilidad de supervivencia tiende a aumentar con la temperatura, lo que indica una relación directa entre estas variables. En conjunto, estos resultados implican que tanto la altitud como la temperatura actúan como factores de selección, favoreciendo a los morfos mejor adaptados a ciertas condiciones, lo que refleja un proceso de adaptación local en la población.



## Punto 4: Entrega

Se comprimió la carpeta con todos los archivos generados en el clúster:
[zip](https://drive.google.com/file/d/1H0yIqhxXAs7BLx5WVvu_1SjLucq7hxFf/view?usp=sharing)

```bash
zip -r parcial2_gabrielarivas.zip parcial2_gabrielarivas/
```
