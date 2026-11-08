# Análisis de Expresión Diferencial y Enriquecimiento

## 1. Estudiante
**Álvaro Cintado Chaves** 

## 2. Lenguaje Elegido
**Python**

## 3. Instrucciones del Proyecto

### Creación del entorno 
Abrir una terminal en la carpeta raíz del proyecto y ejecutar los siguientes comandos:
(Se precisa de lenguaje ubuntu y miniconda3)

```bash
# Crear el entorno virtual .venv en la raíz del proyecto (pydeseq2) a partir de environment.yml: 
conda env create -p ./.venv -f environment.yml
```

### Ejecución del análisis
Previamente a ejecutar el análisis, es necesario que tanto los datos como los metadatos a utilizar se encuentren en la siguiente localización y con estos nombres específicos:

data/raw/counts.csv.gz
data/raw/metadata.csv.gz

```bash
#Otorgar permisos por primera vez para ejecutar ./main.sh:  
chmod +x main.sh
```
```bash
# Ejecutar desde la raíz del proyecto main.sh para realizar el análisis:
./main.sh 
```

### Localización de resultados
Una vez realizado el análisis mediante el script, los resultados se generarán en la careta results, tanto los
archivos tsv como los plots que aparecen a lo largo del código. En la carpeta logs aparecerán dos archivos: 
run_summary.txt con el resumen del análisis y session_info.txt con las versiones de python y librerías utilizadas
en este proyecto. 

## 4. Incidencias y decisiones
Los resultados pueden presentar variaciones respecto de un análisis similar llevado a cabo en lenguaje R,
ya que en el momento del enriquecimiento nos basamos en un filtro tal que (["Adjusted P-value"] < alpha = 0.05),
diferente de posibles umbrales cutoff de R. Además, también pueden existir diferencias debido a la lista de genes
que se usa en el enrich como términos GO y Reactome (en este proceso se utiliza GO_Biological_Process_2026 y 
Reactome_Pathways_2024).
 
## 5. Bibliografía
-Enlace del repositorio utilizado: https://github.com/alvarocintado19-dev/Bioinformatica-Analisis-de-expresion-y-enriquecimiento-Python

-Librerías utilizadas en el script:
PyDESeq2: https://pydeseq2.readthedocs.io
gseapy: https://gseapy.readthedocs.io
pandas: https://pandas.pydata.org/docs
matplotlib: https://matplotlib.org/stable