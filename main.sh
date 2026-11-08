# Detener el script en caso de error
set -e
echo "1. Activando el entorno..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ./.venv
echo "2. Ejecutando script..."
python src/01_run_deseq2.py
echo "Ejecutado con éxito"