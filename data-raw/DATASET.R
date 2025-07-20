usethis::use_data_raw()

usethis::use_data(DATASET, overwrite = TRUE)

# --- CÓDIGO PARA PREPARAR Y GUARDAR TUS DATASETS ---

# --- CÓDIGO PARA PREPARAR Y GUARDAR TUS DATASETS ESPECÍFICOS ---

# 1. Cargar las librerías necesarias
library(readxl)
library(stringr) # Para limpiar los nombres de los ficheros

# 2. Definir las rutas a las carpetas donde están tus datos
ruta_viejos2 <- "C:/Users/Usuario/Documents/toledo23/paper beneficio v2/viejos/viejos2"
ruta_viejos <- "C:/Users/Usuario/Documents/toledo23/paper beneficio v2/viejos"

# 3. Crear la lista EXACTA de ficheros a cargar
nombres_ficheros_viejos2 <- c(
  "AURA3_1A.xlsx", "BOLERO3_2.xlsx", "CONCUR_2A.xlsx", "COUAA202_2.xlsx",
  "IMELDA_2A.xlsx", "KEYNOTE024_1A.xlsx", "LUXLung8_2A.xlsx", "MASTER.DATA.xlsx",
  "MASTER.DATA2.xlsx", "MONALEESA2_1.xlsx", "PALOMA2_1A.xlsx", "PROFILE1014_1A.xlsx",
  "RECOURSE_1A.xlsx", "toga.xlsx"
)

nombres_ficheros_viejos <- c(
  "NSABPB40_3B.xlsx", "NSABPB40_3A.xlsx", "NSABPB35_3.xlsx", "NSABPB35_2.xlsx",
  "NOAH_2B.xlsx", "NOAH_2A.xlsx", "NO16968_2C.xlsx", "NO16968_2A.xlsx",
  "MA17R_1A.xlsx", "GETUG12_2.xlsx", "GEICAM2003_2A.xlsx", "DART0105_2A.xlsx"
)

# Combinar las rutas con los nombres para obtener las rutas completas
ficheros_a_cargar <- c(
  file.path(ruta_viejos2, nombres_ficheros_viejos2),
  file.path(ruta_viejos, nombres_ficheros_viejos)
)

# 4. Bucle para leer cada fichero y asignarlo a un objeto en R
for (fichero in ficheros_a_cargar) {
  # Crear un nombre de objeto válido en R
  # "MASTER.DATA.xlsx" se convertirá en "MASTER_DATA" para ser un nombre válido
  nombre_objeto <- basename(fichero) |>
    str_remove("\\.xlsx$") |>
    str_replace_all("\\.", "_") # Reemplaza puntos por guiones bajos

  # Leer el fichero Excel y asignarlo al nuevo nombre de objeto
  assign(nombre_objeto, read_excel(fichero))

  cat("Cargado:", basename(fichero), "como ->", nombre_objeto, "\n")
}

# 5. Guardar TODOS los objetos creados en la carpeta /data del paquete
nombres_objetos_creados <- basename(ficheros_a_cargar) |>
  str_remove("\\.xlsx$") |>
  str_replace_all("\\.", "_")

# Usamos usethis::use_data() para guardarlos todos a la vez
usethis::use_data(list = nombres_objetos_creados, overwrite = TRUE)

cat("\n¡Proceso completado! Los", length(ficheros_a_cargar), "datasets especificados han sido guardados.\n")
