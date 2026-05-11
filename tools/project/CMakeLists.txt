## -*- mode: cmake -*-
## This is petscfem/applications/advdif/CMakeLists.txt

# 1. Definir los fuentes de la LIBRERÍA (Sin main.cpp)
# Agrega aquí todos los .cpp que conforman la lógica de advdif
set(PROJECT_LIB_SOURCES project.cpp)

# 2. Definir los fuentes del EJECUTABLE (Solo el que tiene main())
set(PROJECT_APP_SOURCES project4.cpp)

# --- PASO A: CREAR LA LIBRERÍA (libproject.so) ---
add_library(project_lib ${PROJECT_LIB_SOURCES})

# Configurar para que el archivo de salida se llame 'libproject.so' 
# (en lugar de libproject_lib.so)
set_target_properties(project_lib PROPERTIES OUTPUT_NAME "project")

# Vincular la librería local con el core
target_link_libraries(project_lib PUBLIC petscfem ann_lib)

# --- PASO B: CREAR EL EJECUTABLE ---
add_executable(project_app ${PROJECT_APP_SOURCES})

# Configurar el nombre final del ejecutable
set_target_properties(project_app PROPERTIES OUTPUT_NAME "project")

# Vincular el ejecutable con TU librería nueva (y hereda el core)
target_link_libraries(project_app PRIVATE project_lib)

file(GLOB PROJECT_HEADERS "*.h")

# Instalar la librería (libproject.so)
install(TARGETS project_lib
    EXPORT petscfemTargets       # <--- AGREGAR ESTA LÍNEA
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Instalar el ejecutable (opcional exportarlo, generalmente NO se exporta)
install(TARGETS project_app
    RUNTIME DESTINATION bin
)

install(FILES ${PROJECT_HEADERS}
    DESTINATION include/applications/project
)
