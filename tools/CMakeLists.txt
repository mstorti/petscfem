## -*- mode: cmake -*-
## This is petscfem/applications/advdif/CMakeLists.txt

# 2. Definir los fuentes del EJECUTABLE (Solo el que tiene main())
set(HEXASPLIT_APP_SOURCES hexasplit.cpp)

# --- PASO B: CREAR EL EJECUTABLE ---
add_executable(hexasplit_app ${HEXASPLIT_APP_SOURCES})

# Vincular el ejecutable con TU librería nueva (y hereda el core)
target_link_libraries(hexasplit_app PRIVATE project_lib PUBLIC petscfem ann_lib)

# Configurar el nombre final del ejecutable
set_target_properties(hexasplit_app PROPERTIES OUTPUT_NAME "hexasplit")

# Instalar el ejecutable (opcional exportarlo, generalmente NO se exporta)
install(TARGETS hexasplit_app RUNTIME DESTINATION bin)
