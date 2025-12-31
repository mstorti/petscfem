## -*- mode: cmake -*-

# 1. Definir los fuentes de la LIBRERÍA (Sin main.cpp)
# Agrega aquí todos los .cpp que conforman la lógica de advdif
set(ADVDIF_LIB_SOURCES
  advdif.cpp advdif_bubbly.cpp advdif_chim.cpp advdif_fsi.cpp		
  compdate.cpp advdifglb.cpp advdife.cpp advdifeo.cpp ffadvfm2_old.cpp	
  ffburg.cpp ffswfm2t.cpp wallswt.cpp ffadvfm2.cpp ffburg_old.cpp		
  enthalpy.cpp genload.cpp linhff.cpp diff.cpp aquifer.cpp stream.cpp		
  blessf.cpp bubbly_full.cpp bubbly_asm.cpp advec.cpp gasflow.cpp		
  advdfgth.cpp advhookf.cpp smoke.cpp ctff.cpp streamsw1d.cpp id.cpp		
  gaschem.cpp nonlres.cpp fm2funm.cpp bubbly_asm_mix.cpp lagmul.cpp		
  advabso.cpp gasflwgth.cpp svenant1d.cpp advpenal.cpp absolay.cpp		
  advdife-gcl.cpp advdife-ale.cpp streamsw2dt.cpp streamsw2d.cpp		
  strat_sw2d.cpp advdife_dual_time.cpp gasflow_preco.cpp streamsw2d.cpp	
  strat_sw2d.cpp advabsow.cpp advdif_dual_time.cpp bccadvfm2.cpp		
  bccadvfm2_old.cpp volabso.cpp advdife-bdf.cpp advdif-bdf.cpp		
  mmvforce.cpp chimera.cpp)

# 2. Definir los fuentes del EJECUTABLE (Solo el que tiene main())
set(ADVDIF_APP_SOURCES main.cpp)

# --- PASO A: CREAR LA LIBRERÍA (libadvdif.so) ---
add_library(advdif_lib ${ADVDIF_LIB_SOURCES})

# Configurar para que el archivo de salida se llame 'libadvdif.so' 
# (en lugar de libadvdif_lib.so)
set_target_properties(advdif_lib PROPERTIES OUTPUT_NAME "advdif")

# Vincular la librería local con el core
target_link_libraries(advdif_lib PUBLIC petscfem)

# --- PASO B: CREAR EL EJECUTABLE ---
add_executable(advdif_app ${ADVDIF_APP_SOURCES})

# Configurar el nombre final del ejecutable
set_target_properties(advdif_app PROPERTIES OUTPUT_NAME "advdif")

# Vincular el ejecutable con TU librería nueva (y hereda el core)
target_link_libraries(advdif_app PRIVATE advdif_lib)

file(GLOB ADVDIF_HEADERS "*.h")

install(TARGETS advdif_lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(TARGETS advdif_app
    RUNTIME DESTINATION bin
)

install(FILES ${ADVDIF_HEADERS}
    DESTINATION include/applications/advdif
)
