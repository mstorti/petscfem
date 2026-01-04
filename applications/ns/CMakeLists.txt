## -*- mode: cmake -*-

# 1. Definir los fuentes de la LIBRERÍA (Sin main.cpp)
# Agrega aquí todos los .cpp que conforman la lógica de advdif
set(NS_LIB_SOURCES
  ns.cpp mmvmain.cpp mmv2main.cpp nsstruct.cpp ns_fsi.cpp compdate.cpp	
  nsglb.cpp nsitetlesfm2.cpp bccnsfm2.cpp wall.cpp nsvol.cpp walldata.cpp	
  nsilesther.cpp bccnstherfm2.cpp nsikeps.cpp nsikepsrot.cpp		
  wallke.cpp nonlr.cpp lagmul.cpp nssup.cpp nssupg.cpp nssupr.cpp nsid.cpp	
  nsres.cpp adaptor.cpp elast.cpp adappg.cpp qharm.cpp qharmm.cpp		
  nsirot.cpp blessf.cpp fracstep.cpp fstepfm2.cpp embgath.cpp		
  mmove.cpp fm2funm.cpp genload.cpp nsgath.cpp nshookf.cpp invcoupl.cpp	
  nullvort.cpp nsilesasm_2.cpp bccnsasmfm2.cpp bccfstep.cpp nslagmul.cpp	
  condwall.cpp condwallpen.cpp nspenal.cpp nsroblgml.cpp fstepfm2cw.cpp	
  nsitetlesc.cpp condwall-new.cpp mmove2.cpp nsitetasm.cpp		
  bccnsitetasm.cpp nsitetasm_avgvol.cpp ns_gasflow.cpp bccns_gasflow.cpp	
  mmoveopt.cpp elast2.cpp  bccnsitetasm_avgvol.cpp elastld.cpp		
  errestim.cpp mmoveopt3.cpp bubblyqint.cpp nsitetlesf.cpp nsitetlesfbf.cpp	
  poiboltz.cpp mmoveopt2.cpp truss.cpp nsitetlesfd.cpp nodeload.cpp		
  nsitetlesls.cpp electrophoresis.cpp  electrophoresisM.cpp flowrev.cpp	
  renorm.cpp renorm2.cpp renorm3.cpp nsitlfctx2.cpp					
  electrophoresisM2.cpp electrophoresisM.cpp electrophoresis_mov.cpp	
  poisson.cpp charge_cons.cpp pot_grad.cpp lubrication.cpp)

# 2. Definir los fuentes del EJECUTABLE (Solo el que tiene main())
set(NS_APP_SOURCES main.cpp)

# --- PASO A: CREAR LA LIBRERÍA (libadvdif.so) ---
add_library(ns_lib ${NS_LIB_SOURCES})

# Configurar para que el archivo de salida se llame 'libns.so' 
# (en lugar de libns_lib.so)
set_target_properties(ns_lib PROPERTIES OUTPUT_NAME "ns")

# Vincular la librería local con el core
target_link_libraries(ns_lib PUBLIC petscfem)

# --- PASO B: CREAR EL EJECUTABLE ---
add_executable(ns_app ${NS_APP_SOURCES})

# Configurar el nombre final del ejecutable
set_target_properties(ns_app PROPERTIES OUTPUT_NAME "ns")

# Vincular el ejecutable con TU librería nueva (y hereda el core)
target_link_libraries(ns_app PRIVATE ns_lib)

file(GLOB NS_HEADERS "*.h")

# Instalar la librería (libns.so)
install(TARGETS ns_lib
    EXPORT petscfemTargets       # <--- AGREGAR ESTA LÍNEA
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Instalar el ejecutable (opcional exportarlo, generalmente NO se exporta)
install(TARGETS ns_app
    RUNTIME DESTINATION bin
)

install(FILES ${NS_HEADERS}
    DESTINATION include/applications/ns)
