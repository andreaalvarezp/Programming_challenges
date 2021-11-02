## ASSIGNMENT 2: Intensive integration using Web APIs
## Author: Andrea Álvarez Pérez

## Input file: co-expressed gene list we wanna know if are known to bind to one another.

## Objetivo: Encontrar todas las redes de interacción proteína-proteína que involucran
## a los miembros de esa lista de genes y determinar qué miembros de la lista de genes interactúan entre sí.

## UTILIZAR FILTROS EN EL CÓDIGO, POR ESPECIE, CALIDAD...

## BD a utilizar:
#   - BAR database from UToronto: devuelve coincidencias utilizando códigos AGI Locus.
# web: http://bar.utoronto.ca/interactions2/
# De alguna manera tengo que acceder a esa web por código, meter la lista de AGI y sacar la tabla con las interacciones.

#   - KEGG REST: anotar el KEGG Pathway del que forman parte los miembros, el KEGG ID y el Pathway name.
# web: https://www.genome.jp/kegg/pathway.html
# desde TOGO: 'http://togows.dbcls.jp/ENTRY O SEARCH CUAL ES LA DIFERENCIA/kegg-pathway' database = kegg-pathway, field: entry_id, pathways
## busquedas: ath (prefijo para arabidopsis thaliana)

# meter ahí los AGI Locus que forman parte de la red y sacar las rutas de las que forman parte (entry + name)

#   - Gene Ontology: anotar términos GO (biological_process) asociados al total de genes que forman la red, el GO:ID y el GO Term Name.
# web: http://geneontology.org/
# Puedo encontrar el GO term en togows uniprot/dr


# http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS
# web para acceder al REST API de cada BD

# 1. SACAR INTERACTION NETWORK DE BAR

require 'rest-client'

def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end

data = fetch('http://bar.utoronto.ca:9090/psicquic/webservices/current/search/query/tair:At2g13360?firstResult=0&format=tab27')

puts data.body

# De la info que me da BAR puedo coger el uniprotkb para buscar la gene ontology y el KEGG pathway desde TOGOWS
# PARA LA BUSQUEDA EN KEGG FILTRAR POR ESPECIE: codigo de arabidopsis "ath".

## DUDAS: COMO METER MAS DE UN GEN A LA VEZ A UN URL













