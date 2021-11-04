
# Bioinformatic programming challenges - Assignment2
# Author: Andrea Álvarez Pérez

## USAGE: $ruby Assignment2_main.rb ArabidopsisSubNetwork_GeneList.txt output_assignment2.txt

## Modules

require 'net/http'
require 'rest-client'
require './Gene.rb'
require './Interaction.rb'
require './Network.rb'
require 'json'
input = ARGV

## Function to fetch from URI

def fetch(url)
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s})
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

## Function to write records in a file:

def write_file(networks, out)
    # meter la info de GO y KEGG de la función Gene.annotate_data
    
    f = File.open(out, "a")
    networks.each do |id, net_feature|
      $MEM_GENES << net_feature.members.count
      f.puts "\r\nNETWORK #{id}\t\tNumber of nodes: #{net_feature.num_nodes}"
      f.puts "----------"
      
      net_feature.members.each do |gene_id, gene_obj|
        f.puts "GENE ID: #{gene_id}"
        i = 0
        while i < gene_obj.go_id.length
          f.puts "\tGO ID: #{gene_obj.go_id[i]}\tGO TERM: #{gene_obj.go_term[i]}"
          i += 1
        end
        j = 0
        while j < gene_obj.kegg_id.length
          f.puts "\tKEGG ID: #{gene_obj.kegg_id[j]}\tKEGG PATH: #{gene_obj.kegg_path[j]}"
          j += 1
        end
        f.puts "\n------------------------------------------------"
      end

    end
    
end

## Variables
$max_depth = 2

# Call Gene.rb with input file

## Error control

if input.length != 2
  puts "\r\n-----------------------------------------------------------------------------------"
  puts "\tERROR: Incorrect number of arguments"
  puts "\tUSAGE: $ruby Assignment2_main.rb <input file name> <output file name>"
  puts "-----------------------------------------------------------------------------------";puts
  exit 1
end

file = input[0]
out = input[1]
File.write(out, "\r\nASSIGNMENT 2\r\nAuthor: Andrea Álvarez Pérez\r\n------------------------------------------------\r\n", mode: "w") 

# 1. Load all the genes AGI Locus Codes in the list and create Interaction Objects with them, including IntAct code
Gene.load_file(file, out)
puts "Done!";puts

# 2. Once we have the list of all protein objects with IntAct ID associated in a ruby method, we assign a network ID

puts "Building networks and annotating genes..."

Interaction.return_method.each do |id, feature|
  
  if not feature.network                  # If the object is network = nil
    new_network = Network.new_net         # Create new network object with 2 nodes (the minimun possible)
    Network.build_network(feature, new_network)
  end
  
end

networks = Network.all_net

puts "$INT: #{$INT}"

puts "Protein objects: #{Interaction.return_method}"

# I want to store in an array the genes which are part of a network

$MEM_GENES = Array.new
puts "Writing into output file: output.txt..."
write_file(networks, out)

sum = 0
$MEM_GENES.each do |number|
  sum += number
end

$TOT_GENES = Gene.genes

puts "Done!"

puts "\r\n-----------------------------------------------------------------------------"
puts "\tOUTPUT SUMMARY:";puts
puts "--- Total number of genes loaded: #{$TOT_GENES.length}"
puts "--- Total number of networks: #{networks.length}"
puts "--- Number of genes of the initial list which are part of some network: #{sum}"
puts "\n-----------------------------------------------------------------------------";puts




