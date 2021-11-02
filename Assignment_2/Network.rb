
# Assignment2: Network.rb
# Author: Andrea Álvarez Pérez

## Este es para construir la network con los resultados anteriores de interacciones y escribirlo en el archivo output

require './Gene.rb' #Import of Gene object
require './Interaction.rb'

class Network
  
  attr_accessor :network_id # The number that identifies the network
  attr_accessor :num_nodes  # Number of nodes that the networks has
  attr_accessor :members    # Here we are going to store the gene OBJECTS of each network with all their properties
  
  @@num_net = 0             # Number of networks created
  @@network_obj = Hash.new
  @@gene_net = Array.new
  
  def initialize(params = {})
    
    @network_id = params.fetch(:network_id, "X")
    @num_nodes = params.fetch(:num_nodes, "0")
    @members = params.fetch(:members, Hash.new)
    
    @@num_net += 1                  # each time we initialize is for a new network
    @@network_obj[network_id] = self    # each time we initialize a new net[1], net[2]... is added where we save the instances
         
  end
  
  def self.all_net
    
    # This function retrieves the total of the Network objects
    return @@network_obj
    
  end
  
  def self.new_net
    
    # Network ID is a integer number starting in 1
    network_id = @@num_net + 1
    
    Network.new(
                :network_id => network_id,
                :num_nodes => 2,            # minimun number of nodes 
                :members => Hash.new
                )

    return network_id
      
  end
  
  def self.add_node(network_id)
    
    # 1. Add nodes to network
    @@network_obj[network_id].num_nodes += 1
    
  end
  
  
  def self.add_gene(gene_obj, network_id)
    
    # 2. Add genes and annotate them by filling the attributes KEGG and GO from the Gene objects
    #puts "Gene object GENE ID: #{gene_obj.gene_id}"
    @@network_obj[network_id].members[gene_obj.gene_id] = gene_obj
    
    gene_obj.annotate_data
    
    #puts @@network_obj
        
  end

  def self.build_network(feature, network_id)
      
    # 3. Walk through the arm of the network

    feature.network = network_id                   # Assign netwrok_id to the empty feature network_id of the interaction object
    
    # I have to iterate through all the interactions of that protein, that I saved it in @@interactions array
    $INT.each do |int|
      
      # Two possible cases: first intact_id of the $INT couple is the actual intact_id (int[0]), or the second one (int[1])
      # Also check if I have alredy a network with the other intact id
      if (int[0] == feature.intact_id) && (Interaction.return_method[int[1]].network == nil)
        self.add_node(network_id)
        #puts "Interaction.return_method[#{int[1]}], #{network_id}"     
        self.build_network(Interaction.return_method[int[1]], network_id) # Call the function itself (recursive)
      end
     
      if (int[1] == feature.intact_id) && (Interaction.return_method[int[0]].network == nil)
        self.add_node(network_id)
        #puts "Interaction.return_method[#{int[0]}], #{network_id}"
        self.build_network(Interaction.return_method[int[0]], network_id) # Call the function itself (recursive)
      end
   end
    
    # Once we have added all the proteins, now add the genes.
    # We search for the protein ID to return the identificator of the gene
    if Gene.genes[feature.prot_id]                          # If the protein has an associated gene, add it to the network
      #puts "Parameters: #{Gene.genes[feature.prot_id]}, #{network_id}"
      self.add_gene(Gene.genes[feature.prot_id], network_id)
    end
     
  end 

end