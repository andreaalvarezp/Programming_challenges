
# Assignment2: Interaction_protein.rb
# Author: Andrea Álvarez Pérez


require 'json'
require 'net/http'
require './Gene.rb'

class Interaction
    
    attr_accessor :prot_id        # UniProt ID
    attr_accessor :intact_id      # IntAct ID
    attr_accessor :network        # Network ID
    
    @@prot_intact = Hash.new      # Not all the protein objects are going to have an IntAct ID, we store here the ones with it
    @@interactions = false        # Counts the interactions
    
    def initialize(params = {})
      
      @prot_id = params.fetch(:prot_id, "XXXXXX") 
      @intact_id = params.fetch(:intact_id, nil)
      @network = params.fetch(:network, nil)

      # Only add new object of there is an intact ID present
      if intact_id
        @@prot_intact[intact_id] = self
      end
      
      if @@interactions == false
          @@interactions = Array.new        # Array containing all interactions, we want it to be acumulative
          # we only start this method once, and store all interactions of all genes
      end

      
    end
    
    def self.return_method
        
        # Only to use the methods in the main script. The idea is to save in thi method all the protein object created that have an IntAct ID
        # in order to assign later a network ID to configure
        
        return @@prot_intact
        
    end
    
    def self.get_prot(prot_id, depth, gene_id = nil, intact_id = nil)
      
      # I need a function to create a protein object in order to collect the information for the network and annotation
      # If it i the first search, we are going to have the gene_id, but if it's the second search, we are going to have de intact_id
      
      if not intact_id                          # If the gene_id i given, I search for Intact_id
        intact_id = self.get_intact_id(gene_id)
      end
      
      if intact_id && (depth.to_f < $max_depth)                     # If the IntAct is given, I search for interactions
        new_interaction(intact_id, depth)        # Llamar a get_prot desde otro lado que no sea load_file (siempre va a ser un gene_id)                 
      end

      #puts "@@interactions: #{@@interactions}"
      Interaction.new(
            :prot_id => prot_id,
            :intact_id => intact_id,
            :network => nil 
            )
      #puts "Protein object created: #{prot_id}, #{intact_id}"   
      
      depth += 1 # Go ahead with the network depth
 
    end
    
    def self.prot_exist(intact_id)

        if @@prot_intact.has_key?(intact_id)
          return true
        
        else
          return false
        
        end
          
    end
    
    def self.get_intact_id(gene_id)
        
        
        # Get intact_id from gene_id
        response = fetch("http://togows.org/entry/ebi-uniprot/#{gene_id}/dr.json")
        
        if response
  
          data = JSON.parse(response)
      
          if data[0]['IntAct']
            #puts "IntAct ID: #{data[0]['IntAct'][0][0]}" # If the protein is present
            return data[0]['IntAct'][0][0]
          else
            return nil                     # Empty result
          end
        else
            puts "Web call failed in function get_intact_id - see STDERR for details..."
        end
      
    end
    
    def self.convert_intact_protid(intact_id)
        
        #Get prot_id from intact_id
        response = fetch("http://togows.org/entry/ebi-uniprot/#{intact_id}/entry_id.json")
        
        if response

            data = JSON.parse(response)
            prot_id = data[0]
            
            return prot_id
        else
            puts "Web call failed in convert_intact_protid - see STDERR for details..."
        end
  
    end
    
    def self.new_interaction(intact_id, depth)
        
        # I have to create the proteins with the elements of the network I got in get_interactions
        #puts "@@interactions: #{@@interactions}"
        interactions = get_interactions(intact_id)
        depth += 1
        
        if interactions                 # I get the rest of the interactions
            #puts "Adding protein-protein interactions..."
            intact = intact_id          # I don't want my IntAct ID to change when I create a new protein object
            
            interactions.each do |prot1, prot2|
                # Check if the added protein is the first or the second in the array and if it doesn't exist as an object
                
                if prot1 == intact && (not prot_exist(prot2))
                    prot2_id = convert_intact_protid(prot2)
                    # Create new protein object
                    self.get_prot(prot2_id, depth, gene_id = nil, intact_id = prot2)    # search for interactions at 2º depth level
                    
                elsif prot2 == intact && (not prot_exist(prot1))    
                    prot1_id = convert_intact_protid(prot1)
                    # Create new protein object
                    self.get_prot(prot1_id, depth, gene_id = nil, intact_id = prot1)     # search for interactions at 2º depth level
                    
                end
            end

        end
        #puts "Done!"

    end
  
    def self.get_interactions(intact_id)
        
        # This function retrieves the IntAct URL with the interactions of my protein.
        
        puts "Fetching for interactions for IntAct ID #{intact_id}..."
        response = fetch ("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{intact_id}?format=tab25")
        #puts "@@interactions: #{@@interactions}"
        
        if response

            data = response.to_s.split("\n")
            
            # create an array to save the protein interactions by couples
            interactions = Array.new
          
            data.each do |line|
                          
                uni = line.split("\t")
               
                uniprot1 = uni[0].split(":")[1]                     # First protein 
                uniprot2 = uni[1].split(":")[1]                     # Second protein
                # Ways of measuring protein-protein interaction: the least effective is two hybrid, so in order to reduce the hits, we filter through this
                hybrid = uni[6].split("(")[1]
                score = uni[14].sub(/intact-miscore:/, "").to_f     # Quality score
    
                ## FILTERS
                
                if hybrid == "two hybrid)" || hybrid == "two hybrid array)"    # Measure protein-protein interaction
                    next
                elsif score < 0                                                     # Quality score > 0
                    next
                elsif uniprot1 == uniprot2
                    next
                elsif uni[9] =~ /taxid:3702/ && uni[10] =~ /taxid:3702/              # Species: Arabidopsis
                    if uniprot1 < uniprot2                                      # Alphabetical order
                        if not @@interactions.include?([uniprot1, uniprot2])    # Avoid repetitions
                            @@interactions << [uniprot1, uniprot2]
                            interactions << [uniprot1, uniprot2]
                        end
                    else                                   # No alphabetical order
                        if not @@interactions.include?([uniprot2, uniprot1])    # Avoid repetitions
                            @@interactions << [uniprot2, uniprot1]
                            interactions << [uniprot2, uniprot1]
                        end
                    end
                end
            end
            
            # List in global variables all the interactions and proteins with IntAct ID
            #puts "Done!"
            $INT = @@interactions
           
            return interactions        # proteins which belong to the same network
        else
            puts "Web call failed in function get_interactions - see STDERR for details..."
        end

    end
  
end


